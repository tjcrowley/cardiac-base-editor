"""
Structured, consent-gated query functions over an ingested subject's genomic
data. This is the one place both the CLI (`cbe query`), the web UI's Query
tab, and the NL front door (query/nl.py) all call into — no code path
bypasses subjects.require_consent().

Every function takes an explicit vcf_path, same convention as
genomic_intake.cli's `run` command and the web UI's run route — this module
doesn't invent a new "the subject's VCF" registry; callers point at whichever
file under that subject's genomic_intake.storage.raw_dir() they mean.
"""

from __future__ import annotations

from dataclasses import asdict

from cardiac_base_editor.genomic_intake import subjects
from cardiac_base_editor.genomic_intake.extract import (
    _fetch_cds_genomic_map,
    build_personalized_cds,
    parse_vcf_snvs,
)
from cardiac_base_editor.cancer.mhc_binding import rank_binders
from cardiac_base_editor.cancer.neoantigen import candidate_peptides
from cardiac_base_editor.models.off_target import score_off_target
from cardiac_base_editor.models.protein_consequence import score_substitution
from cardiac_base_editor.pipeline import (
    CODON_TABLE,
    KNOWN_TARGETS,
    fetch_cds,
    run as pipeline_run,
)


def _resolve_transcript(gene: str) -> str:
    return KNOWN_TARGETS.get(gene.upper(), gene)


def _variant_to_protein_change(transcript_id: str, vcf_path: str, genomic_pos: int) -> dict:
    """
    Shared by explain_variant and rank_neoantigens: resolves a genomic SNV to
    its codon-level consequence and the full translated protein sequence.
    """
    cds = fetch_cds(transcript_id)
    genomic_to_cds, strand = _fetch_cds_genomic_map(transcript_id, len(cds))

    variant = next((v for v in parse_vcf_snvs(vcf_path) if v.pos == genomic_pos), None)
    if variant is None:
        raise ValueError(f"No SNV at genomic position {genomic_pos} in {vcf_path}")

    cds_pos = genomic_to_cds.get(genomic_pos)
    if cds_pos is None:
        raise ValueError(f"Genomic position {genomic_pos} is not within {transcript_id}'s CDS")

    alt_base = variant.alt.translate(str.maketrans("ACGT", "TGCA")) if strand == -1 else variant.alt

    codon_idx = (cds_pos - 1) // 3
    base_in_codon = (cds_pos - 1) % 3
    ref_codon = cds[codon_idx * 3: codon_idx * 3 + 3]
    alt_codon = ref_codon[:base_in_codon] + alt_base + ref_codon[base_in_codon + 1:]
    ref_aa = CODON_TABLE.get(ref_codon, "?")
    alt_aa = CODON_TABLE.get(alt_codon, "?")

    full_translation = [CODON_TABLE.get(cds[i:i + 3], "?") for i in range(0, len(cds) - len(cds) % 3, 3)]
    protein_seq = "".join(full_translation).rstrip("*")  # ESM-2's vocabulary has no stop-codon token

    return {
        "variant": variant,
        "codon_idx": codon_idx,
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "protein_seq": protein_seq,
    }


def list_variants(subject_id: str, gene: str, vcf_path: str, operator: str = "query-engine") -> list[dict]:
    """
    All of the subject's SNVs (from vcf_path) that fall within the target
    gene's coding sequence, with their CDS-relative position.
    """
    subjects.require_consent(subject_id, gene, operator=operator)

    transcript_id = _resolve_transcript(gene)
    cds = fetch_cds(transcript_id)
    genomic_to_cds, strand = _fetch_cds_genomic_map(transcript_id, len(cds))

    hits = []
    for v in parse_vcf_snvs(vcf_path):
        cds_pos = genomic_to_cds.get(v.pos)
        if cds_pos is None:
            continue
        hits.append({
            "chrom": v.chrom, "genomic_pos": v.pos, "ref": v.ref, "alt": v.alt,
            "cds_position": cds_pos, "gene": gene, "transcript_id": transcript_id,
        })
    return hits


def rank_guides(subject_id: str, gene: str, vcf_path: str, editor: str = "ABE8e", operator: str = "query-engine") -> list[dict]:
    """
    Thin wrapper around the same extract.build_personalized_cds -> pipeline.run
    path already used by genomic_intake.cli.cmd_run and the web UI's run route.
    """
    subjects.require_consent(subject_id, gene, operator=operator)

    transcript_id = _resolve_transcript(gene)
    personalized_cds = build_personalized_cds(transcript_id, vcf_path)
    guides = pipeline_run(sequence=personalized_cds, editor_name=editor, top_n=10)
    return [asdict(g) for g in guides]


def explain_variant(subject_id: str, gene: str, vcf_path: str, genomic_pos: int, operator: str = "query-engine") -> dict:
    """
    Combines the codon-level consequence pipeline.py already computes with an
    ESM-2 protein-level disruption score for that substitution.
    """
    subjects.require_consent(subject_id, gene, operator=operator)

    transcript_id = _resolve_transcript(gene)
    change = _variant_to_protein_change(transcript_id, vcf_path, genomic_pos)
    variant, codon_idx, ref_aa, alt_aa, protein_seq = (
        change["variant"], change["codon_idx"], change["ref_aa"], change["alt_aa"], change["protein_seq"],
    )

    esm2_score = None
    if ref_aa not in ("?", "*") and alt_aa not in ("?", "*") and codon_idx < len(protein_seq):
        esm2_score = score_substitution(protein_seq, codon_idx + 1, ref_aa, alt_aa)

    return {
        "gene": gene,
        "transcript_id": transcript_id,
        "genomic_pos": genomic_pos,
        "ref": variant.ref,
        "alt": variant.alt,
        "codon_position": codon_idx + 1,
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "is_synonymous": ref_aa == alt_aa,
        "is_nonsense": alt_aa == "*",
        "esm2_disruption_score": esm2_score,
    }


def verify_off_target(subject_id: str, gene: str, vcf_path: str, guide_index: int = 0, editor: str = "ABE8e", operator: str = "query-engine") -> dict:
    """
    Opt-in, one-guide-at-a-time genome-wide off-target check via NCBI BLAST.
    Runs the same rank_guides() path, takes the guide_index'th ranked
    candidate, and submits it for a real BLAST search. Slow (real network
    round trip, often 1+ minutes) — not run automatically for every guide.
    """
    guides = rank_guides(subject_id, gene, vcf_path, editor=editor, operator=operator)
    if not guides:
        raise ValueError(f"No ranked guides available for {gene} to verify")
    if not (0 <= guide_index < len(guides)):
        raise ValueError(f"guide_index {guide_index} out of range (0..{len(guides) - 1})")

    guide = guides[guide_index]
    result = score_off_target(guide["protospacer"], guide["pam_seq"], expected_locus=gene)
    result["guide_index"] = guide_index
    result["protospacer"] = guide["protospacer"]
    return result


def rank_neoantigens(subject_id: str, gene: str, vcf_path: str, genomic_pos: int, hla_alleles: list[str], operator: str = "query-engine") -> list[dict]:
    """
    Somatic variant -> candidate neoantigen peptides (cancer.neoantigen) ->
    MHC-I binding ranking against the given HLA alleles (cancer.mhc_binding).
    hla_alleles are supplied directly by the caller — deriving them from a
    patient's own tumor WES (HLA typing) is out of scope here.
    """
    subjects.require_consent(subject_id, gene, operator=operator)

    transcript_id = _resolve_transcript(gene)
    change = _variant_to_protein_change(transcript_id, vcf_path, genomic_pos)
    codon_idx, protein_seq = change["codon_idx"], change["protein_seq"]

    if codon_idx >= len(protein_seq):
        raise ValueError(f"Variant at codon {codon_idx + 1} falls outside the translated protein (stop codon region)")

    peptides = candidate_peptides(protein_seq, mutated_position=codon_idx + 1)
    if not peptides:
        return []

    return rank_binders([p["peptide"] for p in peptides], hla_alleles)
