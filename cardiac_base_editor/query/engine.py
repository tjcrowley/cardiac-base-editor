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
from cardiac_base_editor.models.protein_consequence import score_substitution
from cardiac_base_editor.pipeline import (
    CODON_TABLE,
    KNOWN_TARGETS,
    fetch_cds,
    run as pipeline_run,
)


def _resolve_transcript(gene: str) -> str:
    return KNOWN_TARGETS.get(gene.upper(), gene)


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
