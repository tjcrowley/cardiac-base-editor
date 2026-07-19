"""
Build a subject's personalized CDS for a target transcript by applying their
called variants (from a VCF) on top of the public reference CDS.

Reuses pipeline.fetch_cds() / pipeline.fetch_gene_info() — no changes to
pipeline.py needed, it already accepts a raw sequence string via run(sequence=...).

Scope (v1): single-nucleotide variants only. Indels change CDS length and
would shift every downstream codon — out of scope until BE-DICT integration
(see pipeline.py roadmap) makes indel-aware scoring meaningful anyway.

Coordinate mapping: genomic VCF coordinates are mapped to CDS-relative
coordinates via the Ensembl REST /map/cds/ endpoint (same API pipeline.py
already depends on for fetch_cds/fetch_gene_info) — no new external service
dependency, no local reference genome download required.
"""

from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass

from cardiac_base_editor.pipeline import ENSEMBL_REST, fetch_cds

COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass
class Variant:
    chrom: str
    pos: int  # 1-based genomic position
    ref: str
    alt: str


def parse_vcf_snvs(vcf_path: str) -> list[Variant]:
    """Minimal VCF parser: yields single-nucleotide variants only, skips
    headers, multi-allelic sites, and indels (len(ref) != 1 or len(alt) != 1)."""
    variants = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt = fields[0], int(fields[1]), fields[2], fields[3], fields[4]
            if "," in alt or len(ref) != 1 or len(alt) != 1:
                continue  # multi-allelic or indel — out of scope for v1
            variants.append(Variant(chrom=chrom, pos=pos, ref=ref, alt=alt))
    return variants


def _fetch_cds_genomic_map(transcript_id: str, cds_length: int) -> tuple[dict[int, int], int]:
    """
    Returns (genomic_pos -> cds_pos map, strand) for the given transcript,
    using Ensembl's /map/cds/ endpoint.
    """
    url = (
        f"{ENSEMBL_REST}/map/cds/{transcript_id}/1..{cds_length}"
        f"?content-type=application/json"
    )
    req = urllib.request.Request(url, headers={"User-Agent": "cardiac-base-editor/1.0"})
    with urllib.request.urlopen(req, timeout=15) as resp:
        data = json.loads(resp.read())

    genomic_to_cds: dict[int, int] = {}
    strand = 1
    cds_cursor = 1
    for mapping in data.get("mappings", []):
        strand = mapping.get("strand", 1)
        g_start, g_end = mapping["start"], mapping["end"]
        span = g_end - g_start + 1
        genomic_positions = (
            range(g_start, g_end + 1) if strand == 1 else range(g_end, g_start - 1, -1)
        )
        for offset, g_pos in enumerate(genomic_positions):
            genomic_to_cds[g_pos] = cds_cursor + offset
        cds_cursor += span
    return genomic_to_cds, strand


def build_personalized_cds(transcript_id: str, vcf_path: str) -> str:
    """
    Fetch the reference CDS, apply all SNVs from vcf_path that fall within
    this transcript's coding sequence, and return the personalized CDS.
    """
    cds = list(fetch_cds(transcript_id))
    genomic_to_cds, strand = _fetch_cds_genomic_map(transcript_id, len(cds))

    variants = parse_vcf_snvs(vcf_path)
    applied = 0
    for v in variants:
        cds_pos = genomic_to_cds.get(v.pos)
        if cds_pos is None:
            continue  # variant not in this transcript's CDS
        alt = v.alt.upper()
        if strand == -1:
            alt = alt.translate(COMPLEMENT)
        cds[cds_pos - 1] = alt
        applied += 1

    print(f"Applied {applied} variant(s) to reference CDS for {transcript_id}")
    return "".join(cds)
