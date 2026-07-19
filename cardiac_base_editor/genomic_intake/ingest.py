"""
FASTQ -> VCF ingest, via minimap2 (alignment) + DeepVariant (variant calling,
run through its official Docker image). If a subject's sequencing vendor
already delivers a VCF directly, skip this module entirely and pass that VCF
straight to extract.py.

Requires on PATH: minimap2, samtools, docker (for DeepVariant).
Reference genome (GRCh38 fasta, indexed) must be supplied by the caller —
not downloaded automatically, since it's a multi-GB file best fetched once
and shared across subjects.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from cardiac_base_editor.genomic_intake.storage import derived_dir

DEEPVARIANT_IMAGE = "google/deepvariant:1.6.1"


def align(fastq_path: str, reference_fasta: str, subject_id: str) -> Path:
    """minimap2 align -> sorted, indexed BAM under the subject's derived dir."""
    out_dir = derived_dir(subject_id)
    out_dir.mkdir(parents=True, exist_ok=True)
    sam_path = out_dir / "aligned.sam"
    bam_path = out_dir / "aligned.sorted.bam"

    subprocess.run(
        ["minimap2", "-ax", "map-ont", reference_fasta, fastq_path, "-o", str(sam_path)],
        check=True,
    )
    subprocess.run(
        ["samtools", "sort", "-o", str(bam_path), str(sam_path)],
        check=True,
    )
    subprocess.run(["samtools", "index", str(bam_path)], check=True)
    sam_path.unlink()
    return bam_path


def call_variants(bam_path: Path, reference_fasta: str, subject_id: str) -> Path:
    """Run DeepVariant (via Docker) against the aligned BAM, producing a VCF."""
    out_dir = derived_dir(subject_id)
    vcf_path = out_dir / "variants.vcf.gz"
    ref_dir = Path(reference_fasta).resolve().parent

    subprocess.run(
        [
            "docker", "run",
            "-v", f"{ref_dir}:/ref",
            "-v", f"{out_dir.resolve()}:/out",
            DEEPVARIANT_IMAGE,
            "/opt/deepvariant/bin/run_deepvariant",
            "--model_type=WGS",
            f"--ref=/ref/{Path(reference_fasta).name}",
            f"--reads=/out/{bam_path.name}",
            "--output_vcf=/out/variants.vcf.gz",
            "--num_shards=8",
        ],
        check=True,
    )
    return vcf_path


def ingest_fastq(fastq_path: str, reference_fasta: str, subject_id: str) -> Path:
    """End-to-end: FASTQ -> aligned BAM -> called VCF, all under the subject's storage dir."""
    bam_path = align(fastq_path, reference_fasta, subject_id)
    return call_variants(bam_path, reference_fasta, subject_id)
