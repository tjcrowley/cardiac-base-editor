"""
HLA typing from RNA-seq via arcasHLA
(https://github.com/RabadanLab/arcasHLA, Orenbuch et al. 2019, Bioinformatics)
— derives a patient's own HLA alleles from their tumor RNA-seq, rather than
requiring them to be supplied directly (cancer/mhc_binding.py and
cancer/tcr_binding.py both currently take hla_alleles as a caller-supplied
argument).

NOT TESTABLE in this environment by construction: arcasHLA genotypes real
RNA-seq alignments (BAM), and there's no synthetic substitute that produces
a meaningful HLA call — a fabricated "toy RNA-seq sample" wouldn't type to
anything real. This is the same situation as genomic_intake/ingest.py's
minimap2/DeepVariant wrappers: real, structurally-correct subprocess code,
verified only by inspection, not a live run — box-only until real tumor
RNA-seq data exists to run it against.

Setup:
    git clone https://github.com/RabadanLab/arcasHLA.git
    # or: conda install -c bioconda arcas-hla
    cd arcasHLA && arcasHLA reference --update
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile


class ArcasHLANotConfigured(Exception):
    pass


def _arcashla_binary(arcashla_dir: str | None) -> str:
    if arcashla_dir:
        binary = os.path.join(arcashla_dir, "arcasHLA")
        if not os.path.exists(binary):
            raise ArcasHLANotConfigured(f"{binary} not found in {arcashla_dir}")
        return binary

    from shutil import which
    binary = which("arcasHLA")
    if not binary:
        raise ArcasHLANotConfigured(
            "arcasHLA not found on PATH and no arcashla_dir given. Setup:\n"
            "  git clone https://github.com/RabadanLab/arcasHLA.git\n"
            "  # or: conda install -c bioconda arcas-hla\n"
            "  cd arcasHLA && arcasHLA reference --update"
        )
    return binary


def type_hla_from_rnaseq(bam_path: str, arcashla_dir: str | None = None) -> list[str]:
    """
    Runs `arcasHLA genotype` against an aligned RNA-seq BAM, returns the
    called HLA-I alleles (e.g. ["A*02:01", "A*11:01", "B*07:02", ...]).
    """
    binary = _arcashla_binary(arcashla_dir)

    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(
            [binary, "genotype", bam_path, "-o", tmp, "--log", os.path.join(tmp, "genotype.log")],
            check=True, capture_output=True, text=True, timeout=1800,
        )

        genotype_files = [f for f in os.listdir(tmp) if f.endswith(".genotype.json")]
        if not genotype_files:
            raise RuntimeError(f"arcasHLA produced no genotype.json output in {tmp}")

        with open(os.path.join(tmp, genotype_files[0])) as f:
            genotype = json.load(f)

        alleles = []
        for gene_alleles in genotype.values():
            alleles.extend(gene_alleles)
        return alleles
