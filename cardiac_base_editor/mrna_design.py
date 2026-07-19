"""
Codon-optimized mRNA sequence design via LinearDesign
(https://github.com/LinearDesignSoftware/LinearDesign, Zhang et al. 2023,
Nature) — the final step in the pipeline's roadmap: turning a target protein
sequence into an actual mRNA payload sequence.

LinearDesign ships as C++ source you compile yourself (`make`), plus a thin
python2 CLI wrapper. That wrapper is trivial — it just forwards
(lambda, verbose, codon_usage_csv) to the compiled binary over stdin/stdout —
so this module calls the compiled binary directly and skips the python2
dependency entirely:

    echo PROTEIN_SEQ | bin/LinearDesign_2D <lambda> <verbose 0|1> <codon_usage.csv>

Set CBE_LINEARDESIGN_DIR to a directory where you've run:
    git clone https://github.com/LinearDesignSoftware/LinearDesign.git
    cd LinearDesign && make
"""

from __future__ import annotations

import os
import re
import subprocess

CODON_USAGE_FILES = {
    "human": "codon_usage_freq_table_human.csv",
    "yeast": "codon_usage_freq_table_yeast.csv",
}


class LinearDesignNotConfigured(Exception):
    pass


def _linear_design_dir() -> str:
    path = os.environ.get("CBE_LINEARDESIGN_DIR")
    if not path or not os.path.isdir(path):
        raise LinearDesignNotConfigured(
            "CBE_LINEARDESIGN_DIR is not set (or doesn't exist). To enable mRNA design:\n"
            "  git clone https://github.com/LinearDesignSoftware/LinearDesign.git\n"
            "  cd LinearDesign && make\n"
            "  export CBE_LINEARDESIGN_DIR=$(pwd)"
        )
    binary = os.path.join(path, "bin", "LinearDesign_2D")
    if not os.path.exists(binary):
        raise LinearDesignNotConfigured(
            f"{binary} not found — run `make` in {path} first."
        )
    return path


_OUTPUT_PATTERN = re.compile(
    r"mRNA sequence:\s*(?P<sequence>[ACGU]+)\s*\n"
    r"mRNA structure:\s*(?P<structure>[().]+)\s*\n"
    r"mRNA folding free energy:\s*(?P<energy>-?\d+\.?\d*)\s*kcal/mol;\s*mRNA CAI:\s*(?P<cai>\d+\.?\d*)"
)


def design_mrna(protein_seq: str, lambda_: float = 0.0, codon_usage: str = "human", verbose: bool = False, timeout_s: int = 900) -> dict:
    """
    Runs LinearDesign on protein_seq, returns
    {"mrna_sequence", "mrna_structure", "folding_free_energy_kcal_mol", "cai"}.

    Runtime scales with protein length — the LinearDesign paper reports ~11
    minutes for the 1273-aa SARS-CoV-2 spike protein at default settings, so
    the default timeout here is generous (15 min) rather than tuned for the
    short toy sequences unit tests use.
    """
    ld_dir = _linear_design_dir()
    binary = os.path.join(ld_dir, "bin", "LinearDesign_2D")
    codon_usage_path = os.path.join(ld_dir, CODON_USAGE_FILES.get(codon_usage, codon_usage))

    result = subprocess.run(
        [binary, str(lambda_), "1" if verbose else "0", codon_usage_path],
        input=protein_seq, capture_output=True, text=True, timeout=timeout_s,
        cwd=ld_dir,  # LinearDesign_2D references its .so via a path relative to cwd, not the binary
    )
    if result.returncode != 0:
        raise RuntimeError(f"LinearDesign_2D exited {result.returncode}: {result.stderr[-500:]}")

    match = _OUTPUT_PATTERN.search(result.stdout)
    if not match:
        raise RuntimeError(f"Couldn't parse LinearDesign output: {result.stdout[-500:]}")

    return {
        "mrna_sequence": match.group("sequence"),
        "mrna_structure": match.group("structure"),
        "folding_free_energy_kcal_mol": float(match.group("energy")),
        "cai": float(match.group("cai")),
    }
