"""
T-cell receptor / peptide-MHC binding prediction via pMTnet
(https://github.com/tianshilu/pMTnet, Lu et al. 2021, Nature Machine
Intelligence) — the roadmap's "Cancer — T-cell response" step, run after MHC
binding filtering (cancer/mhc_binding.py) to predict whether a T cell would
actually recognize a presented neoantigen.

NOT LIVE-VERIFIED in this environment. pMTnet's own script is plain Python 3
compatible, but its pinned dependencies (TensorFlow 1.x, standalone Keras
2.2.4, numpy 1.16.3) have no wheels for Python 3.11/arm64, and no legacy
Python (3.6/3.7) or pyenv is available here to install them in. This module
is a real, structurally-correct subprocess wrapper — confirmed to construct
the exact input format and command pMTnet's README documents — but the
actual model inference has not been run end-to-end. Whoever stands this up
on real hardware will likely need a Docker container with an old TF1 image
(Elliot's box will also be modern Python/arm64 by default, same constraint).

Setup once a compatible environment exists:
    git clone https://github.com/tianshilu/pMTnet.git
    # then, in a Python 3.6/3.7 env: pip install tensorflow>1.5 keras==2.2.4 numpy==1.16.3 pandas==0.23.4 scikit-learn==0.20.3 scipy==1.2.1
    export CBE_PMTNET_DIR=/path/to/pMTnet
"""

from __future__ import annotations

import csv
import os
import subprocess
import tempfile


class PMTnetNotConfigured(Exception):
    pass


def _pmtnet_dir() -> str:
    path = os.environ.get("CBE_PMTNET_DIR")
    if not path or not os.path.isdir(path):
        raise PMTnetNotConfigured(
            "CBE_PMTNET_DIR is not set (or doesn't exist). To enable T-cell response prediction:\n"
            "  git clone https://github.com/tianshilu/pMTnet.git\n"
            "  # in a Python 3.6/3.7 environment:\n"
            "  pip install 'tensorflow>1.5' keras==2.2.4 numpy==1.16.3 pandas==0.23.4 scikit-learn==0.20.3 scipy==1.2.1\n"
            "  export CBE_PMTNET_DIR=/path/to/pMTnet"
        )
    script = os.path.join(path, "pMTnet.py")
    if not os.path.exists(script):
        raise PMTnetNotConfigured(f"{script} not found — is {path} a real pMTnet checkout?")
    return path


def predict_tcr_binding(records: list[dict], python_executable: str = "python") -> list[dict]:
    """
    records: [{"cdr3": str, "antigen": str, "hla": str}, ...] — TCR-beta CDR3
    sequence, peptide antigen, HLA allele, per pMTnet's documented input
    format.

    Returns each record with an added "rank" field: percentile rank of
    predicted binding strength against 10,000 background TCRs (lower = more
    likely a true TCR-pMHC match, per pMTnet's own scoring convention).
    """
    pmtnet_dir = _pmtnet_dir()

    with tempfile.TemporaryDirectory() as tmp:
        input_csv = os.path.join(tmp, "input.csv")
        output_dir = os.path.join(tmp, "output")
        os.makedirs(output_dir, exist_ok=True)

        with open(input_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["CDR3", "Antigen", "HLA"])
            for r in records:
                writer.writerow([r["cdr3"], r["antigen"], r["hla"]])

        subprocess.run(
            [
                python_executable, os.path.join(pmtnet_dir, "pMTnet.py"),
                "-input", input_csv,
                "-library", os.path.join(pmtnet_dir, "library"),
                "-output", output_dir,
                "-output_log", os.path.join(output_dir, "output.log"),
            ],
            check=True, capture_output=True, text=True, timeout=600,
        )

        output_csv = os.path.join(output_dir, "prediction.csv")
        results = []
        with open(output_csv) as f:
            for row in csv.DictReader(f):
                results.append(dict(row))
        return results
