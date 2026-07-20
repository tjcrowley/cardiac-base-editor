"""
LNP (lipid nanoparticle) delivery-efficacy prediction via LiON
(Witten et al. 2024, Nat Biotech), trained on the public LNPDB dataset
(https://github.com/evancollins1/LNPDB, Nat Comms 2026).

Unlike everything else in this package, this module does NOT touch subject
genomic data — it predicts delivery efficacy for a candidate LNP
*formulation* (lipid structure + mixing ratios + assay conditions), so it
isn't consent-gated and isn't part of query/engine.py's subject-scoped
function table. It's exposed as a standalone `cbe lnp-predict` CLI verb.

LiON is a standard `chemprop` (message-passing GNN) regression model, but the
exact version it was trained with — chemprop==1.7.0 — requires Python
3.7/3.8, incompatible with this package's own Python 3.10+ requirement. It
needs its own isolated environment:

    # one-time setup
    brew install pyenv          # or your platform's pyenv equivalent
    pyenv install 3.8.18
    ~/.pyenv/versions/3.8.18/bin/python3.8 -m venv ~/lion-venv
    ~/lion-venv/bin/pip install chemprop==1.7.0

    git clone https://github.com/evancollins1/LNPDB.git

    export CBE_LNPDB_DIR=/path/to/LNPDB
    export CBE_LION_VENV_DIR=~/lion-venv

Note on the predicted value: LNPDB's "Experiment_value" is the dataset's own
normalized/pooled delivery-efficacy metric across heterogeneous assay types
(luminescence, flow cytometry, barcode sequencing, etc.) — a real, correctly
computed prediction, but interpreting the number requires the LNPDB paper's
methodology, not a raw physical unit.
"""

from __future__ import annotations

import csv
import os
import subprocess
import tempfile
from functools import lru_cache

from cardiac_base_editor.plugins import ExternalTool, register_tool, require_configured

CHECKPOINT_SUBPATH = "data/LNPDB_for_LiON/single_split/trained_model_checkpoints"
SCHEMA_CSV_SUBPATH = "data/LNPDB_for_LiON/single_split/train_extra_x.csv"

NUMERIC_COLUMNS = [
    "IL_to_nucleicacid_massratio", "IL_molratio", "HL_molratio",
    "CHL_molratio", "PEG_molratio", "Dose_ug_nucleicacid",
]

# Order matters: "cargo_type" must be checked before "cargo" since
# "Cargo_type_*" columns also start with "Cargo_".
CATEGORICAL_GROUPS = [
    ("cargo_type", "Cargo_type_"),
    ("hl_name", "HL_name_"),
    ("chl_name", "CHL_name_"),
    ("peg_name", "PEG_name_"),
    ("aqueous_buffer", "Aqueous_buffer_"),
    ("dialysis_buffer", "Dialysis_buffer_"),
    ("mixing_method", "Mixing_method_"),
    ("model_type", "Model_type_"),
    ("model_target", "Model_target_"),
    ("route_of_administration", "Route_of_administration_"),
    ("cargo", "Cargo_"),
    ("experiment_batching", "Experiment_batching_"),
]


def _lnpdb_ready() -> bool:
    path = os.environ.get("CBE_LNPDB_DIR")
    if not path or not os.path.isdir(path):
        return False
    return os.path.isdir(os.path.join(path, CHECKPOINT_SUBPATH))


def _lion_venv_ready() -> bool:
    venv = os.environ.get("CBE_LION_VENV_DIR")
    if not venv or not os.path.isdir(venv):
        return False
    return os.path.exists(os.path.join(venv, "bin", "chemprop_predict"))


TOOL = register_tool(ExternalTool(
    name="LiON",
    env_vars=["CBE_LNPDB_DIR", "CBE_LION_VENV_DIR"],
    setup_instructions=(
        "git clone https://github.com/evancollins1/LNPDB.git\n"
        "  export CBE_LNPDB_DIR=$(pwd)/LNPDB\n"
        "  pyenv install 3.8.18\n"
        "  ~/.pyenv/versions/3.8.18/bin/python3.8 -m venv ~/lion-venv\n"
        "  ~/lion-venv/bin/pip install chemprop==1.7.0\n"
        "  export CBE_LION_VENV_DIR=~/lion-venv"
    ),
    check=lambda: _lnpdb_ready() and _lion_venv_ready(),
))


def _lnpdb_dir() -> str:
    return os.environ["CBE_LNPDB_DIR"]


def _chemprop_predict_binary() -> str:
    return os.path.join(os.environ["CBE_LION_VENV_DIR"], "bin", "chemprop_predict")


@lru_cache(maxsize=1)
def _load_feature_schema() -> list[str]:
    """Reads the exact one-hot/numeric column order LiON was trained on."""
    schema_path = os.path.join(_lnpdb_dir(), SCHEMA_CSV_SUBPATH)
    with open(schema_path) as f:
        return next(csv.reader(f))


def _build_extra_x_row(formulation: dict, header: list[str]) -> dict:
    """
    Expands a human-friendly formulation dict into the full one-hot/numeric
    feature row LiON expects, in header order. Unrecognized/unspecified
    categorical values default to all-zero for that category; missing
    numeric ratios default to 0.0.
    """
    row = {}
    for col in header:
        if col in NUMERIC_COLUMNS:
            key = col.lower()
            row[col] = formulation.get(key, 0.0)
            continue

        for group_key, prefix in CATEGORICAL_GROUPS:
            if col.startswith(prefix):
                suffix = col[len(prefix):]
                row[col] = 1 if str(formulation.get(group_key, "")) == suffix else 0
                break
        else:
            row[col] = 0  # column doesn't match any known group — leave unset
    return row


def predict_delivery_efficacy(formulations: list[dict], no_cuda: bool = True) -> list[dict]:
    """
    formulations: list of dicts, each with at minimum "il_smiles", plus any
    of the formulation parameters described in this module's docstring
    (numeric ratios + categorical group keys). Returns each formulation with
    an added "predicted_experiment_value" field.
    """
    require_configured(TOOL)
    lnpdb_dir = _lnpdb_dir()
    chemprop_predict = _chemprop_predict_binary()
    header = _load_feature_schema()

    with tempfile.TemporaryDirectory() as tmp:
        test_csv = os.path.join(tmp, "test.csv")
        extra_x_csv = os.path.join(tmp, "extra_x.csv")
        preds_csv = os.path.join(tmp, "preds.csv")

        with open(test_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["IL_SMILES"])
            for formulation in formulations:
                writer.writerow([formulation["il_smiles"]])

        with open(extra_x_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            for formulation in formulations:
                writer.writerow(_build_extra_x_row(formulation, header))

        cmd = [
            chemprop_predict,
            "--checkpoint_dir", os.path.join(lnpdb_dir, CHECKPOINT_SUBPATH),
            "--test_path", test_csv,
            "--features_path", extra_x_csv,
            "--preds_path", preds_csv,
        ]
        if no_cuda:
            cmd.append("--no_cuda")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            raise RuntimeError(f"chemprop_predict exited {result.returncode}: {result.stderr[-1000:]}")

        with open(preds_csv) as f:
            preds = list(csv.DictReader(f))

    return [
        {**formulation, "predicted_experiment_value": float(pred["Experiment_value"])}
        for formulation, pred in zip(formulations, preds)
    ]
