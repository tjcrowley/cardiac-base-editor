"""
Real (non-mocked) tests against the actual pretrained LiON checkpoint.
Requires CBE_LNPDB_DIR (a clone of github.com/evancollins1/LNPDB) and
CBE_LION_VENV_DIR (a Python 3.7/3.8 venv with chemprop==1.7.0 installed —
see lnp_delivery.py's module docstring for setup). Skipped if either isn't
configured, same pattern as the LinearDesign/Ollama/BLAST/mhcflurry live
tests elsewhere in this suite.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cardiac_base_editor import lnp_delivery

TEST_SMILES = "CCCCCCCCCCCCCCCOC(=O)CCN(CCCN(CCO)CCO)CCC(=O)OCCCCCCCCCCCCCCC"

# The full real formulation this SMILES corresponds to in LNPDB.csv row 1
# (Akinc et al. 2008 siRNA/HeLa formulation) — transcribed from the live
# dataset, not invented. Running this exact formulation through
# predict_delivery_efficacy() should reproduce LNPDB's own bundled
# test_results.csv value for this SMILES, hand-verified to floating-point
# precision before this test suite was written.
KNOWN_FORMULATION = {
    "il_smiles": TEST_SMILES,
    "il_to_nucleicacid_massratio": 5.0,
    "il_molratio": 42.0,
    "hl_molratio": 0.0,
    "chl_molratio": 48.0,
    "peg_molratio": 10.0,
    "dose_ug_nucleicacid": 0.9,
    "hl_name": "None",
    "chl_name": "Cholesterol",
    "peg_name": "DMG-PEG2000",
    "aqueous_buffer": "acetate",
    "dialysis_buffer": "None",
    "mixing_method": "handmixed",
    "model_type": "HeLa",
    "model_target": "in_vitro",
    "route_of_administration": "in_vitro",
    "cargo": "siRNA",
    "cargo_type": "FLuc",
    "experiment_batching": "individual",
}
EXPECTED_VALUE = -0.0858398


def _lion_available() -> bool:
    try:
        lnp_delivery._lnpdb_dir()
        lnp_delivery._chemprop_predict_binary()
        return True
    except lnp_delivery.LiONNotConfigured:
        return False


def test_raises_clear_error_when_lnpdb_dir_not_configured(monkeypatch):
    monkeypatch.delenv("CBE_LNPDB_DIR", raising=False)
    with pytest.raises(lnp_delivery.LiONNotConfigured):
        lnp_delivery.predict_delivery_efficacy([{"il_smiles": TEST_SMILES}])


# ── Feature schema construction (no external tools required) ─────────────

def test_build_extra_x_row_maps_numeric_and_categorical_fields():
    header = [
        "IL_molratio", "HL_molratio",
        "HL_name_DSPC", "HL_name_DOPE", "HL_name_None",
        "Cargo_mRNA", "Cargo_siRNA",
        "Cargo_type_FLuc", "Cargo_type_GFP",
    ]
    formulation = {
        "il_smiles": TEST_SMILES,
        "il_molratio": 42.0,
        "hl_name": "DSPC",
        "cargo": "mRNA",
        "cargo_type": "FLuc",
    }
    row = lnp_delivery._build_extra_x_row(formulation, header)

    assert row["IL_molratio"] == 42.0
    assert row["HL_molratio"] == 0.0  # unspecified numeric defaults to 0.0
    assert row["HL_name_DSPC"] == 1 and row["HL_name_DOPE"] == 0 and row["HL_name_None"] == 0
    assert row["Cargo_mRNA"] == 1 and row["Cargo_siRNA"] == 0
    assert row["Cargo_type_FLuc"] == 1 and row["Cargo_type_GFP"] == 0


def test_build_extra_x_row_all_zero_for_unrecognized_category():
    header = ["HL_name_DSPC", "HL_name_None"]
    formulation = {"il_smiles": TEST_SMILES, "hl_name": "SomeNewLipidNotInSchema"}
    row = lnp_delivery._build_extra_x_row(formulation, header)
    assert row["HL_name_DSPC"] == 0 and row["HL_name_None"] == 0


# ── Real inference against the pretrained checkpoint ──────────────────────

@pytest.mark.skipif(not _lion_available(), reason="CBE_LNPDB_DIR/CBE_LION_VENV_DIR not configured")
def test_predict_delivery_efficacy_matches_known_value():
    results = lnp_delivery.predict_delivery_efficacy([KNOWN_FORMULATION])
    assert len(results) == 1
    assert results[0]["predicted_experiment_value"] == pytest.approx(EXPECTED_VALUE, abs=1e-5)


@pytest.mark.skipif(not _lion_available(), reason="CBE_LNPDB_DIR/CBE_LION_VENV_DIR not configured")
def test_predict_delivery_efficacy_handles_multiple_formulations():
    results = lnp_delivery.predict_delivery_efficacy([
        {"il_smiles": TEST_SMILES},
        {"il_smiles": "CCCCCCCCCCCCCCNC(=O)CCN(CCCN(CCO)CCO)CCC(=O)NCCCCCCCCCCCCCC"},
    ])
    assert len(results) == 2
    assert all("predicted_experiment_value" in r for r in results)
