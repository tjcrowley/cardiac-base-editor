"""
Real (non-mocked) test against the actual downloaded mhcflurry pretrained
models. Requires `mhcflurry-downloads fetch` to have been run once (see
README) — skipped if the models aren't available so this doesn't fail CI
runs that haven't fetched the ~500MB weights.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cardiac_base_editor.cancer import mhc_binding


def _models_available() -> bool:
    try:
        mhc_binding._load_predictor()
        return True
    except Exception:
        return False


@pytest.mark.skipif(not _models_available(), reason="mhcflurry pretrained models not downloaded (run: mhcflurry-downloads fetch)")
def test_known_strong_binder_scores_better_than_random_peptide():
    """
    NLVPMVATV / HLA-A*02:01 is the classic CMV pp65 epitope used throughout
    mhcflurry's own reference examples — a well-established strong binder.
    Compared against an arbitrary peptide with no particular affinity for
    A*02:01, to sanity-check the wrapper produces sensible relative rankings
    before trusting it on synthetic variant-derived peptides.
    """
    results = mhc_binding.rank_binders(
        peptides=["NLVPMVATV", "AAAAAAAAA"],
        hla_alleles=["A0201"],
    )
    by_peptide = {r["peptide"]: r for r in results}

    assert by_peptide["NLVPMVATV"]["affinity_percentile"] < by_peptide["AAAAAAAAA"]["affinity_percentile"]
    assert by_peptide["NLVPMVATV"]["affinity_nm"] < 500  # well-established strong binder territory


@pytest.mark.skipif(not _models_available(), reason="mhcflurry pretrained models not downloaded (run: mhcflurry-downloads fetch)")
def test_rank_binders_sorted_ascending_by_percentile():
    results = mhc_binding.rank_binders(
        peptides=["NLVPMVATV", "AAAAAAAAA", "SIINFEKL"],
        hla_alleles=["A0201"],
    )
    percentiles = [r["affinity_percentile"] for r in results]
    assert percentiles == sorted(percentiles)
