"""
Real (non-mocked) test against the actual pretrained BEDICT-V2 checkpoint.
Requires CBE_BEDICTV2_DIR pointing at a clone of
github.com/uzh-dqbm-cmi/BEDICT-V2 — skipped if not configured, same pattern
as the other live-model tests in this suite.
"""

import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cardiac_base_editor.models import be_dict

# A real 20nt protospacer + NGG PAM (24nt total) with a targetable A in the
# ABE8e editing window (positions 4-8), same convention pipeline.py uses.
TEST_PROTOSPACER = "AAAGTAATTCACTTACAGTC"[:20]
TEST_PAM = "TGGC"
TARGET_POSITIONS = [6]  # 1-indexed position of the targetable A (protospacer[5] == 'A')


def _bedictv2_available() -> bool:
    return be_dict._bedictv2_dir() is not None


def test_falls_back_to_heuristic_when_not_configured(monkeypatch):
    monkeypatch.delenv(be_dict.BEDICTV2_DIR_ENV_VAR, raising=False)
    be_dict._load_predictor.cache_clear()

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        score = be_dict.score_efficiency(TEST_PROTOSPACER, "ABE8e", 0.5, False)

    expected = be_dict._heuristic_score(TEST_PROTOSPACER, 0.5, False)
    assert score == expected


def test_falls_back_to_heuristic_for_unsupported_editor(monkeypatch):
    monkeypatch.setenv(be_dict.BEDICTV2_DIR_ENV_VAR, "/tmp/doesnt-matter-for-this-test")
    score = be_dict.score_efficiency(
        TEST_PROTOSPACER, "ABE7.10", 0.5, False,
        pam_seq=TEST_PAM, target_positions=TARGET_POSITIONS,
    )
    expected = be_dict._heuristic_score(TEST_PROTOSPACER, 0.5, False)
    assert score == expected


def test_falls_back_to_heuristic_when_missing_pam_or_positions(monkeypatch):
    monkeypatch.setenv(be_dict.BEDICTV2_DIR_ENV_VAR, "/tmp/doesnt-matter-for-this-test")
    score = be_dict.score_efficiency(TEST_PROTOSPACER, "ABE8e", 0.5, False)  # no pam_seq/target_positions
    expected = be_dict._heuristic_score(TEST_PROTOSPACER, 0.5, False)
    assert score == expected


@pytest.mark.skipif(not _bedictv2_available(), reason="CBE_BEDICTV2_DIR not configured")
def test_real_model_score_in_valid_range_and_differs_from_heuristic():
    be_dict._load_predictor.cache_clear()
    be_dict._bedict_modules.cache_clear()

    score = be_dict.score_efficiency(
        TEST_PROTOSPACER, "ABE8e", 0.5, False,
        pam_seq=TEST_PAM, target_positions=TARGET_POSITIONS,
    )
    heuristic = be_dict._heuristic_score(TEST_PROTOSPACER, 0.5, False)

    assert 0.0 <= score <= 1.0
    # Not a strict requirement that they differ by construction, but for this
    # input they should — if they're identical, the real model path silently
    # wasn't exercised (e.g. fell back without our noticing).
    assert score != heuristic


@pytest.mark.skipif(not _bedictv2_available(), reason="CBE_BEDICTV2_DIR not configured")
def test_real_model_score_with_both_targetable_adenines():
    """protospacer positions 6 and 7 are both targetable A's within ABE8e's
    editing window (4-8) for this sequence — scoring with both should give a
    score >= scoring with just one (strictly more haplotypes count as
    "on-target edited")."""
    be_dict._load_predictor.cache_clear()
    be_dict._bedict_modules.cache_clear()

    score_one = be_dict.score_efficiency(
        TEST_PROTOSPACER, "ABE8e", 0.5, False, pam_seq=TEST_PAM, target_positions=[6],
    )
    score_both = be_dict.score_efficiency(
        TEST_PROTOSPACER, "ABE8e", 0.5, False, pam_seq=TEST_PAM, target_positions=[6, 7],
    )
    assert 0.0 <= score_one <= 1.0
    assert 0.0 <= score_both <= 1.0
    assert score_both >= score_one
