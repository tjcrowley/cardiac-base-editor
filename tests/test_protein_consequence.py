"""
Real test against the actual downloaded ESM-2 model (facebook/esm2_t6_8M_UR50D,
8M params) — no mocking. Confirms a disruptive substitution scores lower than
a conservative one, and that stop-codon consequences short-circuit to 0.0
without needing a valid ESM-2 token for '*'.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cardiac_base_editor.models import protein_consequence as pc

# A short, unremarkable protein context to substitute within.
PROTEIN = "MGKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQ"


def test_stop_codon_short_circuits_without_model_call():
    assert pc.score_substitution(PROTEIN, 10, "V", "*") == 0.0
    assert pc.score_substitution(PROTEIN, 10, "*", "V") == 0.0


def test_position_out_of_range_raises():
    with pytest.raises(ValueError):
        pc.score_substitution(PROTEIN, len(PROTEIN) + 5, "V", "A")


def test_synonymous_like_substitution_scores_higher_than_charge_reversal():
    # Position 10 is 'V' (Valine). Compare a conservative swap (V->I, both
    # small hydrophobic) against a disruptive one (V->D, hydrophobic->charged).
    conservative = pc.score_substitution(PROTEIN, 10, "V", "I")
    disruptive = pc.score_substitution(PROTEIN, 10, "V", "D")
    assert conservative > disruptive
