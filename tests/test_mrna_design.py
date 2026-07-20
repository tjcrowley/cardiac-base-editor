"""
Real (non-mocked) test against a compiled LinearDesign binary. Requires
CBE_LINEARDESIGN_DIR to point at a directory where you've run:
    git clone https://github.com/LinearDesignSoftware/LinearDesign.git
    cd LinearDesign && make
Skipped if that env var isn't set/valid — this is a real external tool this
environment doesn't ship, same as the Ollama/BLAST/mhcflurry live tests.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cardiac_base_editor import mrna_design
from cardiac_base_editor.pipeline import CODON_TABLE
from cardiac_base_editor.plugins import ToolNotConfigured

PROTEIN = "MGKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQ"


def _linear_design_available() -> bool:
    return mrna_design.TOOL.check()


def test_raises_clear_error_when_not_configured(monkeypatch):
    monkeypatch.delenv("CBE_LINEARDESIGN_DIR", raising=False)
    with pytest.raises(ToolNotConfigured):
        mrna_design.design_mrna(PROTEIN)


@pytest.mark.skipif(not _linear_design_available(), reason="CBE_LINEARDESIGN_DIR not configured with a compiled LinearDesign")
def test_design_mrna_produces_valid_translatable_sequence():
    result = mrna_design.design_mrna(PROTEIN)

    assert set(result["mrna_sequence"]) <= set("ACGU")
    assert len(result["mrna_sequence"]) == len(PROTEIN) * 3
    assert result["cai"] > 0
    assert isinstance(result["folding_free_energy_kcal_mol"], float)

    # Translate the returned mRNA back and confirm it encodes the same protein
    dna_seq = result["mrna_sequence"].replace("U", "T")
    translated = "".join(
        CODON_TABLE.get(dna_seq[i:i + 3], "?")
        for i in range(0, len(dna_seq), 3)
    )
    assert translated == PROTEIN
