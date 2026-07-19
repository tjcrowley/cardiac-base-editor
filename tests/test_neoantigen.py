import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cardiac_base_editor.cancer import neoantigen

PROTEIN = "MGKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQ"


def test_candidate_peptides_all_contain_mutation():
    mutated_position = 10  # 'V' in this toy protein
    peptides = neoantigen.candidate_peptides(PROTEIN, mutated_position)
    assert peptides
    for p in peptides:
        zero_idx = mutated_position - 1
        peptide_start = p["start"] - 1
        assert peptide_start <= zero_idx < peptide_start + p["length"]
        assert p["mutated_index"] == zero_idx - peptide_start


def test_candidate_peptides_respects_requested_lengths():
    peptides = neoantigen.candidate_peptides(PROTEIN, 10, lengths=(9,))
    assert all(p["length"] == 9 for p in peptides)
    assert len(peptides) == 9  # 9 windows of length 9 contain a given position


def test_candidate_peptides_handles_near_terminus():
    # Mutation right at the start of the protein — windows can't extend before position 1.
    peptides = neoantigen.candidate_peptides(PROTEIN, 1, lengths=(9,))
    assert peptides
    assert all(p["start"] == 1 for p in peptides)


def test_candidate_peptides_raises_out_of_range():
    with pytest.raises(ValueError):
        neoantigen.candidate_peptides(PROTEIN, len(PROTEIN) + 10)
