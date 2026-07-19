"""
Candidate neoantigen peptide generation: given a protein sequence and the
position of a somatic missense mutation, slide windows of each MHC-I-typical
length across that position to produce candidate peptides for binding
prediction (cancer/mhc_binding.py).

Only the mutated-residue-containing windows are returned — peptides that
don't actually span the mutation aren't neoantigen candidates (they're
identical to a self peptide the immune system already tolerates).
"""

from __future__ import annotations

DEFAULT_LENGTHS = (8, 9, 10, 11)


def candidate_peptides(protein_seq: str, mutated_position: int, lengths: tuple[int, ...] = DEFAULT_LENGTHS) -> list[dict]:
    """
    mutated_position is 1-indexed into protein_seq, matching the
    codon_position convention used elsewhere (pipeline.py, query/engine.py).

    Returns a list of {"peptide": str, "length": int, "start": int, "mutated_index": int}
    — start is the 1-indexed protein position of the peptide's first residue,
    mutated_index is the 0-indexed position of the mutation within the peptide.
    """
    if not (1 <= mutated_position <= len(protein_seq)):
        raise ValueError(f"mutated_position {mutated_position} out of range for sequence of length {len(protein_seq)}")

    zero_idx = mutated_position - 1
    candidates = []
    for length in lengths:
        # Every window of this length that contains the mutated residue.
        earliest_start = max(0, zero_idx - length + 1)
        latest_start = min(len(protein_seq) - length, zero_idx)
        for start in range(earliest_start, latest_start + 1):
            peptide = protein_seq[start:start + length]
            if len(peptide) != length or "?" in peptide or "*" in peptide:
                continue
            candidates.append({
                "peptide": peptide,
                "length": length,
                "start": start + 1,
                "mutated_index": zero_idx - start,
            })
    return candidates
