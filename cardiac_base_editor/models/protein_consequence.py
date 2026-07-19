"""
ESM-2 masked-marginal variant-effect scoring.

Given a protein sequence and a single amino-acid substitution, scores how
disruptive that substitution is using the standard lightweight ESM
variant-effect approach: mask the position, read the model's log-probability
for the reference amino acid and for the alternate amino acid at that masked
position, and take the difference (alt - ref). A large negative score means
the model finds the alternate amino acid much less likely than the reference
at that position — i.e. a more disruptive substitution.

Model: facebook/esm2_t6_8M_UR50D — the smallest public ESM-2 checkpoint
(8M params), chosen so this runs on CPU without requiring the inference
box's GPU. Swap ESM2_MODEL_NAME for a larger checkpoint once running on
hardware with more headroom, if higher accuracy is worth the latency.
"""

from __future__ import annotations

import os
from functools import lru_cache

import torch

ESM2_MODEL_NAME = os.environ.get("CBE_ESM2_MODEL", "facebook/esm2_t6_8M_UR50D")


@lru_cache(maxsize=1)
def _load_model():
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    tokenizer = AutoTokenizer.from_pretrained(ESM2_MODEL_NAME)
    model = AutoModelForMaskedLM.from_pretrained(ESM2_MODEL_NAME)
    model.eval()
    return tokenizer, model


def score_substitution(protein_seq: str, position: int, ref_aa: str, alt_aa: str) -> float:
    """
    position is 1-indexed into protein_seq (matches pipeline.py's codon_position
    convention in AminoAcidConsequence).

    Returns the masked-marginal log-probability difference log P(alt) - log P(ref)
    at that position. Negative = alt is less likely than ref (more disruptive).
    Returns 0.0 for stop-codon consequences ("*") — not a substitution ESM-2's
    vocabulary scores meaningfully.
    """
    if ref_aa == "*" or alt_aa == "*":
        return 0.0
    if not (1 <= position <= len(protein_seq)):
        raise ValueError(f"position {position} out of range for sequence of length {len(protein_seq)}")

    tokenizer, model = _load_model()

    masked_seq = protein_seq[: position - 1] + tokenizer.mask_token + protein_seq[position:]
    inputs = tokenizer(masked_seq, return_tensors="pt")
    mask_token_index = (inputs["input_ids"][0] == tokenizer.mask_token_id).nonzero(as_tuple=True)[0]

    with torch.no_grad():
        logits = model(**inputs).logits

    log_probs = torch.log_softmax(logits[0, mask_token_index[0]], dim=-1)

    ref_id = tokenizer.convert_tokens_to_ids(ref_aa)
    alt_id = tokenizer.convert_tokens_to_ids(alt_aa)
    if ref_id == tokenizer.unk_token_id or alt_id == tokenizer.unk_token_id:
        raise ValueError(f"amino acid not in ESM-2 vocabulary: ref={ref_aa} alt={alt_aa}")

    return float(log_probs[alt_id] - log_probs[ref_id])
