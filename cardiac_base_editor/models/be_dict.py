"""
Pluggable local BE-DICT guide-efficiency scorer.

BE-DICT (BEDICT-V2, Schwank Lab / uzh-dqbm-cmi) is a real, published deep
learning model for predicting base-editing outcomes, trained on 38k+ measured
guide RNAs. It has no publicly downloadable pretrained checkpoint as of this
writing — the maintained version is served as a hosted tool at bedict.app,
and running your own copy means training or otherwise obtaining weights on
the actual inference box, not something fetchable from a dev environment.

This module is the seam that model plugs into once available: point
CBE_BEDICT_CHECKPOINT at a local checkpoint file and score_efficiency() will
load and use it. Until then, it falls back to the existing heuristic already
in pipeline.score_guides() — same numbers callers get today, just routed
through one function so the swap-in later is a one-line change (set the env
var), not a code change.
"""

from __future__ import annotations

import os
import warnings
from functools import lru_cache

CHECKPOINT_ENV_VAR = "CBE_BEDICT_CHECKPOINT"

_warned = False


def _checkpoint_path() -> str | None:
    path = os.environ.get(CHECKPOINT_ENV_VAR)
    if path and os.path.exists(path):
        return path
    return None


@lru_cache(maxsize=1)
def _load_model(checkpoint_path: str):
    import torch

    return torch.load(checkpoint_path, map_location="cpu")


def score_efficiency(protospacer: str, editor_name: str, gc_content: float, bystander_risk: bool) -> float:
    """
    Returns a guide efficiency score in [0, 1]. Uses a local BE-DICT
    checkpoint if CBE_BEDICT_CHECKPOINT is set and exists; otherwise falls
    back to the heuristic pipeline.score_guides() already uses, so behavior
    is unchanged until a real checkpoint is provided.
    """
    checkpoint = _checkpoint_path()
    if checkpoint is None:
        global _warned
        if not _warned:
            warnings.warn(
                f"{CHECKPOINT_ENV_VAR} not set (or path doesn't exist) — "
                "falling back to heuristic guide scoring. See models/be_dict.py "
                "docstring for how to plug in a real checkpoint on the inference box.",
                stacklevel=2,
            )
            _warned = True
        return _heuristic_score(protospacer, gc_content, bystander_risk)

    model = _load_model(checkpoint)
    return _run_model(model, protospacer, editor_name)


def _heuristic_score(protospacer: str, gc_content: float, bystander_risk: bool) -> float:
    """Mirrors pipeline.score_guides()'s heuristic exactly, so callers see
    identical numbers whether they go through pipeline.py directly or through
    this seam."""
    gc_penalty = abs(gc_content - 0.525) * 1.8
    bystander_penalty = 0.25 if bystander_risk else 0.0
    polyt_penalty = 0.30 if "TTTT" in protospacer else 0.0
    return max(0.0, 1.0 - gc_penalty - bystander_penalty - polyt_penalty)


def _run_model(model, protospacer: str, editor_name: str) -> float:  # pragma: no cover
    """Not exercised until a real checkpoint exists — the exact input/output
    format depends on which BE-DICT checkpoint Elliot ends up training or
    obtaining, so this is intentionally left as the integration point rather
    than guessed at."""
    raise NotImplementedError(
        "A checkpoint is configured but model inference isn't wired up yet — "
        "implement this once the actual BEDICT-V2 checkpoint format is known."
    )
