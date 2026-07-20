"""
BE-DICT guide-efficiency scorer, backed by BEDICT-V2
(https://github.com/uzh-dqbm-cmi/BEDICT-V2, Schwank Lab, MIT licensed) — a
real, trained deep learning model, not an aspirational placeholder. An
earlier version of this module claimed no public checkpoint existed; that
was never actually verified against the maintained repo and turned out to
be wrong (same mistake as this project initially made about the LNP
delivery predictor, see lnp_delivery.py).

BEDICT-V2 ships real trained PyTorch checkpoints directly in the repo (no
separate download step) for several base-editor variants. Unlike
mrna_design.py/lnp_delivery.py/cancer/tcr_binding.py, this runs IN-PROCESS —
BEDICT-V2's dependencies (torch, pandas, scipy, etc.) are all current-Python
compatible, no legacy-runtime problem, and torch is already a
cardiac_base_editor dependency.

Set CBE_BEDICTV2_DIR to a directory where you've run:
    git clone https://github.com/uzh-dqbm-cmi/BEDICT-V2.git

BEDICT-V2's own code imports matplotlib/seaborn at module load time (only
actually used by its visualization helpers, not by the inference path this
module calls) — not worth adding as a core cardiac_base_editor dependency
for a rarely-exercised path, so it's on you to `pip install matplotlib
seaborn` once alongside the other requirements in BEDICT-V2's own
requirements.txt if CBE_BEDICTV2_DIR is set. Without them, score_efficiency()
degrades to the heuristic (caught by the same try/except that handles any
other real-model failure) rather than crashing guide ranking.

Editor coverage: BEDICT-V2 has trained models for ABE8e-SpCas9 (maps to this
project's "ABE8e", both NGG PAM) and ABE8e-NG (exact name match). It has no
model for "ABE7.10" — that editor always falls back to the heuristic below,
same as when BEDICT-V2 isn't configured at all.

Scoring — real signal, with a real limitation, both worth understanding
before trusting this number:

BEDICT-V2 ships two different models. This module uses the "proportion"
model, which predicts the *relative* distribution among possible bystander-
editing outcomes (which combination of a guide's targetable A's actually got
converted) — not the *absolute* probability that editing happens at all.
That's a different model in the same repo (a CNN, config-file/argparse-
driven rather than a clean importable API, and dependent on extra engineered
features — RNA secondary structure/MFE, melting temperature — that would
need their own external tooling to reproduce correctly). Wrapping that one
was investigated and deliberately deferred; not attempted here.

The practical consequence: for a guide with only one targetable A in its
editing window (the common case), the proportion model's combinatorial
haplotype enumeration only ever produces one possible outcome — the score
below trivially resolves to 1.0 for those guides, regardless of how
favorable the actual sequence context is. It's genuinely informative for
**bystander-risk guides** (multiple targetable A's), where the relative
probability of hitting the intended position vs. an unintended bystander is
exactly what this model is trained to predict. For everything else, this is
real model output, just not distinguishing information — no worse than the
heuristic it replaces, but don't read "1.0" as "the model is confident this
guide works," because for single-target guides it isn't measuring that.
"""

from __future__ import annotations

import os
import sys
import warnings
from functools import lru_cache

BEDICTV2_DIR_ENV_VAR = "CBE_BEDICTV2_DIR"

# Our EDITORS dict keys (pipeline.py) -> BEDICT-V2's trained model directory
# names + (target_nucleotide, conversion_nucleotide) pair it expects.
EDITOR_MODEL_MAP = {
    "ABE8e": ("ABE8e-SpCas9", ("A", "G")),
    "ABE8e-NG": ("ABE8e-NG", ("A", "G")),
    # "ABE7.10" has no BEDICT-V2 model — falls back to heuristic.
}

SEQ_LEN = 24  # 20nt protospacer + 4nt PAM window, matches BEDICT-V2's trained config

_warned_keys = set()


def _bedictv2_dir() -> str | None:
    path = os.environ.get(BEDICTV2_DIR_ENV_VAR)
    if path and os.path.isdir(path):
        return path
    return None


@lru_cache(maxsize=1)
def _bedict_modules():
    """Imports BEDICT-V2's Python modules from the configured checkout.
    Cached — sys.path is mutated once per process, not per call."""
    bedict_dir = _bedictv2_dir()
    if bedict_dir not in sys.path:
        sys.path.insert(0, bedict_dir)

    from utils.sequence_process import SeqProcessConfig, HaplotypeSeqProcessor
    from proportion_model.src.predict_model import BEDICT_EncEnc_HaplotypeModel

    return SeqProcessConfig, HaplotypeSeqProcessor, BEDICT_EncEnc_HaplotypeModel


@lru_cache(maxsize=8)
def _load_predictor(editor_name: str):
    import torch

    model_subdir, conv_nucl = EDITOR_MODEL_MAP[editor_name]
    SeqProcessConfig, HaplotypeSeqProcessor, BEDICT_EncEnc_HaplotypeModel = _bedict_modules()

    seqconfig = SeqProcessConfig(SEQ_LEN, (1, SEQ_LEN), (1, 20), 1)
    seq_processor = HaplotypeSeqProcessor(model_subdir, conv_nucl, seqconfig)
    predictor = BEDICT_EncEnc_HaplotypeModel(seq_processor, seqconfig, torch.device("cpu"))

    model_dir = os.path.join(
        _bedictv2_dir(), "proportion_model", "output",
        "experiment_run_proportions_encenc_two_model",
        f"{model_subdir}_proportions_encenc_two_model",
        "protospacer_PAM", "exp_version_0", "train_val", "run_0",
    )
    return predictor, model_dir


def score_efficiency(
    protospacer: str,
    editor_name: str,
    gc_content: float,
    bystander_risk: bool,
    pam_seq: str | None = None,
    target_positions: list[int] | None = None,
) -> float:
    """
    Returns a guide efficiency score in [0, 1]. Uses the real BEDICT-V2
    proportion model if CBE_BEDICTV2_DIR is configured and editor_name has a
    trained model available; otherwise falls back to the heuristic
    pipeline.score_guides() used before this integration existed.

    See this module's docstring for an important caveat: the real-model path
    is genuinely informative for bystander_risk=True guides (multiple
    targetable A's) but trivially returns 1.0 for single-target guides by
    construction — not a confidence signal for those, just a non-informative
    real number.
    """
    bedict_dir = _bedictv2_dir()
    if bedict_dir is None:
        _warn_once("not_configured",
            f"{BEDICTV2_DIR_ENV_VAR} not set (or path doesn't exist) — "
            "falling back to heuristic guide scoring. See models/be_dict.py "
            "docstring for setup.")
        return _heuristic_score(protospacer, gc_content, bystander_risk)

    if editor_name not in EDITOR_MODEL_MAP or pam_seq is None or not target_positions:
        _warn_once(f"unsupported_{editor_name}",
            f"BEDICT-V2 has no trained model for editor '{editor_name}' "
            "(or missing pam_seq/target_positions) — falling back to heuristic.")
        return _heuristic_score(protospacer, gc_content, bystander_risk)

    try:
        return _run_bedict(protospacer, editor_name, pam_seq, target_positions)
    except Exception as e:  # real model call failed — degrade, don't crash guide ranking
        _warn_once(f"error_{editor_name}", f"BEDICT-V2 scoring failed ({e}) — falling back to heuristic.")
        return _heuristic_score(protospacer, gc_content, bystander_risk)


def _run_bedict(protospacer: str, editor_name: str, pam_seq: str, target_positions: list[int]) -> float:
    import pandas as pd

    predictor, model_dir = _load_predictor(editor_name)

    input_seq = (protospacer + pam_seq).upper()
    if len(input_seq) != SEQ_LEN:
        raise ValueError(f"protospacer+PAM length {len(input_seq)} != expected {SEQ_LEN}")

    df = pd.DataFrame({"ID": ["guide"], "protospacer_PAM": [input_seq]})
    pred_df = predictor.predict_from_dataframe(
        df, ["ID", "protospacer_PAM"], model_dir,
        outpseq_col=None, outcome_col=None, renormalize=True, batch_size=1,
    )

    # target_positions are 1-indexed protospacer positions (pipeline.py's
    # GuideRNA.edit_window_adenines convention) — 0-indexed here to match string slicing.
    target_idx = [p - 1 for p in target_positions]

    # Sums probability mass over haplotypes that edit at least one of our
    # target positions. Note: BEDICT-V2's combinatorial haplotype generator
    # only enumerates outcomes where SOME targetable A got edited (the fully
    # -unedited reference isn't itself a haplotype option) — so when a guide
    # has exactly one targetable A in the whole 24nt window, every generated
    # haplotype edits it by construction, and this trivially sums to 1.0.
    # See this module's docstring for what that does and doesn't tell you.
    on_target_mass = 0.0
    for _, row in pred_df.iterrows():
        inp, outp = row["Inp_seq"], row["Outp_seq"]
        if any(inp[i] != outp[i] for i in target_idx if i < len(inp)):
            on_target_mass += row["pred_score"]

    return min(1.0, max(0.0, float(on_target_mass)))


def _heuristic_score(protospacer: str, gc_content: float, bystander_risk: bool) -> float:
    """The original scoring formula from pipeline.score_guides(), kept here
    as the fallback so both call sites produce identical numbers."""
    gc_penalty = abs(gc_content - 0.525) * 1.8
    bystander_penalty = 0.25 if bystander_risk else 0.0
    polyt_penalty = 0.30 if "TTTT" in protospacer else 0.0
    return max(0.0, 1.0 - gc_penalty - bystander_penalty - polyt_penalty)


def _warn_once(key: str, message: str) -> None:
    if key not in _warned_keys:
        warnings.warn(message, stacklevel=3)
        _warned_keys.add(key)
