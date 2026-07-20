"""
Cardiac Base Editing Pipeline — Inference Prototype

Given a target gene (by Ensembl transcript ID or raw sequence), find and rank
guide RNAs for adenine base editing (ABE), predict editing outcomes, and flag
amino acid consequences.

Full flow:
  Ensembl transcript ID
    → fetch CDS via REST API  (or pass sequence directly)
    → guide RNA candidates (PAM scan, both strands)
    → filter: editing window must contain targetable A
    → score: efficiency + off-target risk
    → annotate: codon consequence of each A→G edit
    → safety filter: drop guides that introduce stop codons
    → rank and report

Dependencies:
    pip install requests
"""

from __future__ import annotations

import re
import time
import urllib.request
import json
from dataclasses import dataclass, field
from typing import Literal

from cardiac_base_editor.models import be_dict


# ─────────────────────────────────────────────
# Ensembl REST API
# ─────────────────────────────────────────────

ENSEMBL_REST = "https://rest.ensembl.org"

def fetch_cds(transcript_id: str) -> str:
    """
    Fetch the coding sequence (CDS) for an Ensembl transcript.
    Returns the sequence as a plain uppercase DNA string.

    Examples:
        PCSK9  → ENST00000302118
        LDLR   → ENST00000558013
        ANGPTL3 → ENST00000264027

    Rate limit: Ensembl allows 15 req/s unauthenticated.
    """
    url = (
        f"{ENSEMBL_REST}/sequence/id/{transcript_id}"
        f"?content-type=application/json&type=cds"
    )
    req = urllib.request.Request(url, headers={"User-Agent": "cardiac-base-editor/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read())
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"Ensembl API error {e.code} for {transcript_id}: {e.reason}")

    seq = data.get("seq", "")
    if not seq:
        raise RuntimeError(f"No sequence returned for {transcript_id}")
    return seq.upper()


def fetch_gene_info(transcript_id: str) -> dict:
    """Return basic metadata for a transcript (gene name, description, chromosome)."""
    url = (
        f"{ENSEMBL_REST}/lookup/id/{transcript_id}"
        f"?content-type=application/json&expand=1"
    )
    req = urllib.request.Request(url, headers={"User-Agent": "cardiac-base-editor/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError:
        return {}


# ─────────────────────────────────────────────
# Base editor configurations
# ─────────────────────────────────────────────

@dataclass
class BaseEditor:
    name: str
    pam: str       # IUPAC string, e.g. "NGG" or "NG"
    window: range  # 1-indexed protospacer positions where A gets edited
    edit: tuple[str, str]

# ABE8e: the editor used in VERVE-102. Converts A→G at positions 4–8.
# ABE8e-NG: relaxed PAM (NG instead of NGG), wider target space.
EDITORS = {
    "ABE8e":    BaseEditor("ABE8e",    pam="NGG", window=range(4, 9), edit=("A", "G")),
    "ABE8e-NG": BaseEditor("ABE8e-NG", pam="NG",  window=range(4, 9), edit=("A", "G")),
    "ABE7.10":  BaseEditor("ABE7.10",  pam="NGG", window=range(4, 8), edit=("A", "G")),
}


# ─────────────────────────────────────────────
# Data types
# ─────────────────────────────────────────────

@dataclass
class AminoAcidConsequence:
    codon_position: int
    ref_codon: str
    alt_codon: str
    ref_aa: str
    alt_aa: str
    is_synonymous: bool
    is_nonsense: bool
    is_lof: bool

@dataclass
class GuideRNA:
    protospacer: str
    pam_seq: str
    strand: Literal["+", "-"]
    genomic_position: int
    editor: str
    edit_window_adenines: list[int]
    gc_content: float
    context_after_pam: str = ""  # one flanking base beyond the PAM, when available — see be_dict.py
    efficiency_score: float | None = None
    off_target_score: float | None = None
    bystander_risk: bool = False
    consequences: list[AminoAcidConsequence] = field(default_factory=list)


# ─────────────────────────────────────────────
# Step 1: Guide discovery
# ─────────────────────────────────────────────

def _pam_regex(pam: str) -> re.Pattern:
    iupac = {"N": "[ACGT]", "A": "A", "C": "C", "G": "G", "T": "T",
             "R": "[AG]", "Y": "[CT]", "W": "[AT]", "S": "[GC]"}
    pattern = "".join(iupac.get(c, c) for c in pam)
    # Trailing [ACGT]? captures one optional base of context after the PAM —
    # BEDICT-V2 (models/be_dict.py) was trained on 24nt windows (20nt
    # protospacer + PAM + 1 flanking base), one more than the PAM alone.
    # Optional so guides at the very end of a sequence (no base available
    # after the PAM) are still found, just without that extra context.
    return re.compile(r"(?=([ACGT]{20}" + pattern + r"[ACGT]?))", re.IGNORECASE)

def _reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def find_guides(sequence: str, editor: BaseEditor) -> list[GuideRNA]:
    """Scan both strands for PAM sites; return candidates with A in editing window."""
    seq = sequence.upper()
    rc  = _reverse_complement(seq)
    pat = _pam_regex(editor.pam)
    guides: list[GuideRNA] = []

    for strand, s, offset_fn in [
        ("+", seq, lambda pos: pos),
        ("-", rc,  lambda pos: len(seq) - pos - 20 - len(editor.pam)),
    ]:
        for match in pat.finditer(s):
            hit   = match.group(1)
            proto = hit[:20]
            pam_s = hit[20:20 + len(editor.pam)]
            context_after_pam = hit[20 + len(editor.pam):20 + len(editor.pam) + 1]

            targetable = [
                i + 1 for i, nt in enumerate(proto)
                if nt == editor.edit[0] and (i + 1) in editor.window
            ]
            if not targetable:
                continue

            gc = (proto.count("G") + proto.count("C")) / 20
            guides.append(GuideRNA(
                protospacer=proto,
                pam_seq=pam_s,
                strand=strand,
                genomic_position=offset_fn(match.start()),
                editor=editor.name,
                edit_window_adenines=targetable,
                gc_content=gc,
                context_after_pam=context_after_pam,
                bystander_risk=len(targetable) > 1,
            ))

    return guides


# ─────────────────────────────────────────────
# Step 2: Scoring
# ─────────────────────────────────────────────

def score_guides(guides: list[GuideRNA]) -> list[GuideRNA]:
    """
    Efficiency scoring via models.be_dict.score_efficiency() — uses the real
    BEDICT-V2 model when CBE_BEDICTV2_DIR is configured and the editor has a
    trained model available, otherwise the heuristic formula that used to
    live inline here (now in be_dict._heuristic_score(), single source of
    truth for both call sites).
    """
    for g in guides:
        g.efficiency_score = be_dict.score_efficiency(
            g.protospacer, g.editor, g.gc_content, g.bystander_risk,
            pam_seq=g.pam_seq + g.context_after_pam, target_positions=g.edit_window_adenines,
        )

        # Seed region = PAM-proximal 12 nt
        seed = g.protospacer[8:]
        homopolymer_runs = sum(
            1 for i in range(len(seed) - 3)
            if len(set(seed[i:i+4])) == 1
        )
        g.off_target_score = min(1.0, homopolymer_runs * 0.2)

    return sorted(guides, key=lambda g: (
        -(g.efficiency_score or 0),
         (g.off_target_score or 1),
    ))


# ─────────────────────────────────────────────
# Step 3: Amino acid consequence
# ─────────────────────────────────────────────

CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

def annotate_consequences(guide: GuideRNA, cds: str, cds_offset: int = 0) -> GuideRNA:
    """
    For each targetable A in the editing window, compute the amino acid
    consequence of A→G in the CDS.

    cds_offset: where the guide's protospacer starts within cds.
    Minus-strand guides: the edit is on the complementary strand, so the CDS
    sees a T→C change — handled by the reverse-complement scan in find_guides.
    """
    guide.consequences = []
    for pos in guide.edit_window_adenines:
        genomic_idx   = cds_offset + (pos - 1)
        codon_idx     = genomic_idx // 3
        base_in_codon = genomic_idx % 3

        if codon_idx * 3 + 3 > len(cds):
            continue

        ref_codon = cds[codon_idx*3 : codon_idx*3 + 3]
        alt_codon = ref_codon[:base_in_codon] + "G" + ref_codon[base_in_codon+1:]
        ref_aa    = CODON_TABLE.get(ref_codon, "?")
        alt_aa    = CODON_TABLE.get(alt_codon, "?")

        guide.consequences.append(AminoAcidConsequence(
            codon_position=codon_idx + 1,
            ref_codon=ref_codon,
            alt_codon=alt_codon,
            ref_aa=ref_aa,
            alt_aa=alt_aa,
            is_synonymous=ref_aa == alt_aa,
            is_nonsense=alt_aa == "*",
            is_lof=alt_aa == "*",
        ))
    return guide


# ─────────────────────────────────────────────
# Step 4: Safety filter
# ─────────────────────────────────────────────

def filter_unsafe(guides: list[GuideRNA]) -> tuple[list[GuideRNA], list[GuideRNA]]:
    """Hard filter: drop guides that introduce stop codons or have high off-target score."""
    safe, rejected = [], []
    for g in guides:
        if any(c.is_nonsense for c in g.consequences):
            rejected.append(g)
        elif (g.off_target_score or 0) > 0.6:
            rejected.append(g)
        else:
            safe.append(g)
    return safe, rejected


# ─────────────────────────────────────────────
# Pipeline entry point
# ─────────────────────────────────────────────

def run(
    transcript_id: str | None = None,
    sequence: str | None = None,
    editor_name: str = "ABE8e",
    top_n: int = 10,
) -> list[GuideRNA]:
    """
    End-to-end cardiac base editing pipeline.

    transcript_id: Ensembl transcript ID (e.g. 'ENST00000302118' for PCSK9).
                   If provided, CDS is fetched automatically from Ensembl REST.
    sequence:      Raw CDS string. Used directly if transcript_id is not given.
    editor_name:   'ABE8e' | 'ABE8e-NG' | 'ABE7.10'
    top_n:         Number of final ranked guides to return.
    """
    if transcript_id:
        print(f"Fetching CDS from Ensembl: {transcript_id} ... ", end="", flush=True)
        cds = fetch_cds(transcript_id)
        print(f"{len(cds)} bp")
    elif sequence:
        cds = sequence.upper()
        print(f"Using provided sequence: {len(cds)} bp")
    else:
        raise ValueError("Provide either transcript_id or sequence.")

    editor = EDITORS[editor_name]
    print(f"Editor: {editor.name}  PAM: {editor.pam}  Window: {list(editor.window)}\n")

    guides = find_guides(cds, editor)
    print(f"Candidates found: {len(guides)}")

    guides = score_guides(guides)

    guides = [annotate_consequences(g, cds, g.genomic_position) for g in guides]

    safe, rejected = filter_unsafe(guides)
    if rejected:
        print(f"Filtered out:     {len(rejected)} (stop codons or high off-target)")
    print(f"Passing safety:   {len(safe)}\n")

    return safe[:top_n]


# ─────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────

KNOWN_TARGETS = {
    "PCSK9":   "ENST00000302118",  # primary LDL-lowering target (VERVE-102)
    "LDLR":    "ENST00000558013",  # LDL receptor — loss-of-function causes FH
    "ANGPTL3": "ENST00000264027",  # triglyceride regulator, pan-lipid target
    "APOB":    "ENST00000233242",  # ApoB — alternative LDL-lowering target
}

if __name__ == "__main__":
    import sys

    # Usage: python pipeline.py [GENE_NAME_OR_TRANSCRIPT_ID] [EDITOR]
    # e.g.:  python pipeline.py PCSK9
    #        python pipeline.py ENST00000302118 ABE8e-NG
    target_arg = sys.argv[1] if len(sys.argv) > 1 else "PCSK9"
    editor_arg = sys.argv[2] if len(sys.argv) > 2 else "ABE8e"

    transcript = KNOWN_TARGETS.get(target_arg.upper(), target_arg)
    print(f"Target: {target_arg}  ({transcript})")

    results = run(transcript_id=transcript, editor_name=editor_arg, top_n=10)

    if not results:
        print("No safe guides found.")
        sys.exit(0)

    print(f"{'#':<3} {'Protospacer':<22} {'Str':<4} {'GC%':<5} "
          f"{'Eff':<5} {'OT':<5} {'A@pos':<12} Consequence")
    print("─" * 85)

    for i, g in enumerate(results, 1):
        cons_parts = []
        for c in g.consequences:
            tag = f"{c.ref_aa}{c.codon_position}{c.alt_aa}"
            if c.is_synonymous:
                tag += "(syn)"
            cons_parts.append(tag)
        cons_str = ", ".join(cons_parts) if cons_parts else "—"

        print(
            f"{i:<3} {g.protospacer:<22} {g.strand:<4} "
            f"{g.gc_content*100:>3.0f}%  "
            f"{g.efficiency_score:>4.2f}  "
            f"{g.off_target_score:>4.2f}  "
            f"{str(g.edit_window_adenines):<12} {cons_str}"
        )

    print()
    print("Real-model upgrades available via `cbe query` (see README roadmap):")
    print("  efficiency_score: BEDICT-V2 if CBE_BEDICTV2_DIR is set (falls back to heuristic otherwise)")
    print("  off-target verification: `cbe query ... \"check off-target risk for guide N\"` (real NCBI BLAST)")
    print("  protein consequence: ESM-2 via `explain_variant`")
    print("  mRNA payload design: LinearDesign via `design_mrna_payload`")
