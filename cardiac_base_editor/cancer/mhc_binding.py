"""
MHC-I binding prediction over candidate neoantigen peptides, via mhcflurry
(a real, local, pretrained model — no patient WES or HLA-typing tool needed
here, since alleles are supplied directly by the caller).

Requires `mhcflurry-downloads fetch` to have been run once to pull the
~500MB pretrained weights (documented in the top-level README).
"""

from __future__ import annotations

from functools import lru_cache


@lru_cache(maxsize=1)
def _load_predictor():
    from mhcflurry import Class1PresentationPredictor

    return Class1PresentationPredictor.load()


def rank_binders(peptides: list[str], hla_alleles: list[str]) -> list[dict]:
    """
    Predicts MHC-I binding affinity for every peptide against the given HLA
    alleles, returns results sorted by affinity_percentile ascending
    (tighter/more-likely-presented binders first).
    """
    predictor = _load_predictor()
    df = predictor.predict_affinity(
        peptides=peptides,
        alleles={"subject": hla_alleles},
        verbose=0,
    )
    df = df.sort_values("affinity_percentile")

    return [
        {
            "peptide": row["peptide"],
            "best_allele": row["best_allele"],
            "affinity_nm": float(row["affinity"]),
            "affinity_percentile": float(row["affinity_percentile"]),
        }
        for _, row in df.iterrows()
    ]
