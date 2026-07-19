"""
Tests for query/engine.py: consent gating and correct structured output.
No network calls — Ensembl-dependent functions are monkeypatched with
synthetic data, same approach as test_genomic_intake.py.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cardiac_base_editor.genomic_intake import audit, storage, subjects
from cardiac_base_editor.genomic_intake import extract
from cardiac_base_editor.query import engine


@pytest.fixture(autouse=True)
def isolated_data_root(tmp_path, monkeypatch):
    monkeypatch.setattr(storage, "DATA_ROOT", tmp_path)
    monkeypatch.setattr(subjects, "REGISTRY_PATH", tmp_path / "subjects.json")
    monkeypatch.setattr(audit, "LOG_PATH", tmp_path / "audit.jsonl")
    yield tmp_path


@pytest.fixture
def synthetic_transcript(monkeypatch):
    """A 30bp toy CDS on the plus strand, genomic pos 1000+i <-> cds pos i+1."""
    reference_cds = "ATGGGGCCCAAATTTTAGATGCCCGGGTTT"
    monkeypatch.setattr(extract, "fetch_cds", lambda transcript_id: reference_cds)
    monkeypatch.setattr(engine, "fetch_cds", lambda transcript_id: reference_cds)
    monkeypatch.setattr(
        extract, "_fetch_cds_genomic_map",
        lambda t, n: ({1000 + i: i + 1 for i in range(n)}, 1),
    )
    monkeypatch.setattr(
        engine, "_fetch_cds_genomic_map",
        lambda t, n: ({1000 + i: i + 1 for i in range(n)}, 1),
    )
    return reference_cds


@pytest.fixture
def synthetic_vcf(tmp_path):
    vcf = tmp_path / "sample.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "1\t1003\t.\tG\tA\n"
    )
    return str(vcf)


def test_list_variants_blocked_without_consent(synthetic_transcript, synthetic_vcf):
    with pytest.raises(subjects.ConsentError):
        engine.list_variants("elliot", "PCSK9", synthetic_vcf)


def test_list_variants_returns_variants_in_transcript(synthetic_transcript, synthetic_vcf):
    subjects.grant_consent("elliot", scope=["PCSK9"], operator="tester")
    hits = engine.list_variants("elliot", "PCSK9", synthetic_vcf)
    assert len(hits) == 1
    assert hits[0]["genomic_pos"] == 1003
    assert hits[0]["cds_position"] == 4


def test_rank_guides_blocked_without_consent(synthetic_transcript, synthetic_vcf):
    with pytest.raises(subjects.ConsentError):
        engine.rank_guides("elliot", "PCSK9", synthetic_vcf)


def test_rank_guides_returns_guide_dicts(synthetic_transcript, synthetic_vcf):
    subjects.grant_consent("elliot", scope=["PCSK9"], operator="tester")
    guides = engine.rank_guides("elliot", "PCSK9", synthetic_vcf)
    assert isinstance(guides, list)
    if guides:
        assert "protospacer" in guides[0]
        assert "efficiency_score" in guides[0]


def test_explain_variant_blocked_without_consent(synthetic_transcript, synthetic_vcf):
    with pytest.raises(subjects.ConsentError):
        engine.explain_variant("elliot", "PCSK9", synthetic_vcf, genomic_pos=1003)


def test_explain_variant_returns_consequence(monkeypatch, synthetic_transcript, synthetic_vcf):
    subjects.grant_consent("elliot", scope=["PCSK9"], operator="tester")
    # Avoid downloading ESM-2 in this test — that's covered separately in
    # test_protein_consequence.py.
    monkeypatch.setattr(engine, "score_substitution", lambda *a, **k: -1.23)

    result = engine.explain_variant("elliot", "PCSK9", synthetic_vcf, genomic_pos=1003)
    assert result["genomic_pos"] == 1003
    assert result["ref"] == "G" and result["alt"] == "A"
    assert result["codon_position"] == 2  # cds pos 4 -> codon index 1 -> position 2
    assert result["esm2_disruption_score"] == -1.23


def test_explain_variant_raises_for_missing_variant(synthetic_transcript, synthetic_vcf):
    subjects.grant_consent("elliot", scope=["PCSK9"], operator="tester")
    with pytest.raises(ValueError):
        engine.explain_variant("elliot", "PCSK9", synthetic_vcf, genomic_pos=9999)
