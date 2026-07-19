"""
Tests for genomic_intake: consent gate, variant application, purge/erasure,
and audit logging. No network calls — Ensembl-dependent functions in
extract.py are monkeypatched with synthetic data.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from genomic_intake import audit, storage, subjects
from genomic_intake import extract


@pytest.fixture(autouse=True)
def isolated_data_root(tmp_path, monkeypatch):
    """Point all storage/registry/audit paths at a fresh tmp dir per test."""
    monkeypatch.setattr(storage, "DATA_ROOT", tmp_path)
    monkeypatch.setattr(subjects, "REGISTRY_PATH", tmp_path / "subjects.json")
    monkeypatch.setattr(audit, "LOG_PATH", tmp_path / "audit.jsonl")
    yield tmp_path


# ── Consent gate ──────────────────────────────────────────────────────────

def test_no_consent_blocks_access():
    with pytest.raises(subjects.ConsentError):
        subjects.require_consent("elliot", "PCSK9", operator="tester")


def test_consent_grant_allows_scoped_gene():
    subjects.grant_consent("elliot", scope=["PCSK9"], operator="tester")
    record = subjects.require_consent("elliot", "PCSK9", operator="tester")
    assert record.subject_id == "elliot"


def test_consent_does_not_cover_out_of_scope_gene():
    subjects.grant_consent("elliot", scope=["PCSK9"], operator="tester")
    with pytest.raises(subjects.ConsentError):
        subjects.require_consent("elliot", "LDLR", operator="tester")


def test_wildcard_scope_covers_any_gene():
    subjects.grant_consent("elliot", scope=["*"], operator="tester")
    subjects.require_consent("elliot", "ANGPTL3", operator="tester")  # should not raise


def test_revoked_consent_blocks_access():
    subjects.grant_consent("elliot", scope=["*"], operator="tester")
    subjects.revoke_consent("elliot", operator="tester")
    with pytest.raises(subjects.ConsentError):
        subjects.require_consent("elliot", "PCSK9", operator="tester")


# ── Purge / right to erasure ──────────────────────────────────────────────

def test_revoke_purges_subject_data(isolated_data_root):
    subjects.grant_consent("elliot", scope=["*"], operator="tester")
    subject_path = storage.subject_dir("elliot")
    (subject_path / "raw" / "sample.vcf").write_text("fake vcf data")
    assert subject_path.exists()

    subjects.revoke_consent("elliot", operator="tester", purge_data=True)

    assert not subject_path.exists()


def test_revoke_with_keep_data_preserves_files(isolated_data_root):
    subjects.grant_consent("elliot", scope=["*"], operator="tester")
    subject_path = storage.subject_dir("elliot")
    (subject_path / "raw" / "sample.vcf").write_text("fake vcf data")

    subjects.revoke_consent("elliot", operator="tester", purge_data=False)

    assert subject_path.exists()


# ── Audit log ──────────────────────────────────────────────────────────────

def test_audit_log_records_denied_access():
    with pytest.raises(subjects.ConsentError):
        subjects.require_consent("elliot", "PCSK9", operator="tester")

    entries = audit.read_for_subject("elliot")
    assert any(e["action"] == "access_denied" and e["allowed"] is False for e in entries)


def test_audit_log_records_granted_access():
    subjects.grant_consent("elliot", scope=["*"], operator="tester")
    subjects.require_consent("elliot", "PCSK9", operator="tester")

    entries = audit.read_for_subject("elliot")
    actions = [e["action"] for e in entries]
    assert "consent_granted" in actions
    assert "access_granted" in actions


# ── Variant application (extract.py), no network ─────────────────────────

def test_vcf_parser_skips_indels_and_multiallelic(tmp_path):
    vcf = tmp_path / "sample.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "1\t100\t.\tA\tG\n"       # valid SNV
        "1\t200\t.\tAT\tA\n"      # indel, skip
        "1\t300\t.\tC\tG,T\n"     # multi-allelic, skip
    )
    variants = extract.parse_vcf_snvs(str(vcf))
    assert len(variants) == 1
    assert variants[0].pos == 100 and variants[0].ref == "A" and variants[0].alt == "G"


def test_build_personalized_cds_applies_snv(monkeypatch, tmp_path):
    reference_cds = "ATGGGGCCCAAATTTTAG"  # 18bp toy CDS, plus strand

    monkeypatch.setattr(extract, "fetch_cds", lambda transcript_id: reference_cds)

    def fake_map(transcript_id, cds_length):
        # genomic pos 1000+i maps 1:1 to cds pos i+1, plus strand
        return {1000 + i: i + 1 for i in range(cds_length)}, 1

    monkeypatch.setattr(extract, "_fetch_cds_genomic_map", fake_map)

    vcf = tmp_path / "sample.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "1\t1003\t.\tG\tA\n"   # cds position 4 (1-indexed): reference_cds[3] == 'G' -> 'A'
    )

    personalized = extract.build_personalized_cds("ENST_FAKE", str(vcf))
    expected = list(reference_cds)
    expected[3] = "A"
    assert personalized == "".join(expected)


def test_build_personalized_cds_ignores_variant_outside_transcript(monkeypatch, tmp_path):
    reference_cds = "ATGGGGCCCAAATTTTAG"
    monkeypatch.setattr(extract, "fetch_cds", lambda transcript_id: reference_cds)
    monkeypatch.setattr(extract, "_fetch_cds_genomic_map", lambda t, n: ({1000 + i: i + 1 for i in range(n)}, 1))

    vcf = tmp_path / "sample.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "1\t9999\t.\tG\tA\n"  # not in the mapped range
    )

    personalized = extract.build_personalized_cds("ENST_FAKE", str(vcf))
    assert personalized == reference_cds
