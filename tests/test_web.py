"""
Tests for the genomic_intake web UI. Uses FastAPI's TestClient, no running
server process required. Ensembl-dependent calls (fetch_cds, coordinate
mapping) are monkeypatched — same synthetic-CDS approach as test_genomic_intake.py.
"""

import io
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cardiac_base_editor.genomic_intake import audit, extract, storage, subjects
from cardiac_base_editor.web import app as web_app


@pytest.fixture(autouse=True)
def isolated_data_root(tmp_path, monkeypatch):
    monkeypatch.setattr(storage, "DATA_ROOT", tmp_path)
    monkeypatch.setattr(subjects, "REGISTRY_PATH", tmp_path / "subjects.json")
    monkeypatch.setattr(audit, "LOG_PATH", tmp_path / "audit.jsonl")
    yield tmp_path


@pytest.fixture
def client():
    return TestClient(web_app.app)


def test_dashboard_loads_empty(client):
    resp = client.get("/")
    assert resp.status_code == 200
    assert "No subjects yet" in resp.text


def test_create_subject_via_form(client):
    resp = client.post("/subjects", data={"subject_id": "elliot", "scope": "PCSK9", "retention_days": "365"})
    assert resp.status_code == 200
    assert "elliot" in resp.text
    assert subjects.get_consent("elliot") is not None


def test_subject_detail_shows_consent(client):
    client.post("/subjects", data={"subject_id": "elliot", "scope": "PCSK9", "retention_days": "365"})
    resp = client.get("/subjects/elliot")
    assert resp.status_code == 200
    assert "PCSK9" in resp.text
    assert "active" in resp.text.lower()


def test_revoke_purges_data_via_web(client):
    client.post("/subjects", data={"subject_id": "elliot", "scope": "*", "retention_days": "365"})
    subject_path = storage.subject_dir("elliot")
    (subject_path / "raw" / "sample.vcf").write_text("fake")
    assert subject_path.exists()

    resp = client.post(f"/subjects/elliot/revoke", data={"purge_data": "true"})
    assert resp.status_code == 200
    assert not subject_path.exists()
    assert "revoked" in resp.text.lower()


def test_run_blocked_without_consent(client, tmp_path):
    vcf_bytes = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n1\t100\t.\tG\tA\n"
    resp = client.post(
        "/subjects/elliot/run",
        data={"gene": "PCSK9", "editor": "ABE8e"},
        files={"vcf": ("sample.vcf", io.BytesIO(vcf_bytes), "text/plain")},
    )
    assert resp.status_code == 200
    assert "Blocked" in resp.text


def test_run_with_consent_returns_ranked_guides(client, monkeypatch):
    reference_cds = "ATGGGGCCCAAATTTTAG" * 5  # toy CDS with a PAM-friendly region
    monkeypatch.setattr(extract, "fetch_cds", lambda transcript_id: reference_cds)
    monkeypatch.setattr(
        extract, "_fetch_cds_genomic_map",
        lambda t, n: ({1000 + i: i + 1 for i in range(n)}, 1),
    )

    client.post("/subjects", data={"subject_id": "elliot", "scope": "PCSK9", "retention_days": "365"})

    vcf_bytes = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n1\t1003\t.\tG\tA\n"
    resp = client.post(
        "/subjects/elliot/run",
        data={"gene": "PCSK9", "editor": "ABE8e"},
        files={"vcf": ("sample.vcf", io.BytesIO(vcf_bytes), "text/plain")},
    )
    assert resp.status_code == 200
    assert "Blocked" not in resp.text

    entries = audit.read_for_subject("elliot")
    assert any(e["action"] == "access_granted" for e in entries)
