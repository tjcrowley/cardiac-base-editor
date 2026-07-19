"""
Tests for models/off_target.py. HTTP calls to NCBI's BLAST API are mocked for
deterministic, fast CI tests; one live smoke test exercises the real service
and is skipped if it's unreachable or the request fails (rate limiting,
transient outage, etc. — NCBI's public API has no uptime guarantee).
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cardiac_base_editor.models import off_target as ot

SUBMIT_RESPONSE_TEXT = """
    Some HTML noise
    RID = ABC123XYZ
    RTOE = 15
    more noise
"""

READY_STATUS_TEXT = "QBlastInfoBegin\n    Status=READY\nQBlastInfoEnd"
WAITING_STATUS_TEXT = "QBlastInfoBegin\n    Status=WAITING\nQBlastInfoEnd"

SAMPLE_XML = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_accession>NM_174936</Hit_accession>
          <Hit_def>Homo sapiens PCSK9 mRNA</Hit_def>
          <Hit_hsps>
            <Hsp>
              <Hsp_identity>23</Hsp_identity>
              <Hsp_align-len>23</Hsp_align-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_accession>NM_000001</Hit_accession>
          <Hit_def>Homo sapiens unrelated gene mRNA</Hit_def>
          <Hit_hsps>
            <Hsp>
              <Hsp_identity>20</Hsp_identity>
              <Hsp_align-len>23</Hsp_align-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


def _mock_response(text):
    resp = MagicMock()
    resp.text = text
    resp.raise_for_status = MagicMock()
    return resp


def test_submit_blast_query_parses_rid_and_rtoe():
    with patch.object(ot.requests, "get", return_value=_mock_response(SUBMIT_RESPONSE_TEXT)):
        rid, rtoe = ot.submit_blast_query("ACGTACGTACGT")
    assert rid == "ABC123XYZ"
    assert rtoe == 15


def test_submit_blast_query_raises_without_rid():
    with patch.object(ot.requests, "get", return_value=_mock_response("no rid here")):
        with pytest.raises(ot.BlastError):
            ot.submit_blast_query("ACGTACGTACGT")


def test_poll_blast_status_parses_status():
    with patch.object(ot.requests, "get", return_value=_mock_response(READY_STATUS_TEXT)):
        assert ot.poll_blast_status("ABC123XYZ") == "READY"
    with patch.object(ot.requests, "get", return_value=_mock_response(WAITING_STATUS_TEXT)):
        assert ot.poll_blast_status("ABC123XYZ") == "WAITING"


def test_fetch_blast_hits_parses_xml():
    with patch.object(ot.requests, "get", return_value=_mock_response(SAMPLE_XML)):
        hits = ot.fetch_blast_hits("ABC123XYZ")
    assert len(hits) == 2
    assert hits[0]["accession"] == "NM_174936"
    assert hits[0]["identity_pct"] == 100.0
    assert hits[1]["identity_pct"] == pytest.approx(87.0, abs=0.1)


def test_score_off_target_filters_by_expected_locus(monkeypatch):
    monkeypatch.setattr(ot, "submit_blast_query", lambda seq, **k: ("ABC123XYZ", 0))
    monkeypatch.setattr(ot, "poll_blast_result", lambda rid, rtoe, **k: [
        {"accession": "NM_174936", "definition": "Homo sapiens PCSK9 mRNA", "identity_pct": 100.0, "align_length": 23},
        {"accession": "NM_000001", "definition": "Homo sapiens unrelated gene mRNA", "identity_pct": 87.0, "align_length": 23},
    ])

    result = ot.score_off_target("ACGTACGTACGTACGTACGT", "AGG", expected_locus="PCSK9")
    assert result["total_hits"] == 2
    assert result["off_target_hit_count"] == 1
    assert result["off_target_hits"][0]["accession"] == "NM_000001"


def _ncbi_reachable() -> bool:
    try:
        requests.get(ot.BLAST_URL, timeout=5).raise_for_status()
        return True
    except requests.RequestException:
        return False


@pytest.mark.skipif(not _ncbi_reachable(), reason="NCBI BLAST service not reachable")
def test_live_smoke_submit_real_blast_query():
    """
    Live test against the real NCBI BLAST API: submits a short, distinctive
    sequence and confirms a parseable RID comes back. Does not wait for full
    completion (that can take minutes) — polling/parsing behavior once READY
    is already covered by the mocked tests above.
    """
    rid, rtoe = ot.submit_blast_query("ACGTACGTACGTACGTACGTACG")
    assert rid
    assert rtoe >= 0
