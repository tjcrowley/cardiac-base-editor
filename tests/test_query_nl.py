"""
Tests for query/nl.py. Ollama calls are mocked for the routing/parsing/error
tests (deterministic, no live model needed for CI); one live smoke test
exercises the actual local Ollama instance and is skipped if it's unreachable.
"""

import sys
from pathlib import Path

import pytest
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cardiac_base_editor.genomic_intake import audit, storage, subjects
from cardiac_base_editor.query import nl


@pytest.fixture(autouse=True)
def isolated_data_root(tmp_path, monkeypatch):
    monkeypatch.setattr(storage, "DATA_ROOT", tmp_path)
    monkeypatch.setattr(subjects, "REGISTRY_PATH", tmp_path / "subjects.json")
    monkeypatch.setattr(audit, "LOG_PATH", tmp_path / "audit.jsonl")
    yield tmp_path


@pytest.fixture
def synthetic_vcf(tmp_path):
    vcf = tmp_path / "sample.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "1\t1003\t.\tG\tA\n"
    )
    return str(vcf)


def _fake_chat(responses):
    calls = iter(responses)

    def _chat(messages, model=nl.QUERY_MODEL):
        return next(calls)
    return _chat


def test_answer_reports_unreachable_ollama(monkeypatch, synthetic_vcf):
    def _raise(*a, **k):
        raise requests.RequestException("connection refused")
    monkeypatch.setattr(nl, "_ollama_chat", _raise)

    result = nl.answer("elliot", "does elliot have any PCSK9 variants?", synthetic_vcf)
    assert "Couldn't route" in result


def test_answer_reports_unparseable_routing_response(monkeypatch, synthetic_vcf):
    monkeypatch.setattr(nl, "_ollama_chat", _fake_chat(["not json at all"]))

    result = nl.answer("elliot", "what's the weather like?", synthetic_vcf)
    assert "Couldn't route" in result


def test_answer_reports_unsupported_function(monkeypatch, synthetic_vcf):
    monkeypatch.setattr(nl, "_ollama_chat", _fake_chat(['{"function": null, "args": {}}']))

    result = nl.answer("elliot", "what's the weather like?", synthetic_vcf)
    assert "couldn't map that question" in result.lower()


def test_answer_blocked_without_consent(monkeypatch, synthetic_vcf):
    monkeypatch.setattr(
        nl, "_ollama_chat",
        _fake_chat(['{"function": "list_variants", "args": {"gene": "PCSK9"}}']),
    )

    result = nl.answer("elliot", "does elliot have any PCSK9 variants?", synthetic_vcf)
    assert "query failed" in result.lower()


def test_answer_full_flow_with_mocked_model(monkeypatch, synthetic_vcf):
    subjects.grant_consent("elliot", scope=["PCSK9"], operator="tester")

    from cardiac_base_editor.query import engine
    monkeypatch.setattr(engine, "list_variants", lambda **kwargs: [{"genomic_pos": 1003, "ref": "G", "alt": "A"}])
    monkeypatch.setattr(nl, "FUNCTIONS", {"list_variants": engine.list_variants})

    monkeypatch.setattr(
        nl, "_ollama_chat",
        _fake_chat([
            '{"function": "list_variants", "args": {"gene": "PCSK9"}}',
            "Elliot has one variant in PCSK9 at position 1003 (G>A).",
        ]),
    )

    result = nl.answer("elliot", "does elliot have any PCSK9 variants?", synthetic_vcf)
    assert result == "Elliot has one variant in PCSK9 at position 1003 (G>A)."


def _ollama_reachable() -> bool:
    try:
        requests.get(f"{nl.OLLAMA_URL}/api/tags", timeout=2).raise_for_status()
        return True
    except requests.RequestException:
        return False


@pytest.mark.skipif(not _ollama_reachable(), reason="Ollama not reachable on this host")
def test_live_smoke_against_local_ollama(monkeypatch, synthetic_vcf):
    """
    Live end-to-end test against whatever model is actually installed
    locally (llama3.1:8b isn't pulled in this dev environment — override to
    a model that is, e.g. qwen2.5:7b, via CBE_QUERY_MODEL).
    """
    subjects.grant_consent("elliot", scope=["PCSK9"], operator="tester")

    from cardiac_base_editor.query import engine
    monkeypatch.setattr(engine, "list_variants", lambda **kwargs: [{"genomic_pos": 1003, "ref": "G", "alt": "A"}])
    monkeypatch.setattr(nl, "FUNCTIONS", {"list_variants": engine.list_variants})

    result = nl.answer("elliot", "does this subject have any variants in PCSK9?", synthetic_vcf)
    assert "query failed" not in result.lower()
    assert "couldn't" not in result.lower()
