"""
Natural-language front door over query/engine.py, backed by a local Ollama
model — no cloud LLM call, consistent with keeping subject genomic data on
the box.

Privacy invariant: this module never reads subject files directly. It only
ever sees (a) the user's question and (b) whatever small structured dict/list
query/engine.py's functions return. Raw VCF bytes and full subject sequences
never enter an LLM prompt.

Two round trips to Ollama's /api/chat:
  1. "route" call — ask the model to pick one of the three engine functions
     and its arguments, as JSON
  2. "phrase" call — ask the model to turn the structured result into a
     plain-English answer

If Ollama is unreachable, or the model's routing response doesn't parse into
a known function + valid arguments, this returns an explicit "couldn't route
that question" message rather than guessing.
"""

from __future__ import annotations

import inspect
import json
import os

import requests

from cardiac_base_editor.query import engine

OLLAMA_URL = os.environ.get("CBE_OLLAMA_URL", "http://localhost:11434")
QUERY_MODEL = os.environ.get("CBE_QUERY_MODEL", "llama3.1:8b")

TOOL_SCHEMA = """\
You may call exactly one of these functions to answer the user's question about a
genomic subject's data. Respond with ONLY a JSON object: {"function": "<name>", "args": {...}}

- list_variants(gene: str) -> list of the subject's variants found in that gene's coding sequence
- rank_guides(gene: str, editor: str = "ABE8e") -> ranked guide RNA candidates for that gene
- explain_variant(gene: str, genomic_pos: int) -> detailed consequence of one specific variant
- verify_off_target(gene: str, guide_index: int = 0, editor: str = "ABE8e") -> real genome-wide BLAST
  off-target check for one ranked guide (slow, real network call - only use when explicitly asked
  to verify/check off-target risk for a specific guide)
- rank_neoantigens(gene: str, genomic_pos: int, hla_alleles: list[str]) -> ranked candidate neoantigen
  peptides by predicted MHC-I binding, for a somatic variant and a list of HLA alleles (e.g. ["A0201", "B0702"])
- design_mrna_payload(gene: str) -> codon-optimized mRNA sequence for the subject's personalized
  version of that gene (final step: turns the edited/target sequence into an actual mRNA payload)

If the question doesn't clearly map to one of these, respond with:
{"function": null, "args": {}}
"""

FUNCTIONS = {
    "list_variants": engine.list_variants,
    "rank_guides": engine.rank_guides,
    "explain_variant": engine.explain_variant,
    "verify_off_target": engine.verify_off_target,
    "rank_neoantigens": engine.rank_neoantigens,
    "design_mrna_payload": engine.design_mrna_payload,
}


class QueryRoutingError(Exception):
    pass


def _ollama_chat(messages: list[dict], model: str = QUERY_MODEL) -> str:
    resp = requests.post(
        f"{OLLAMA_URL}/api/chat",
        json={"model": model, "messages": messages, "stream": False},
        timeout=30,
    )
    resp.raise_for_status()
    return resp.json()["message"]["content"]


def _route(question: str) -> tuple[str | None, dict]:
    try:
        raw = _ollama_chat([
            {"role": "system", "content": TOOL_SCHEMA},
            {"role": "user", "content": question},
        ])
    except requests.RequestException as e:
        raise QueryRoutingError(f"Could not reach local LLM at {OLLAMA_URL}: {e}") from e

    try:
        json_start = raw.index("{")
        json_end = raw.rindex("}") + 1
        parsed = json.loads(raw[json_start:json_end])
    except (ValueError, json.JSONDecodeError) as e:
        raise QueryRoutingError(f"Couldn't parse routing response: {raw!r}") from e

    return parsed.get("function"), parsed.get("args", {})


def answer(subject_id: str, question: str, vcf_path: str) -> str:
    """
    Answer a free-text question about subject_id's already-ingested genome.
    vcf_path points at the subject's stored VCF, same convention query/engine.py
    functions already use.
    """
    try:
        func_name, args = _route(question)
    except QueryRoutingError as e:
        return f"Couldn't route that question: {e}"

    if func_name not in FUNCTIONS:
        return "I couldn't map that question to a supported query (list variants, rank guides, or explain a specific variant)."

    func = FUNCTIONS[func_name]
    valid_params = set(inspect.signature(func).parameters) - {"subject_id", "vcf_path", "operator"}
    args = {k: v for k, v in args.items() if k in valid_params}  # ignore extra/hallucinated keys

    try:
        result = func(subject_id=subject_id, vcf_path=vcf_path, **args)
    except TypeError as e:
        return f"Couldn't call {func_name} with the arguments the model gave ({args}): {e}"
    except Exception as e:  # consent errors, missing-variant errors, etc. — surface plainly
        return f"That query failed: {e}"

    try:
        phrased = _ollama_chat([
            {"role": "system", "content": "Answer the user's question in plain English using only the JSON data provided. Be concise."},
            {"role": "user", "content": f"Question: {question}\n\nData: {json.dumps(result, default=str)}"},
        ])
        return phrased
    except requests.RequestException:
        # Routing succeeded and we have a real result — degrade to raw JSON
        # rather than losing the answer just because the phrasing call failed.
        return f"(raw result, phrasing call failed) {json.dumps(result, default=str)}"
