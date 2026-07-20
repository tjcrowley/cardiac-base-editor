"""
Single source of truth for which functions in query/engine.py are callable
as a "query" — by query/nl.py's LLM router, and by web/app.py's generic
structured-query route. Before this existed, query/nl.py hand-maintained a
TOOL_SCHEMA description string and a FUNCTIONS dict separately from
engine.py's actual function list, and they drifted out of sync (the
fallback "I couldn't map that question..." message referenced only 3 of the
6 functions that existed by the time this was written).

Adding a new query function: write it in query/engine.py, decorate it with
@query_function. Nothing else needs updating — query/nl.py's TOOL_SCHEMA and
FUNCTIONS, and web/app.py's structured-query route, all derive from this
registry automatically.
"""

from __future__ import annotations

import inspect
from dataclasses import dataclass
from typing import Callable


@dataclass
class QueryFunctionSpec:
    func: Callable
    hint: str  # one-line routing hint for the LLM


QUERY_FUNCTIONS: dict[str, QueryFunctionSpec] = {}


def query_function(hint: str | None = None):
    """
    hint is the routing description the LLM sees in query/nl.py's
    TOOL_SCHEMA — kept as an explicit decorator argument rather than always
    scraped from the docstring, since some functions need LLM-routing-
    specific guidance (e.g. "slow, only use when explicitly asked") that
    doesn't belong in a docstring written for a Python caller. Defaults to
    the docstring's first line when not given.
    """
    def decorator(func: Callable) -> Callable:
        resolved_hint = hint or (func.__doc__ or "").strip().split("\n")[0]
        QUERY_FUNCTIONS[func.__name__] = QueryFunctionSpec(func=func, hint=resolved_hint)
        return func
    return decorator


def describe_for_llm() -> str:
    """Renders the current registry as the function-list section of
    query/nl.py's TOOL_SCHEMA prompt."""
    lines = []
    for name, spec in QUERY_FUNCTIONS.items():
        params = [
            p for p in inspect.signature(spec.func).parameters
            if p not in ("subject_id", "vcf_path", "operator")
        ]
        lines.append(f"- {name}({', '.join(params)}) -> {spec.hint}")
    return "\n".join(lines)
