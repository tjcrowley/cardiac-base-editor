import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cardiac_base_editor.query import registry


def test_query_function_registers_with_explicit_hint():
    snapshot = dict(registry.QUERY_FUNCTIONS)
    try:
        @registry.query_function(hint="does a thing")
        def my_query(subject_id: str, gene: str, vcf_path: str):
            """Not used as the hint since one was given explicitly."""

        assert registry.QUERY_FUNCTIONS["my_query"].hint == "does a thing"
        assert registry.QUERY_FUNCTIONS["my_query"].func is my_query
    finally:
        registry.QUERY_FUNCTIONS.clear()
        registry.QUERY_FUNCTIONS.update(snapshot)


def test_query_function_defaults_hint_to_docstring_first_line():
    snapshot = dict(registry.QUERY_FUNCTIONS)
    try:
        @registry.query_function()
        def my_other_query(subject_id: str, gene: str, vcf_path: str):
            """First line of the docstring.

            More detail here that shouldn't appear in the hint.
            """

        assert registry.QUERY_FUNCTIONS["my_other_query"].hint == "First line of the docstring."
    finally:
        registry.QUERY_FUNCTIONS.clear()
        registry.QUERY_FUNCTIONS.update(snapshot)


def test_describe_for_llm_lists_registered_functions_excluding_boilerplate_params():
    snapshot = dict(registry.QUERY_FUNCTIONS)
    try:
        registry.QUERY_FUNCTIONS.clear()

        @registry.query_function(hint="ranks things")
        def rank_things(subject_id: str, gene: str, vcf_path: str, top_n: int = 5, operator: str = "x"):
            pass

        description = registry.describe_for_llm()
        assert "rank_things(gene, top_n)" in description
        assert "ranks things" in description
        # boilerplate params shared by every query function shouldn't clutter the LLM prompt
        assert "subject_id" not in description
        assert "vcf_path" not in description
        assert "operator" not in description
    finally:
        registry.QUERY_FUNCTIONS.clear()
        registry.QUERY_FUNCTIONS.update(snapshot)


def test_real_engine_functions_are_registered():
    """Importing query/engine.py should register all 6 @query_function-decorated functions."""
    import cardiac_base_editor.query.engine  # noqa: F401

    for name in (
        "list_variants", "rank_guides", "explain_variant",
        "verify_off_target", "rank_neoantigens", "design_mrna_payload",
    ):
        assert name in registry.QUERY_FUNCTIONS, f"{name} not registered"
