import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cardiac_base_editor import plugins


def test_register_tool_adds_to_registry():
    registry_snapshot = dict(plugins.TOOL_REGISTRY)
    try:
        tool = plugins.register_tool(plugins.ExternalTool(
            name="TestTool", env_vars=["CBE_TEST_VAR"],
            setup_instructions="do the thing", check=lambda: True,
        ))
        assert plugins.TOOL_REGISTRY["TestTool"] is tool
    finally:
        plugins.TOOL_REGISTRY.clear()
        plugins.TOOL_REGISTRY.update(registry_snapshot)


def test_require_configured_passes_when_check_true():
    tool = plugins.ExternalTool(name="Ready", env_vars=[], setup_instructions="", check=lambda: True)
    plugins.require_configured(tool)  # should not raise


def test_require_configured_raises_when_check_false():
    tool = plugins.ExternalTool(
        name="NotReady", env_vars=["CBE_X"], setup_instructions="set CBE_X", check=lambda: False,
    )
    with pytest.raises(plugins.ToolNotConfigured) as exc_info:
        plugins.require_configured(tool)

    assert exc_info.value.tool is tool
    assert "NotReady" in str(exc_info.value)
    assert "set CBE_X" in str(exc_info.value)


def test_real_integrations_are_registered():
    """Sanity check that importing the actual model modules registers them —
    catches the case where a future module forgets to call register_tool()."""
    import cardiac_base_editor.mrna_design  # noqa: F401
    import cardiac_base_editor.lnp_delivery  # noqa: F401
    import cardiac_base_editor.cancer.tcr_binding  # noqa: F401
    import cardiac_base_editor.cancer.hla_typing  # noqa: F401

    for name in ("LinearDesign", "LiON", "pMTnet", "arcasHLA"):
        assert name in plugins.TOOL_REGISTRY, f"{name} not registered"
