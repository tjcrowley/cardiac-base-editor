"""
Shared configuration pattern for external model integrations that have no
fallback (mrna_design.py, lnp_delivery.py, cancer/tcr_binding.py,
cancer/hla_typing.py) — each needs one or more env vars pointing at a real
external tool/checkout, and each used to hand-roll its own exception class +
env-var-check function + setup-instructions string to express that. This
module is that pattern, written once.

models/be_dict.py deliberately does NOT use this: it always has a working
fallback (the original heuristic), so "not configured" there is a warning +
degrade, not a hard failure — a genuinely different design, not an
inconsistency to paper over by forcing it in here.

Adding a new hard-fail external-model integration:
    from cardiac_base_editor.plugins import ExternalTool, register_tool, require_configured

    TOOL = register_tool(ExternalTool(
        name="MyModel",
        env_vars=["CBE_MYMODEL_DIR"],
        setup_instructions="git clone ...\\n  export CBE_MYMODEL_DIR=...",
        check=lambda: bool(os.environ.get("CBE_MYMODEL_DIR")),
    ))

    def my_function(...):
        require_configured(TOOL)
        ...

`cbe plugins` (cli.py) lists every registered tool's configured/not status.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable


@dataclass
class ExternalTool:
    name: str
    env_vars: list[str]
    setup_instructions: str
    check: Callable[[], bool]


class ToolNotConfigured(Exception):
    def __init__(self, tool: ExternalTool):
        self.tool = tool
        super().__init__(f"{tool.name} not configured. {tool.setup_instructions}")


TOOL_REGISTRY: dict[str, ExternalTool] = {}


def register_tool(tool: ExternalTool) -> ExternalTool:
    TOOL_REGISTRY[tool.name] = tool
    return tool


def require_configured(tool: ExternalTool) -> None:
    if not tool.check():
        raise ToolNotConfigured(tool)
