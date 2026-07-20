"""
Structural tests only for cancer/hla_typing.py (arcasHLA wrapper) — no live
run against real RNA-seq data, per the module's own docstring: there's no
synthetic BAM that would produce a meaningful HLA call. These tests confirm
the wrapper resolves the binary correctly and fails clearly when unconfigured.
"""

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cardiac_base_editor.cancer import hla_typing
from cardiac_base_editor.plugins import ToolNotConfigured


def test_raises_clear_error_when_not_on_path_and_no_dir_given(monkeypatch):
    monkeypatch.setattr(hla_typing, "which", lambda name: None)
    with pytest.raises(ToolNotConfigured):
        hla_typing.type_hla_from_rnaseq("/fake/path.bam")


def test_raises_when_arcashla_dir_missing_binary(tmp_path):
    with pytest.raises(ToolNotConfigured):
        hla_typing.type_hla_from_rnaseq("/fake/path.bam", arcashla_dir=str(tmp_path))


def test_uses_binary_from_arcashla_dir_and_parses_genotype_json(tmp_path):
    arcashla_dir = tmp_path / "arcasHLA"
    arcashla_dir.mkdir()
    binary = arcashla_dir / "arcasHLA"
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o755)

    def fake_run(cmd, check, capture_output, text, timeout):
        output_dir = cmd[cmd.index("-o") + 1]
        with open(f"{output_dir}/sample.genotype.json", "w") as f:
            json.dump({"A": ["A*02:01", "A*11:01"], "B": ["B*07:02"]}, f)
        return MagicMock(returncode=0)

    with patch.object(hla_typing.subprocess, "run", side_effect=fake_run):
        alleles = hla_typing.type_hla_from_rnaseq("/fake/path.bam", arcashla_dir=str(arcashla_dir))

    assert set(alleles) == {"A*02:01", "A*11:01", "B*07:02"}
