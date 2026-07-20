"""
Structural tests only for cancer/tcr_binding.py (pMTnet wrapper) — no live
model inference, per the module's own docstring: this environment has no
legacy TF1/Keras 2.2.4/numpy 1.16.3-compatible Python available. These tests
confirm the wrapper builds the correct input format and command, and fails
clearly when unconfigured.
"""

import csv
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cardiac_base_editor.cancer import tcr_binding
from cardiac_base_editor.plugins import ToolNotConfigured


def test_raises_clear_error_when_not_configured(monkeypatch):
    monkeypatch.delenv("CBE_PMTNET_DIR", raising=False)
    with pytest.raises(ToolNotConfigured):
        tcr_binding.predict_tcr_binding([{"cdr3": "CASSVASSGNIQYF", "antigen": "TQPPSGFR", "hla": "A*11:01"}])


def test_raises_when_configured_dir_missing_script(monkeypatch, tmp_path):
    monkeypatch.setenv("CBE_PMTNET_DIR", str(tmp_path))  # exists but no pMTnet.py inside
    with pytest.raises(ToolNotConfigured):
        tcr_binding.predict_tcr_binding([{"cdr3": "X", "antigen": "Y", "hla": "Z"}])


def test_writes_correct_input_csv_and_reads_output(monkeypatch, tmp_path):
    pmtnet_dir = tmp_path / "pMTnet"
    pmtnet_dir.mkdir()
    (pmtnet_dir / "pMTnet.py").write_text("# stub")
    monkeypatch.setenv("CBE_PMTNET_DIR", str(pmtnet_dir))

    captured_input_csv = {}

    def fake_run(cmd, check, capture_output, text, timeout):
        # Locate the -input arg to inspect what was actually written
        input_csv_path = cmd[cmd.index("-input") + 1]
        with open(input_csv_path) as f:
            captured_input_csv["rows"] = list(csv.DictReader(f))

        output_dir = cmd[cmd.index("-output") + 1]
        with open(f"{output_dir}/prediction.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["CDR3", "Antigen", "HLA", "Rank"])
            writer.writerow(["CASSVASSGNIQYF", "TQPPSGFR", "A*11:01", "0.05"])
        return MagicMock(returncode=0)

    with patch.object(tcr_binding.subprocess, "run", side_effect=fake_run):
        results = tcr_binding.predict_tcr_binding([
            {"cdr3": "CASSVASSGNIQYF", "antigen": "TQPPSGFR", "hla": "A*11:01"},
        ])

    assert captured_input_csv["rows"] == [{"CDR3": "CASSVASSGNIQYF", "Antigen": "TQPPSGFR", "HLA": "A*11:01"}]
    assert results == [{"CDR3": "CASSVASSGNIQYF", "Antigen": "TQPPSGFR", "HLA": "A*11:01", "Rank": "0.05"}]
