"""
Append-only audit log.

Every consent check, ingest, extract, and pipeline run gets a JSONL entry —
including rejected attempts (no consent on file). The log is append-only by
convention: nothing in this module ever rewrites or deletes existing lines.
purge_subject() in storage.py deliberately leaves this log untouched.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from genomic_intake.storage import DATA_ROOT

LOG_PATH = DATA_ROOT / "audit.jsonl"


def record(subject_id: str, action: str, operator: str, detail: str = "", allowed: bool = True) -> None:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "subject_id": subject_id,
        "operator": operator,
        "action": action,
        "allowed": allowed,
        "detail": detail,
    }
    with open(LOG_PATH, "a") as f:
        f.write(json.dumps(entry) + "\n")


def read_all() -> list[dict]:
    if not LOG_PATH.exists():
        return []
    with open(LOG_PATH) as f:
        return [json.loads(line) for line in f if line.strip()]


def read_for_subject(subject_id: str) -> list[dict]:
    return [e for e in read_all() if e["subject_id"] == subject_id]
