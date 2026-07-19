"""
Subject-scoped storage layout.

All subject genomic data lives under DATA_ROOT/subjects/<subject_id>/. This
module assumes DATA_ROOT itself sits on an encrypted volume (e.g. a LUKS
container on the inference box) — it does not perform its own encryption.

Layout per subject:
    subjects/<subject_id>/
        raw/            # original FASTQ/VCF as received
        derived/        # alignment/variant-calling intermediates (ingest.py output)
        RETENTION.json   # {"delete_by": ISO date or null}
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

DATA_ROOT = Path(os.environ.get("GENOMIC_INTAKE_DATA_ROOT", "./data"))


def subject_dir(subject_id: str) -> Path:
    return DATA_ROOT / "subjects" / subject_id


def raw_dir(subject_id: str) -> Path:
    return subject_dir(subject_id) / "raw"


def derived_dir(subject_id: str) -> Path:
    return subject_dir(subject_id) / "derived"


def init_subject_storage(subject_id: str) -> None:
    raw_dir(subject_id).mkdir(parents=True, exist_ok=True)
    derived_dir(subject_id).mkdir(parents=True, exist_ok=True)


def purge_subject(subject_id: str) -> bool:
    """
    Right-to-erasure: permanently delete all stored data for a subject.

    Returns True if a directory existed and was removed, False if there was
    nothing to purge. Does not touch the audit log or the subject registry
    entry — callers (subjects.py) are responsible for marking the subject
    record as erased.
    """
    d = subject_dir(subject_id)
    if not d.exists():
        return False
    shutil.rmtree(d)
    return True
