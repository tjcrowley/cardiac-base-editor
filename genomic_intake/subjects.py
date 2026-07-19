"""
Subject registry + consent gate.

A "subject" is one genome donor (Elliot, or anyone else whose genome is
processed). Nothing in genomic_intake touches a subject's data without an
active consent record. This is a flat JSON registry — appropriate for one
operator managing a handful of subjects on a single box, not a multi-tenant
service with logins.

Registry file: DATA_ROOT/subjects.json
{
  "<subject_id>": {
    "granted_at": "2026-07-19T...",
    "revoked_at": null,
    "scope": ["PCSK9", "LDLR"],   # genes/purposes consented to, or ["*"] for all
    "retention_days": 365
  }
}
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path

from genomic_intake import audit
from genomic_intake.storage import DATA_ROOT, init_subject_storage, purge_subject

REGISTRY_PATH = DATA_ROOT / "subjects.json"


class ConsentError(Exception):
    pass


@dataclass
class ConsentRecord:
    subject_id: str
    granted_at: str
    revoked_at: str | None
    scope: list[str]
    retention_days: int

    @property
    def active(self) -> bool:
        return self.revoked_at is None

    def covers(self, gene: str) -> bool:
        return self.active and ("*" in self.scope or gene.upper() in [s.upper() for s in self.scope])


def _load_registry() -> dict:
    if not REGISTRY_PATH.exists():
        return {}
    with open(REGISTRY_PATH) as f:
        return json.load(f)


def _save_registry(registry: dict) -> None:
    REGISTRY_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(REGISTRY_PATH, "w") as f:
        json.dump(registry, f, indent=2)


def grant_consent(subject_id: str, scope: list[str], retention_days: int = 365, operator: str = "unknown") -> ConsentRecord:
    registry = _load_registry()
    record = {
        "granted_at": datetime.now(timezone.utc).isoformat(),
        "revoked_at": None,
        "scope": scope,
        "retention_days": retention_days,
    }
    registry[subject_id] = record
    _save_registry(registry)
    init_subject_storage(subject_id)
    audit.record(subject_id, "consent_granted", operator, detail=f"scope={scope}")
    return ConsentRecord(subject_id=subject_id, **record)


def revoke_consent(subject_id: str, operator: str = "unknown", purge_data: bool = True) -> None:
    registry = _load_registry()
    if subject_id not in registry:
        raise ConsentError(f"No consent record for subject '{subject_id}'")
    registry[subject_id]["revoked_at"] = datetime.now(timezone.utc).isoformat()
    _save_registry(registry)
    audit.record(subject_id, "consent_revoked", operator)
    if purge_data:
        purged = purge_subject(subject_id)
        audit.record(subject_id, "data_purged", operator, detail=f"purged={purged}")


def get_consent(subject_id: str) -> ConsentRecord | None:
    registry = _load_registry()
    if subject_id not in registry:
        return None
    return ConsentRecord(subject_id=subject_id, **registry[subject_id])


def list_all() -> dict[str, ConsentRecord]:
    registry = _load_registry()
    return {
        subject_id: ConsentRecord(subject_id=subject_id, **fields)
        for subject_id, fields in registry.items()
    }


def require_consent(subject_id: str, gene: str, operator: str = "unknown") -> ConsentRecord:
    """
    Raise ConsentError if the subject has no active consent covering this
    gene. Every call — allowed or denied — is written to the audit log.
    """
    record = get_consent(subject_id)
    if record is None or not record.covers(gene):
        audit.record(subject_id, "access_denied", operator, detail=f"gene={gene}", allowed=False)
        raise ConsentError(f"No active consent for subject '{subject_id}' covering gene '{gene}'")
    audit.record(subject_id, "access_granted", operator, detail=f"gene={gene}")
    return record


def retention_expired(record: ConsentRecord) -> bool:
    granted = datetime.fromisoformat(record.granted_at)
    return datetime.now(timezone.utc) > granted + timedelta(days=record.retention_days)
