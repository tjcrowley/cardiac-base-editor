# Genomic Intake

Data-handling layer in front of `pipeline.py`. `pipeline.py` itself only ever
sees a plain CDS sequence string — it has no idea whether that sequence came
from the public Ensembl reference or from a real person's genome. This module
is what's responsible for handling the latter safely.

## Why this exists

`pipeline.py` was built to rank guide RNAs against a *reference* transcript.
Running it against a real subject's genome means the input is now personal
genomic data — arguably the most sensitive category of personal data there
is — so a few things have to be true before any of that data touches disk:

1. **No data without consent.** Every subject has an explicit consent record
   (`subjects.py`) naming which genes/purposes they've consented to. Nothing
   in this module reads or writes subject data outside that scope.
2. **Everything is attributable.** Every consent grant, revocation, ingest,
   extract, and pipeline run is written to an append-only audit log
   (`audit.py`) — including denied attempts. If someone ever asks "who ran
   what on my genome and when," the answer is in one file.
3. **Data can be permanently deleted.** Revoking consent (`subjects.revoke_consent`)
   purges all stored data for that subject by default (`storage.purge_subject`).
   Right-to-erasure isn't a stretch goal bolted on later — it's the default
   behavior of consent revocation.
4. **At-rest encryption is the operator's box, not this code's job.**
   `DATA_ROOT` is expected to live on an encrypted volume (LUKS/FileVault-style,
   unlocked by a passphrase the operator holds). This module doesn't manage
   keys itself — see the phase-2 note below if that box's threat model
   changes.

## Flow

```
consent grant  →  ingest (FASTQ→VCF, optional)  →  extract (personalized CDS)
                                                        │
                                                        ▼
                                          pipeline.run(sequence=...)  [unchanged]
                                                        │
                                                        ▼
                                                  audit log entry
```

## Usage

```bash
# One-time: grant consent for a subject, scoped to specific genes
python -m genomic_intake consent grant elliot --scope PCSK9,LDLR --days 365

# Run the existing guide-ranking pipeline against that subject's real variants
python -m genomic_intake run elliot PCSK9 --vcf /path/to/elliot.vcf

# Revoke consent — purges all stored data for the subject by default
python -m genomic_intake consent revoke elliot
```

If starting from raw sequencer output instead of a VCF, run `ingest.ingest_fastq()`
first (requires `minimap2`, `samtools`, and `docker` for DeepVariant) to produce
the VCF that `run` expects.

## Web UI

A local, browser-based front end over the same consent/audit/extract/pipeline logic —
no separate service, no new business logic, just friendlier forms and tables:

```bash
pip install -r requirements.txt
uvicorn genomic_intake.web.app:app --reload
```

Open `http://127.0.0.1:8000`. Binds to localhost only by default — this is a
single-operator tool for the box it runs on, not a hosted service. From the dashboard
you can grant/revoke consent per subject, upload a subject's VCF and run it against a
target gene, and review that subject's full audit history in one place.

## Scope limitations (v1)

- **SNVs only.** Indels aren't applied — they'd shift every downstream codon,
  and are lower priority until BE-DICT integration (see `pipeline.py`'s
  roadmap) makes indel-aware scoring meaningful.
- **Single operator, no multi-tenant auth.** This is a CLI for one person
  (Elliot) managing subject records on one box — not a hosted service with
  logins. If this ever needs to run for a *team* of operators, that's a
  deliberate phase-2 scope change, not an oversight.
- **Key management is OS-level, not per-subject.** Currently relies on
  whole-disk encryption on the host. If the box's custody model changes
  (e.g. shared infrastructure, cloud hosting), revisit with per-subject
  encrypted archives instead.

## Testing

See `tests/test_genomic_intake.py` for the consent-gate, variant-application,
and purge-deletion checks. Run with `pytest tests/`.
