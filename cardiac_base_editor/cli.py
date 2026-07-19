"""
Merged CLI entry point: consent management, guide-ranking pipeline, and the
subject genome query engine.

Usage:
    cbe consent grant <subject_id> --scope PCSK9,LDLR [--days 365]
    cbe consent revoke <subject_id>
    cbe run <subject_id> <gene> --vcf path/to/subject.vcf [--editor ABE8e]
    cbe query <subject_id> "<question>" --vcf path/to/subject.vcf

Every `run`/`query` invocation checks consent first (subjects.require_consent).
`run` builds the subject's personalized CDS (extract.py) and hands it to the
existing, unmodified pipeline.run() from pipeline.py. `query` routes a
free-text question through query/nl.py's local-LLM front door, which itself
only ever calls the same consent-gated functions in query/engine.py.
"""

from __future__ import annotations

import argparse
import getpass
import sys

from cardiac_base_editor.genomic_intake import subjects
from cardiac_base_editor.genomic_intake.extract import build_personalized_cds
from cardiac_base_editor.pipeline import KNOWN_TARGETS, run as pipeline_run
from cardiac_base_editor.query import nl as query_nl


def cmd_consent_grant(args: argparse.Namespace) -> None:
    scope = args.scope.split(",") if args.scope else ["*"]
    subjects.grant_consent(
        subject_id=args.subject_id,
        scope=scope,
        retention_days=args.days,
        operator=getpass.getuser(),
    )
    print(f"Consent granted for '{args.subject_id}', scope={scope}, retention={args.days}d")


def cmd_consent_revoke(args: argparse.Namespace) -> None:
    subjects.revoke_consent(args.subject_id, operator=getpass.getuser(), purge_data=not args.keep_data)
    print(f"Consent revoked for '{args.subject_id}'" + ("" if args.keep_data else " and data purged"))


def cmd_run(args: argparse.Namespace) -> None:
    operator = getpass.getuser()
    try:
        subjects.require_consent(args.subject_id, args.gene, operator=operator)
    except subjects.ConsentError as e:
        print(f"BLOCKED: {e}", file=sys.stderr)
        sys.exit(1)

    transcript_id = KNOWN_TARGETS.get(args.gene.upper(), args.gene)
    personalized_cds = build_personalized_cds(transcript_id, args.vcf)

    results = pipeline_run(sequence=personalized_cds, editor_name=args.editor, top_n=args.top_n)
    if not results:
        print("No safe guides found.")
        return

    print(f"\nTop guide candidates for subject '{args.subject_id}', gene {args.gene}:\n")
    for i, g in enumerate(results, 1):
        print(f"{i}. {g.protospacer} ({g.strand}) eff={g.efficiency_score:.2f} ot={g.off_target_score:.2f}")


def cmd_query(args: argparse.Namespace) -> None:
    print(query_nl.answer(args.subject_id, args.question, args.vcf))


def main() -> None:
    parser = argparse.ArgumentParser(prog="cbe")
    sub = parser.add_subparsers(dest="command", required=True)

    consent = sub.add_parser("consent")
    consent_sub = consent.add_subparsers(dest="consent_action", required=True)

    grant = consent_sub.add_parser("grant")
    grant.add_argument("subject_id")
    grant.add_argument("--scope", default="*", help="Comma-separated gene list, or '*' for all")
    grant.add_argument("--days", type=int, default=365)
    grant.set_defaults(func=cmd_consent_grant)

    revoke = consent_sub.add_parser("revoke")
    revoke.add_argument("subject_id")
    revoke.add_argument("--keep-data", action="store_true", help="Revoke consent without purging stored data")
    revoke.set_defaults(func=cmd_consent_revoke)

    run_p = sub.add_parser("run")
    run_p.add_argument("subject_id")
    run_p.add_argument("gene")
    run_p.add_argument("--vcf", required=True)
    run_p.add_argument("--editor", default="ABE8e")
    run_p.add_argument("--top-n", type=int, default=10, dest="top_n")
    run_p.set_defaults(func=cmd_run)

    query_p = sub.add_parser("query")
    query_p.add_argument("subject_id")
    query_p.add_argument("question")
    query_p.add_argument("--vcf", required=True)
    query_p.set_defaults(func=cmd_query)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
