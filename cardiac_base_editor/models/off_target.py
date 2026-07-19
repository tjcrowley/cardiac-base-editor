"""
Genome-wide off-target scoring via NCBI's public BLAST Common URL API
(https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html) — no local genome
download needed.

This is deliberately NOT wired into pipeline.score_guides()'s default path:
BLAST jobs against the human genome routinely take a minute or more to
complete, and NCBI's usage policy caps request volume. It's an opt-in
verification step for one guide at a time (query/engine.py's
verify_off_target), used after the fast local heuristic has already ranked
candidates — not a bulk replacement for it.

NCBI usage policy requires identifying email/tool parameters and asks
callers not to poll more than once every ~10s or submit >100 searches/24h.
"""

from __future__ import annotations

import os
import re
import time
import xml.etree.ElementTree as ET

import requests

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
CONTACT_EMAIL = os.environ.get("CBE_CONTACT_EMAIL", "anonymous@example.com")
TOOL_NAME = "cardiac-base-editor"

# Human-restricted nucleotide search via ENTREZ_QUERY rather than a specific
# internal genome database name, which NCBI doesn't guarantee stays stable —
# "nt" + organism filter is the documented, portable way to search "the human
# genome" through this API.
DEFAULT_DATABASE = "nt"
DEFAULT_ENTREZ_QUERY = "Homo sapiens[Organism]"

MIN_POLL_INTERVAL_S = 10


class BlastError(Exception):
    pass


def submit_blast_query(sequence: str, database: str = DEFAULT_DATABASE, entrez_query: str = DEFAULT_ENTREZ_QUERY) -> tuple[str, int]:
    """Submit a short-sequence megablast search. Returns (RID, RTOE_seconds)."""
    resp = requests.get(BLAST_URL, params={
        "CMD": "Put",
        "PROGRAM": "blastn",
        "MEGABLAST": "on",
        "DATABASE": database,
        "QUERY": sequence,
        "ENTREZ_QUERY": entrez_query,
        "HITLIST_SIZE": 20,
        "email": CONTACT_EMAIL,
        "tool": TOOL_NAME,
    }, timeout=30)
    resp.raise_for_status()

    rid_match = re.search(r"RID = (\S+)", resp.text)
    rtoe_match = re.search(r"RTOE = (\d+)", resp.text)
    if not rid_match:
        raise BlastError(f"No RID returned by BLAST submission: {resp.text[:500]}")

    return rid_match.group(1), int(rtoe_match.group(1)) if rtoe_match else 20


def poll_blast_status(rid: str) -> str:
    """Returns 'WAITING', 'READY', 'FAILED', or 'UNKNOWN'."""
    resp = requests.get(BLAST_URL, params={
        "CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid,
    }, timeout=30)
    resp.raise_for_status()
    match = re.search(r"Status=(\w+)", resp.text)
    return match.group(1) if match else "UNKNOWN"


def fetch_blast_hits(rid: str) -> list[dict]:
    """Fetch and parse XML results for a completed search."""
    resp = requests.get(BLAST_URL, params={
        "CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid,
    }, timeout=60)
    resp.raise_for_status()

    root = ET.fromstring(resp.text)
    hits = []
    for hit in root.iter("Hit"):
        hit_def = hit.findtext("Hit_def", default="")
        hit_accession = hit.findtext("Hit_accession", default="")
        for hsp in hit.iter("Hsp"):
            identity = int(hsp.findtext("Hsp_identity", default="0"))
            align_len = int(hsp.findtext("Hsp_align-len", default="1"))
            hits.append({
                "accession": hit_accession,
                "definition": hit_def,
                "identity_pct": round(100 * identity / align_len, 1),
                "align_length": align_len,
            })
    return hits


def poll_blast_result(rid: str, rtoe: int = 20, max_wait_s: int = 300) -> list[dict]:
    """Blocks until the search is READY (or raises on FAILED/timeout), then
    returns parsed hits. Honors NCBI's suggested minimum poll interval."""
    time.sleep(min(rtoe, max_wait_s))
    waited = rtoe
    while waited < max_wait_s:
        status = poll_blast_status(rid)
        if status == "READY":
            return fetch_blast_hits(rid)
        if status == "FAILED":
            raise BlastError(f"BLAST search {rid} failed")
        if status == "UNKNOWN":
            raise BlastError(f"BLAST search {rid} expired or RID unknown")
        time.sleep(MIN_POLL_INTERVAL_S)
        waited += MIN_POLL_INTERVAL_S

    raise BlastError(f"BLAST search {rid} did not complete within {max_wait_s}s")


def score_off_target(protospacer: str, pam_seq: str, expected_locus: str | None = None) -> dict:
    """
    Submits protospacer+PAM as a single query, waits for genome-wide BLAST
    results, and returns a summary: total hit count, and hits that don't
    match expected_locus (a substring to match against hit definitions/
    accessions, e.g. a gene symbol or RefSeq ID) as candidate off-target
    sites. More non-matching hits with high identity = higher off-target risk.
    """
    query_seq = protospacer + pam_seq
    rid, rtoe = submit_blast_query(query_seq)
    hits = poll_blast_result(rid, rtoe)

    off_target_hits = hits
    if expected_locus:
        off_target_hits = [
            h for h in hits
            if expected_locus.lower() not in h["definition"].lower()
            and expected_locus.lower() not in h["accession"].lower()
        ]

    return {
        "rid": rid,
        "total_hits": len(hits),
        "off_target_hit_count": len(off_target_hits),
        "off_target_hits": off_target_hits[:10],
    }
