"""
Local web UI for genomic_intake. Single-operator, localhost-only by default —
this is a friendlier front end over the same consent/audit/extract/pipeline
logic used by genomic_intake.cli, not a hosted multi-tenant service.

Run:
    uvicorn genomic_intake.web.app:app --reload
    (binds to 127.0.0.1:8000 by default via uvicorn's own default host)
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

from fastapi import FastAPI, File, Form, Request, UploadFile
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from genomic_intake import audit, subjects
from genomic_intake.extract import build_personalized_cds
from pipeline import KNOWN_TARGETS, run as pipeline_run

BASE_DIR = Path(__file__).resolve().parent

app = FastAPI(title="genomic_intake")
app.mount("/static", StaticFiles(directory=BASE_DIR / "static"), name="static")
templates = Jinja2Templates(directory=BASE_DIR / "templates")

WEB_OPERATOR = "web-ui"


def _retention_status(record: subjects.ConsentRecord) -> str:
    if not record.active:
        return "revoked"
    if subjects.retention_expired(record):
        return "retention expired"
    return "active"


def _render(request: Request, template: str, context: dict) -> HTMLResponse:
    context["request"] = request
    return templates.TemplateResponse(template, context)


@app.get("/", response_class=HTMLResponse)
def dashboard(request: Request):
    all_subjects = subjects.list_all()
    rows = [
        {"subject_id": sid, "record": rec, "status": _retention_status(rec)}
        for sid, rec in sorted(all_subjects.items())
    ]
    return _render(request, "dashboard.html", {"rows": rows})


@app.post("/subjects", response_class=HTMLResponse)
def create_subject(
    request: Request,
    subject_id: str = Form(...),
    scope: str = Form("*"),
    retention_days: int = Form(365),
):
    scope_list = [s.strip() for s in scope.split(",") if s.strip()] or ["*"]
    subjects.grant_consent(subject_id, scope=scope_list, retention_days=retention_days, operator=WEB_OPERATOR)
    return dashboard(request)


@app.get("/subjects/{subject_id}", response_class=HTMLResponse)
def subject_detail(request: Request, subject_id: str):
    record = subjects.get_consent(subject_id)
    history = list(reversed(audit.read_for_subject(subject_id)))
    return _render(
        request,
        "subject_detail.html",
        {
            "subject_id": subject_id,
            "record": record,
            "status": _retention_status(record) if record else "no consent record",
            "history": history,
            "known_genes": sorted(KNOWN_TARGETS.keys()),
        },
    )


@app.post("/subjects/{subject_id}/consent", response_class=HTMLResponse)
def grant_consent(request: Request, subject_id: str, scope: str = Form("*"), retention_days: int = Form(365)):
    scope_list = [s.strip() for s in scope.split(",") if s.strip()] or ["*"]
    subjects.grant_consent(subject_id, scope=scope_list, retention_days=retention_days, operator=WEB_OPERATOR)
    return subject_detail(request, subject_id)


@app.post("/subjects/{subject_id}/revoke", response_class=HTMLResponse)
def revoke_consent(request: Request, subject_id: str, purge_data: bool = Form(True)):
    subjects.revoke_consent(subject_id, operator=WEB_OPERATOR, purge_data=purge_data)
    return subject_detail(request, subject_id)


@app.post("/subjects/{subject_id}/run", response_class=HTMLResponse)
async def run_pipeline(
    request: Request,
    subject_id: str,
    gene: str = Form(...),
    vcf: UploadFile = File(...),
    editor: str = Form("ABE8e"),
):
    try:
        subjects.require_consent(subject_id, gene, operator=WEB_OPERATOR)
    except subjects.ConsentError as e:
        return _render(request, "run_result.html", {"error": str(e), "subject_id": subject_id, "gene": gene})

    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as tmp:
        shutil.copyfileobj(vcf.file, tmp)
        tmp_path = tmp.name

    try:
        transcript_id = KNOWN_TARGETS.get(gene.upper(), gene)
        personalized_cds = build_personalized_cds(transcript_id, tmp_path)
        results = pipeline_run(sequence=personalized_cds, editor_name=editor, top_n=10)
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    return _render(request, "run_result.html", {"subject_id": subject_id, "gene": gene, "results": results, "error": None})
