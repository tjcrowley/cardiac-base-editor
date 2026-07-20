"""
Local web UI for cardiac_base_editor. Single-operator, localhost-only by
default — a friendlier front end over the same consent/audit/extract/pipeline/
query-engine logic used by cardiac_base_editor.cli, not a hosted
multi-tenant service.

Run:
    cbe-web
    (or: uvicorn cardiac_base_editor.web.app:app --reload)
    binds to 127.0.0.1:8000 by default
"""

from __future__ import annotations

import inspect
import json
import shutil
import tempfile
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from cardiac_base_editor.genomic_intake import audit, subjects
from cardiac_base_editor.genomic_intake.extract import build_personalized_cds
from cardiac_base_editor.pipeline import KNOWN_TARGETS, run as pipeline_run
from cardiac_base_editor.query import nl as query_nl
from cardiac_base_editor.query.registry import QUERY_FUNCTIONS

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
    return templates.TemplateResponse(request, template, context)


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


@app.post("/subjects/{subject_id}/query", response_class=HTMLResponse)
async def query_subject(
    request: Request,
    subject_id: str,
    question: str = Form(...),
    vcf: UploadFile = File(...),
):
    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as tmp:
        shutil.copyfileobj(vcf.file, tmp)
        tmp_path = tmp.name

    try:
        answer_text = query_nl.answer(subject_id, question, tmp_path)
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    return _render(request, "query_result.html", {"question": question, "answer": answer_text})


@app.post("/subjects/{subject_id}/query/{function_name}", response_class=HTMLResponse)
async def query_function_structured(
    request: Request,
    subject_id: str,
    function_name: str,
    vcf: UploadFile = File(...),
    args: str = Form("{}"),
):
    """
    Structured entry point for any @query_function-registered function (see
    query/registry.py) — the generic counterpart to /query's NL router. Every
    current and future query function gets this for free, no per-function
    route needed.
    """
    if function_name not in QUERY_FUNCTIONS:
        raise HTTPException(status_code=404, detail=f"Unknown query function: {function_name}")

    func = QUERY_FUNCTIONS[function_name].func
    try:
        parsed_args = json.loads(args)
    except json.JSONDecodeError:
        return _render(request, "query_function_result.html", {
            "function_name": function_name, "error": f"args is not valid JSON: {args!r}", "result": None,
        })

    valid_params = set(inspect.signature(func).parameters) - {"subject_id", "vcf_path", "operator"}
    parsed_args = {k: v for k, v in parsed_args.items() if k in valid_params}

    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as tmp:
        shutil.copyfileobj(vcf.file, tmp)
        tmp_path = tmp.name

    try:
        result = func(subject_id=subject_id, vcf_path=tmp_path, operator=WEB_OPERATOR, **parsed_args)
        error = None
    except Exception as e:
        result = None
        error = str(e)
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    return _render(request, "query_function_result.html", {
        "function_name": function_name, "result": result, "error": error,
    })


def serve() -> None:
    """Console-script entry point (`cbe-web`)."""
    import uvicorn

    uvicorn.run("cardiac_base_editor.web.app:app", host="127.0.0.1", port=8000)
