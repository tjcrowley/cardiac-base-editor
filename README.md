# Open mRNA Design Pipeline — Cardiac & Cancer

Free, open-source inference pipeline for personalized mRNA therapeutics. No pharma paywall, no API key, no institutional affiliation required. Runs locally at zero cost.

Built with [Biopunk Labs](https://biopunklab.com). Biology advisor: Elliot Roth.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---

## Arms

| Arm | Target | Input | Output |
|---|---|---|---|
| **Cardiac** | PCSK9 base editing | DNA sequence | Ranked guide RNAs |
| **Cancer** | Neoantigen vaccines | Somatic variants | Ranked MHC binders |

**Cardiac:** Designs and ranks guide RNAs for PCSK9 base editing — the same gene target as VERVE-102, the first in-human base editing clinical trial. Trained on 38,000 measured guide RNAs (Arbab et al. 2020, Nature Biotech).

**Cancer:** Predicts MHC binding affinity across patient HLA alleles to rank neoantigen candidates for personalized cancer vaccine design, using IEDB and ClinVar public datasets.

All outputs are deterministic and reproducible.

---

## Install

Packaged as an installable Python package, meant to run on a dedicated local
inference box (not a hosted service):

```bash
pip install -e .
mhcflurry-downloads fetch   # ~500MB MHC binding weights
```

## Run

```bash
python -m cardiac_base_editor.pipeline PCSK9   # reference-sequence guide ranking, as before

cbe consent grant elliot --scope PCSK9 --days 365   # genomic_intake: consent-gated real-genome intake
cbe run elliot PCSK9 --vcf /path/to/elliot.vcf
cbe query elliot "does this subject have any PCSK9 variants?" --vcf /path/to/elliot.vcf
cbe query elliot "check off-target risk for the top guide in PCSK9" --vcf /path/to/elliot.vcf
cbe query elliot "rank neoantigens for the variant at position 55039847 in PCSK9 for HLA alleles A0201 and B0702" --vcf /path/to/elliot.vcf

cbe lnp-predict formulations.json   # LNP delivery-efficacy prediction (LiON) — see setup below, not subject-gated
cbe plugins   # list optional model integrations and whether each is configured

cbe-web   # local web UI at http://127.0.0.1:8000
```

`cbe lnp-predict` needs its own environment (`chemprop==1.7.0` requires Python
3.7/3.8, incompatible with this package's own Python 3.10+ requirement) — see
`cardiac_base_editor/lnp_delivery.py`'s module docstring for the one-time
`pyenv`/venv setup and the `CBE_LNPDB_DIR` / `CBE_LION_VENV_DIR` env vars.

Real BE-DICT guide scoring (`cbe run`/`rank_guides`) runs in-process — clone
`github.com/uzh-dqbm-cmi/BEDICT-V2`, `pip install matplotlib seaborn` (plus
this package's own `pandas<3.0` pin, required for BEDICT-V2's code), and set
`CBE_BEDICTV2_DIR`. See `cardiac_base_editor/models/be_dict.py`'s docstring
for the real-model-vs-heuristic tradeoffs.

See `cardiac_base_editor/genomic_intake/README.md` for the full consent/audit/
retention model behind real-genome runs and queries — every genome is treated
as private: consent-gated, audit-logged, and fully purgeable on revocation.
The query engine (`cbe query` / the web UI's Query tab) is a local-LLM
(Ollama) front end over the same consent-gated functions — the model only
ever sees the question and already-computed structured results, never the
raw genome file.

---

## Adding a new model integration

Two small registries exist so a new integration only needs to touch its own
files — no other module needs manual updates to pick it up.

**If it needs an external tool/checkpoint the user has to set up** (most
model integrations — see `mrna_design.py`, `lnp_delivery.py`,
`cancer/tcr_binding.py`, `cancer/hla_typing.py` for real examples):

```python
from cardiac_base_editor.plugins import ExternalTool, register_tool, require_configured

TOOL = register_tool(ExternalTool(
    name="MyModel",
    env_vars=["CBE_MYMODEL_DIR"],
    setup_instructions="git clone ...\n  export CBE_MYMODEL_DIR=...",
    check=lambda: bool(os.environ.get("CBE_MYMODEL_DIR")),
))

def my_function(...):
    require_configured(TOOL)  # raises plugins.ToolNotConfigured with the setup instructions above
    ...
```
That's it — `cbe plugins` automatically lists it (configured or not, with
setup instructions), and any CLI command calling it should catch
`plugins.ToolNotConfigured` the way `cli.py`'s `cmd_lnp_predict` does.

If the model always has a working fallback instead of a hard requirement
(like `models/be_dict.py`'s heuristic), don't force it into this pattern —
warn and degrade instead, it's a genuinely different situation.

**If it should be queryable per-subject** (i.e. it belongs in
`query/engine.py`, consent-gated like `list_variants`/`rank_guides`/etc.):

```python
from cardiac_base_editor.query.registry import query_function

@query_function(hint="one-line description for the LLM router")
def my_query(subject_id: str, gene: str, vcf_path: str, operator: str = "query-engine") -> dict:
    """Docstring for Python callers — hint above is what the NL router sees instead."""
    subjects.require_consent(subject_id, gene, operator=operator)
    ...
```
`cbe query`'s LLM router and the web UI's generic
`/subjects/{id}/query/{function_name}` route both pick it up automatically —
no changes needed to `query/nl.py` or `web/app.py`. Not every capability
belongs here, though: `lnp_delivery.py`'s `predict_delivery_efficacy`
operates on a candidate formulation, not subject genomic data, so it's a
standalone CLI verb (`cbe lnp-predict`) instead of a query function — don't
force something through the consent gate that isn't actually about a
subject's genome.

---

## Roadmap

1. **Guide scoring:** Replace heuristics with a real BE-DICT model — **shipped** (`cardiac_base_editor/models/be_dict.py`), backed by [BEDICT-V2](https://github.com/uzh-dqbm-cmi/BEDICT-V2) (real pretrained checkpoints, MIT licensed), pluggable via `CBE_BEDICTV2_DIR`. Honest caveat: this wraps BEDICT-V2's "proportion" model, which predicts the *relative* likelihood among bystander-editing outcomes, not absolute efficiency — genuinely informative for bystander-risk guides (multiple targetable A's), but trivially returns 1.0 for single-target guides by construction. BEDICT-V2's separate absolute-efficiency CNN model would fix that but needs extra engineered features (RNA folding/MFE, melting temp) this pass didn't reproduce; falls back to the original heuristic for editors without a trained model (`ABE7.10`) or when unconfigured
2. **Off-target:** Genome-wide BLAST per guide — **shipped** (`cardiac_base_editor/models/off_target.py`), real queries against NCBI's public BLAST API; opt-in per-guide via `cbe query`'s `verify_off_target` (not run automatically — a real genome-wide BLAST search takes real wall-clock time)
3. **Protein consequence:** ESM-2 call on resulting amino acid change — **shipped** (`cardiac_base_editor/models/protein_consequence.py`), used by `cbe query`'s `explain_variant`
4. **Cancer — MHC binding:** **shipped** (`cardiac_base_editor/cancer/`) — somatic variant → candidate neoantigen peptides → MHC-I binding ranking via mhcflurry, exposed as `cbe query`'s `rank_neoantigens`. Still requires HLA alleles to be supplied directly (see #4b)
4b. **Cancer — HLA typing:** **interface shipped** (`cardiac_base_editor/cancer/hla_typing.py`), a real subprocess wrapper around [arcasHLA](https://github.com/RabadanLab/arcasHLA) — not live-verified here, since it genotypes real patient RNA-seq (BAM) and no synthetic sample would produce a meaningful call; box-only until run against real tumor RNA-seq
5. **Cancer — T-cell response:** **interface shipped** (`cardiac_base_editor/cancer/tcr_binding.py`), a real subprocess wrapper around [pMTnet](https://github.com/tianshilu/pMTnet) — not live-verified here, since pMTnet is pinned to TensorFlow 1.x/Keras 2.2.4/numpy 1.16.3, none of which have wheels for this environment's Python 3.11/arm64; needs a compatible runtime (likely Docker w/ an old TF1 image) to actually run
6. **mRNA sequence design:** **shipped** (`cardiac_base_editor/mrna_design.py`) — real [LinearDesign](https://github.com/LinearDesignSoftware/LinearDesign) integration, compiles clean with `make` (Apple clang/GCC, no exotic deps) and runs by calling the compiled binary directly (its bundled python2 CLI wrapper is skipped entirely). Exposed as `cbe query`'s `design_mrna_payload` — verified end-to-end, including translating the returned mRNA back and confirming it encodes the same protein
7. **Delivery:** **shipped** (`cardiac_base_editor/lnp_delivery.py`) — real integration of [LiON](https://www.nature.com/articles/s41587-024-02490-y) (Witten et al. 2024), a pretrained `chemprop` model trained on the public [LNPDB](https://github.com/evancollins1/LNPDB) dataset. Predictions verified against LNPDB's own bundled test results to floating-point precision. Exposed as `cbe lnp-predict`, not part of `cbe query` — it predicts a candidate LNP formulation's delivery efficacy, not anything about subject genomic data, so it isn't consent-gated

---

## Data Sources

| Dataset | What it provides | Access |
|---|---|---|
| Arbab et al. 2020 (Nature Biotech) | 38k guide RNAs + base editing outcomes | Public |
| IEDB | MHC binding measurements | Public (iedb.org) |
| TCGA | Tumor somatic variants + clinical outcomes | dbGaP (controlled) |
| ClinVar | Pathogenic variant catalog | Public |
| VERVE-102 trial data | Clinical base editing outcomes | Locked (pharma) |

---

## Support This Project

We're competing in the [Artizen Season 6 Crescendo Fund Drive](https://artizen.fund/index/p/open-mrna-design-pipeline--cardiac--cancer?season=6). Artifact purchases fund the next milestone and keep this pipeline open forever.

---

## License

This project is licensed under the **GNU General Public License v3.0**.

You are free to use, modify, and distribute this software. Any derivative work must also be released under GPL v3 — improvements stay open, always.

See [LICENSE](./LICENSE) for full terms.

