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

cbe-web   # local web UI at http://127.0.0.1:8000
```

See `cardiac_base_editor/genomic_intake/README.md` for the full consent/audit/
retention model behind real-genome runs and queries — every genome is treated
as private: consent-gated, audit-logged, and fully purgeable on revocation.
The query engine (`cbe query` / the web UI's Query tab) is a local-LLM
(Ollama) front end over the same consent-gated functions — the model only
ever sees the question and already-computed structured results, never the
raw genome file.

---

## Roadmap

1. **Guide scoring:** Replace heuristics with [BE-DICT](https://github.com/hui-liang/BE-DICT) transformer — **interface shipped** (`cardiac_base_editor/models/be_dict.py`), pluggable via `CBE_BEDICT_CHECKPOINT`; falls back to the existing heuristic until a real checkpoint is available (next funded milestone)
2. **Off-target:** Genome-wide BLAST per guide — **shipped** (`cardiac_base_editor/models/off_target.py`), real queries against NCBI's public BLAST API; opt-in per-guide via `cbe query`'s `verify_off_target` (not run automatically — a real genome-wide BLAST search takes real wall-clock time)
3. **Protein consequence:** ESM-2 call on resulting amino acid change — **shipped** (`cardiac_base_editor/models/protein_consequence.py`), used by `cbe query`'s `explain_variant`
4. **Cancer — MHC binding:** **shipped** (`cardiac_base_editor/cancer/`) — somatic variant → candidate neoantigen peptides → MHC-I binding ranking via mhcflurry, exposed as `cbe query`'s `rank_neoantigens`. Still requires HLA alleles to be supplied directly (see #4b)
4b. **Cancer — HLA typing:** Deriving those alleles from a patient's own tumor WES via `optitype` or `arcasHLA`, rather than supplying them directly — not yet built
5. **Cancer — T-cell response:** pMTnet call after MHC binding filter
6. **mRNA sequence design:** [LinearDesign](https://github.com/LinearDesignSoftware/LinearDesign) for codon-optimized mRNA sequence of final payload — real tool exists but is C++, needs compiling from source; bigger lift, not yet built
7. **Delivery:** LNP formulation predictor (Moderna/Inivio published datasets)

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

