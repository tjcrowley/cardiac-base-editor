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

```bash
pip install transformers torch mhcflurry biopython
mhcflurry-downloads fetch   # ~500MB MHC binding weights
```

## Run

```bash
python pipeline.py
```

---

## Roadmap

1. **Guide scoring:** Replace heuristics with [BE-DICT](https://github.com/hui-liang/BE-DICT) transformer (next funded milestone)
2. **Off-target:** Genome-wide BLAST or crisprSFM API call per guide
3. **Protein consequence:** ESM-2 or AlphaMissense call on resulting amino acid change
4. **Cancer — HLA typing:** Patient alleles from `optitype` or `arcasHLA` on tumor WES data
5. **Cancer — T-cell response:** pMTnet call after MHC binding filter
6. **mRNA sequence design:** LinearDesign for codon-optimized mRNA sequence of final payload
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

