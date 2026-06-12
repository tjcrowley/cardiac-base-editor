# mRNA Therapeutic Pipeline — Prototype

Two-arm inference pipeline for personalized mRNA therapeutics.

## Arms

| Arm | Target | Input | Output |
|---|---|---|---|
| **Cardiac** | PCSK9 base editing | DNA sequence | Ranked guide RNAs |
| **Cancer** | Neoantigen vaccines | Somatic variants | Ranked MHC binders |

## Install

```bash
pip install transformers torch mhcflurry biopython
mhcflurry-downloads fetch   # ~500MB MHC binding weights
```

## Run

```bash
python pipeline.py
```

## Production Upgrades (in priority order)

1. **Guide scoring:** Replace heuristics with [BE-DICT](https://github.com/hui-liang/BE-DICT) model call
2. **Off-target:** Add genome-wide BLAST or crisprSFM API call per guide
3. **Protein consequence:** Add ESM-2 or AlphaMissense call on the resulting amino acid change
4. **Cancer — HLA typing:** Patient alleles from `optitype` or `arcasHLA` on tumor WES data
5. **Cancer — T-cell response:** Add pMTnet call after MHC binding filter
6. **mRNA sequence design:** Add LinearDesign for codon-optimized mRNA sequence of final payload
7. **Delivery:** LNP formulation predictor (literature: Moderna/Inivio published datasets)

## Data Sources for Fine-tuning

| Dataset | What it provides | Access |
|---|---|---|
| Arbab et al. 2020 (Nature Biotech) | 38k guide RNAs + base editing outcomes | Public |
| IEDB | MHC binding measurements | Public (iedb.org) |
| TCGA | Tumor somatic variants + clinical outcomes | dbGaP (controlled) |
| ClinVar | Pathogenic variant catalog | Public |
| VERVE-102 trial data | Clinical base editing outcomes | Locked (pharma) |
