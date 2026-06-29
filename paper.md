---
title: 'cardiac-base-editor: An open pipeline for designing base editor interventions in inherited cardiac conditions'
tags:
  - Python
  - bioinformatics
  - base editing
  - CRISPR
  - cardiac genetics
  - genomics
authors:
  - name: Darren Mckeeman
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Elliot Roth
    affiliation: 2
affiliations:
  - name: Independent Researcher, San Francisco, CA, USA
    index: 1
  - name: Biopunk Labs, San Francisco, CA, USA
    index: 2
date: 29 June 2026
bibliography: paper.bib
---

# Summary

`cardiac-base-editor` is an open source pipeline for designing base editor interventions targeting inherited cardiac arrhythmias and cardiomyopathies. Given a pathogenic variant from ClinVar, the pipeline performs guide RNA design, predicts base editing efficiency using BE-DICT [@arbab2023], scores off-target risk, and produces a structured output package (JSON and PDF) suitable for direct use in wet lab validation workflows. The pipeline is implemented in Python, licensed under GPL v3, and designed to run without institutional software dependencies.

Base editing — the programmable conversion of one DNA nucleotide to another without introducing a double-strand break — has emerged as a promising therapeutic strategy for the large class of hereditary cardiac conditions caused by point mutations [@komor2016; @gaudelli2017]. Conditions including Long QT syndrome (SCN5A, KCNQ1), hypertrophic cardiomyopathy (MYH7), and Brugada syndrome involve well-characterized pathogenic variants that are, in principle, correctable by adenine or cytosine base editors delivered via mRNA. The computational design of such interventions — selecting guide RNAs, predicting editing outcomes, and minimizing off-target activity — currently requires chaining together multiple tools that do not interoperate and are often closed source or require institutional licenses.

# Statement of Need

Researchers designing base editor experiments for cardiac applications currently rely on a fragmented set of tools: BE-Hive [@arbab2020] or CRISPRscan [@moreno2015] for guide design, BE-DICT [@arbab2023] for efficiency prediction, and Cas-OFFinder [@bae2014] for off-target analysis. These tools are not integrated, require significant scripting expertise to connect, and produce outputs in incompatible formats. Each laboratory re-implements the same pipeline independently, producing workflows that are difficult to reproduce, compare, or build upon.

No open, unified pipeline exists that takes a clinical variant as input and produces a validated, reproducible design package as output. Commercial platforms that offer integrated base editor design tools are closed source, expensive, and inaccessible to independent researchers and laboratories in lower-resource settings. The absence of such a tool represents a reproducibility gap in a rapidly moving therapeutic area.

`cardiac-base-editor` addresses this gap by providing a single entry point — a ClinVar variant ID or VCF record — and automating each step of the design process through to a structured output that can be passed directly to a wet lab. By releasing under GPL v3, the pipeline ensures that improvements made by any user or institution remain available to the broader community.

# Usage

The pipeline accepts a ClinVar variant identifier and produces a design package:

```bash
python pipeline.py --variant NM_000335.5:c.1129T>C --output ./output/
```

Output includes:

- `design.json` — guide RNA candidates with predicted efficiency scores and off-target risk ratings
- `protocol.pdf` — a formatted wet lab handoff document mapping computational predictions to standard assay steps (plasmid synthesis, cell line selection, editing efficiency assay)
- `summary.tsv` — tabular results suitable for downstream analysis

Wet lab validation of pipeline outputs is ongoing at Biopunk Labs (San Francisco, CA), with initial results expected against a panel of SCN5A and KCNQ1 variants.

# Acknowledgements

The authors thank the developers of BE-DICT, BE-Hive, CRISPRscan, and Cas-OFFinder for maintaining open research tools that made this work possible. Wet lab validation is supported by Biopunk Labs.

# References
