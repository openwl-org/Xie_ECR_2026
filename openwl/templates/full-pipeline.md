---
title: Xie et al. (2026) — Full Pipeline
description: Complete mtDNA-CN depletion transcriptomics analysis — normalization, 6 regressions, pattern classification, GO enrichment, figures, and HTML report
---

# Setup — Add Specialist

:::{prompt}
Add the mtdna-transcriptomics-specialist to the team.
:::

---

# Full Pipeline

Run the complete Xie et al. (2026) analysis from the Arking Lab:

*"Transcriptomic analysis of cells following decreased mitochondrial DNA-copy number reveals distinct transcriptional response patterns"*

:::{note}
**Pipeline steps:**
1. `prepare_expression_data` — TMM normalize, PCA QC, compute Med_log2FC_MT
2. `regress_dose/cn/cn_spline/cn_quadratic/mt` — 5 regression models (parallel)
3. `regress_longitudinal` — per-timepoint Treatment vs Control
4. `classify_response_patterns` — Linear/Switch/Delayed/None via likelihood templates
5. `run_go_enrichment` — stratified enrichGO + 15-way goana + rrvgo treemaps
6. `generate_publication_figures` — 9 figures (PNG)
7. `generate_html_report` — self-contained HTML with embedded figures
:::

**Input files** (in vault):
- `merged.csv` — RNA-seq count matrix
- `AllSamples.Randomized.RNAseq.csv` — sample metadata
- `ChemicalProjectMatchUpDeltaCT.xlsx` — qPCR deltaCT
- `Human.MitoCarta3.0.xlsx` — MitoCarta 3.0
- `mart_export.txt` — Ensembl gene annotations

:::{prompt}
Run the full Xie et al. mtDNA depletion pipeline using:
- countFile: /pfs/promptable-public-test-data/arking-mtdna/merged.csv
- bridgeFile: /pfs/promptable-public-test-data/arking-mtdna/AllSamples.Randomized.RNAseq.csv
- deltaCTFile: /pfs/promptable-public-test-data/arking-mtdna/ChemicalProjectMatchUpDeltaCT.xlsx
- mitoCartaFile: /pfs/promptable-public-test-data/arking-mtdna/Human.MitoCarta3.0.xlsx
- annotationFile: /pfs/promptable-public-test-data/arking-mtdna/mart_export.txt

Run all steps: prepare_expression_data, all regressions (dose, CN, CN spline, CN quadratic, MT, longitudinal), classify_response_patterns, run_go_enrichment, generate_publication_figures, and generate_html_report.
:::

---

# Results

:::{prompt}
Summarize the complete results: how many genes tested, significant per model, pattern classification breakdown, top enriched GO terms, and which figures were generated.
:::
