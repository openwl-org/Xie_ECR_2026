---
title: Xie et al. (2026) — Full Pipeline
description: Complete mtDNA-CN depletion transcriptomics analysis — normalization, 6 regressions, pattern classification, GO enrichment, figures, and HTML report
---

# Setup — Add Specialist

Add the mtDNA Transcriptomics Specialist to your notebook's workflow team. This agent implements the complete Xie et al. (2026) analysis pipeline from the Arking Lab at Johns Hopkins.

:::{prompt}
Add the mtdna-transcriptomics-specialist to the team.
:::

---

# Full Pipeline

Run the complete Xie et al. (2026) analysis in a single execution:

*"Transcriptomic analysis of cells following decreased mitochondrial DNA-copy number reveals distinct transcriptional response patterns"*

:::{note}
**Pipeline steps (11 tools, ~15 min):**
1. `prepare_expression_data` — TMM normalize, PCA QC, compute Med_log2FC_MT
2. `regress_dose/cn/cn_spline/cn_quadratic/mt` — 5 regression models (parallel)
3. `regress_longitudinal` — per-timepoint Treatment vs Control
4. `classify_response_patterns` — Linear/Switch/Delayed/None via likelihood templates
5. `run_go_enrichment` — stratified enrichGO + 15-way goana + rrvgo treemaps
6. `generate_publication_figures` — 9 figures (PNG)
7. `generate_html_report` — self-contained HTML with embedded figures

**Input files** (already in vault):
- `merged.csv` — RNA-seq count matrix (genes x samples)
- `AllSamples.Randomized.RNAseq.csv` — sample metadata bridge file
- `ChemicalProjectMatchUpDeltaCT.xlsx` — qPCR deltaCT measurements
- `Human.MitoCarta3.0.xlsx` — MitoCarta 3.0 database
- `mart_export.txt` — Ensembl gene annotations (BioMart export)
:::

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

After the pipeline completes (~15 min), ask for a summary of the findings. The specialist will report actual values from the output RDS and JSON files.

:::{tip}
The paper found ~1,500 significant genes with distinct Linear, Switch, and Delayed response patterns. Compensatory upregulation of central dogma genes and a glycolytic shift were key findings.
:::

:::{prompt}
Summarize the complete results: how many genes tested, significant per model, pattern classification breakdown, top enriched GO terms, and which figures were generated.
:::
