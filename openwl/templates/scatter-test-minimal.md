---
title: Xie et al. (2026) — Quick Scatter Demo (2 Regressions + Figures)
description: Fast demo of parallel gene-level scatter. Runs data prep, 2 regressions (dose + CN) each split across 2 parallel pods, merge chunks, pattern classification, and publication figures. ~5 min, ~9 pods.
---

# Setup

Add the mtDNA Transcriptomics Specialist to your notebook's workflow team.

:::{prompt}
Add the mtdna-transcriptomics-specialist to the team.
:::

---

# Quick Scatter Pipeline

Run a compact version of the Xie et al. pipeline — 2 regressions with gene-level scatter, pattern classification, and publication figures.

:::{note}
**Pipeline overview (6 tools, ~9 pods, ~5 min):**

**Phase 1 — Data Preparation**
- `prepare_expression_data` — TMM-normalize RNA-seq counts, parse sample metadata, PCA outlier removal

**Phase 2 — Gene-wise Regressions (2 tools x 2 gene chunks = 4 parallel pods)**
- `regress_dose` — log(CPM) ~ Dose (scaled 0-6); OLS + mixed model
- `regress_cn` — log(CPM) ~ deltaCT_Avg; linear + quadratic, OLS + mixed (4 result sets)
- Each regression runs with `scatter: "chunks:2"` — ~11,800 genes split across 2 parallel pods

**Phase 2b — Merge Chunks**
- `merge_chunks` — merge chunked regression outputs (depends on ALL scattered regressions)

**Phase 3 — Pattern Classification**
- `classify_response_patterns` — Classify genes into Linear/Switch/Delayed/None response patterns

**Phase 4 — Figures**
- `generate_publication_figures` — 9 publication-quality PNG figures

**Input files** (all in vault):
- `merged.csv` — RNA-seq count matrix
- `AllSamples.Randomized.RNAseq.csv` — sample metadata bridge file
- `ChemicalProjectMatchUpDeltaCT.xlsx` — qPCR deltaCT measurements
- `Human.MitoCarta3.0.xlsx` — MitoCarta 3.0 database
- `mart_export.txt` — Ensembl gene annotations
:::

:::{prompt}
Run a quick scatter pipeline with these steps:
1. prepare_expression_data
2. regress_dose with scatter: "chunks:2"
3. regress_cn with scatter: "chunks:2"
4. merge_chunks (depends on steps 2-3)
5. classify_response_patterns (depends on step 4)
6. generate_publication_figures (depends on step 5)

Use these input files:
- countFile: /pfs/promptable-public-test-data/arking-mtdna/merged.csv
- bridgeFile: /pfs/promptable-public-test-data/arking-mtdna/AllSamples.Randomized.RNAseq.csv
- deltaCTFile: /pfs/promptable-public-test-data/arking-mtdna/ChemicalProjectMatchUpDeltaCT.xlsx
- mitoCartaFile: /pfs/promptable-public-test-data/arking-mtdna/Human.MitoCarta3.0.xlsx
- annotationFile: /pfs/promptable-public-test-data/arking-mtdna/mart_export.txt
:::
