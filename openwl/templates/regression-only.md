---
title: Xie et al. (2026) — Regression Analysis
description: Data preparation, 6 regression models, and response pattern classification (no enrichment or figures)
---

# Setup

:::{prompt}
Add the mtdna-transcriptomics-specialist to the team.
:::

---

# Step 1 — Prepare Data

:::{prompt}
Prepare the expression data using:
- countFile: /pfs/promptable-public-test-data/arking-mtdna/merged.csv
- bridgeFile: /pfs/promptable-public-test-data/arking-mtdna/AllSamples.Randomized.RNAseq.csv
- deltaCTFile: /pfs/promptable-public-test-data/arking-mtdna/ChemicalProjectMatchUpDeltaCT.xlsx
- annotationFile: /pfs/promptable-public-test-data/arking-mtdna/mart_export.txt
:::

---

# Step 2 — Regressions (Parallel)

:::{prompt}
Run all 6 regressions in parallel: regress_dose, regress_cn, regress_cn_spline, regress_cn_quadratic, regress_mt, and regress_longitudinal.
:::

---

# Step 3 — Classify Patterns

:::{prompt}
Classify response patterns using /pfs/promptable-public-test-data/arking-mtdna/Human.MitoCarta3.0.xlsx and /pfs/promptable-public-test-data/arking-mtdna/mart_export.txt with the regression results from the previous step.
:::

---

# Results

:::{prompt}
Summarize: how many genes significant per model (Holm p<0.05)? Pattern breakdown (Linear/Switch/Delayed/None)? Which central dogma genes had the strongest CN association?
:::
