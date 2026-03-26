---
title: Xie et al. (2026) — Regression Analysis
description: Data preparation, 6 regression models, and response pattern classification (no enrichment or figures)
---

# Setup

Add the mtDNA Transcriptomics Specialist to your team. This template runs only the statistical analysis portion of the pipeline — data prep, regressions, and pattern classification.

:::{prompt}
Add the mtdna-transcriptomics-specialist to the team.
:::

---

# Step 1 — Prepare Data

Read raw RNA-seq counts, parse sample metadata from SampleName, TMM-normalize with edgeR, remove PCA outliers (>4 SD on first 10 PCs), and compute the Med_log2FC_MT proxy for each sample.

:::{note}
**Outputs:** `Experiment_Samples.rds` (sample metadata with Dose, Chemical, Timepoint, deltaCT), `RNA_TMM.rds` (TMM-normalized log-CPM expression matrix), `mtDNA.rds` (qPCR measurements).
:::

:::{prompt}
Prepare the expression data using:
- countFile: /pfs/promptable-public-test-data/arking-mtdna/merged.csv
- bridgeFile: /pfs/promptable-public-test-data/arking-mtdna/AllSamples.Randomized.RNAseq.csv
- deltaCTFile: /pfs/promptable-public-test-data/arking-mtdna/ChemicalProjectMatchUpDeltaCT.xlsx
- annotationFile: /pfs/promptable-public-test-data/arking-mtdna/mart_export.txt
:::

---

# Step 2 — Regressions (Parallel)

Run all 6 gene-wise regression models in parallel. Each fits log(CPM) ~ predictor for every gene, testing the association between mitochondrial perturbation and nuclear gene expression.

:::{important}
**6 models, 5 parallel K8s jobs:**
1. `regress_dose` — GE ~ EtBr Dose (scaled 0-6)
2. `regress_cn` — GE ~ deltaCT_Avg (linear + quadratic, OLS + mixed)
3. `regress_cn_spline` — GE ~ ns(deltaCT_Avg, df=2) and df=3
4. `regress_cn_quadratic` — AIC/LRT model comparison (linear vs quad vs splines)
5. `regress_mt` — GE ~ Med_log2FC_MT (median MT gene fold-change)
6. `regress_longitudinal` — GE ~ Experiment (Control vs EtBr) per timepoint
:::

:::{prompt}
Run all 6 regressions in parallel: regress_dose, regress_cn, regress_cn_spline, regress_cn_quadratic, regress_mt, and regress_longitudinal.
:::

---

# Step 3 — Classify Patterns

Combine all regression results and classify each gene into one of four response trajectory patterns using likelihood templates. The classification uses scaled fold-change at 4 dose-binned timepoints with dnorm product likelihoods.

:::{note}
**Pattern templates:**
- **Linear** `(0, 0.5, 1, 0.5)` — proportional to CN depletion, partial recovery
- **Switch** `(0, 1, 1, 0)` — mirrors MT gene expression, threshold-driven, rapid recovery
- **Delayed** `(0, 0.25, 1, 1)` — lagged onset, persists during recovery (stress/signaling genes)
- **None** — ambiguous (top two likelihoods lack >2x separation)
:::

:::{prompt}
Classify response patterns using /pfs/promptable-public-test-data/arking-mtdna/Human.MitoCarta3.0.xlsx and /pfs/promptable-public-test-data/arking-mtdna/mart_export.txt with the regression results from the previous step.
:::

---

# Results

Review the statistical results. The paper reported ~1,500 significant genes (Holm p<0.05) with distinct response pattern distributions across MitoCarta-annotated mitochondrial genes.

:::{tip}
Ask about specific gene families — central dogma genes (POLG, TWINKLE, TFAM), OXPHOS complexes (I-V), or glycolytic enzymes (SLC2A1, HK1, GPI) — to explore the biology.
:::

:::{prompt}
Summarize: how many genes significant per model (Holm p<0.05)? Pattern breakdown (Linear/Switch/Delayed/None)? Which central dogma genes had the strongest CN association?
:::
