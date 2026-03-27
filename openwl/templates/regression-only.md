---
title: Xie et al. (2026) — Regression Analysis
description: Data preparation, 6 regression models, and response pattern classification — the statistical core of the mtDNA-CN depletion transcriptomics pipeline (no enrichment or figures)
---

# Setup — Add Specialist

Add the mtDNA Transcriptomics Specialist to your notebook's workflow team. This agent implements the statistical analysis portion of the Xie et al. (2026) pipeline from the Arking Lab at Johns Hopkins — "Transcriptomic analysis of cells following decreased mitochondrial DNA-copy number reveals distinct transcriptional response patterns."

This template runs the statistical core only: data preparation, all 6 regression models, and response pattern classification. Use the full pipeline template if you also want GO enrichment, publication figures, and the HTML report.

:::{prompt}
Add the mtdna-transcriptomics-specialist to the team.
:::

---

# Step 1 — Prepare Expression Data

Read raw RNA-seq counts, parse sample metadata, TMM-normalize, remove PCA outliers, and compute the mitochondrial RNA proxy (Med_log2FC_MT).

:::{note}
**What this step does in detail:**

1. **Read raw counts** — `merged.csv` contains genes x 112 samples. Column names encode experimental metadata in the `SampleName` format (Chemical_Dose_Curve_Time_Rep).
2. **Parse metadata** — Extracts Chemical (EtBr, RESV, DMSO, NC), Dose (0-150 ng/mL, scaled 0-6), Curve (Dosage, Treatment, Recovery), Timepoint (0-192h), and Replicate from SampleName.
3. **Bridge file merge** — `AllSamples.Randomized.RNAseq.csv` maps RNA-seq sample IDs to the qPCR sample IDs used in `ChemicalProjectMatchUpDeltaCT.xlsx`. This is how gene expression and mtDNA copy number measurements get linked.
4. **Filter to EtBr** — Keeps only EtBr-treated samples plus their matched controls across Dosage (6 doses x 2 reps x 4 lanes), Treatment (7 timepoints x 2 reps), and Recovery (4 timepoints x 2 reps) arms. Result: 126 samples.
5. **TMM normalization** — edgeR's `calcNormFactors(method="TMM")` on genes with median count > 49 in control samples. Only keeps ENSG-prefixed genes (removes spike-ins, unmapped). Result: ~11,800 genes.
6. **PCA outlier removal** — Runs PCA on first 10 principal components, removes any sample > 4 SD from the mean on any PC. The Xie et al. data is clean — expect 0 outliers removed.
7. **Compute Med_log2FC_MT** — For each sample, computes the median log2 fold-change across the 13 MT protein-coding genes (MT-ND1 through MT-CYB) relative to condition-matched controls. This is used as the continuous mitochondrial RNA proxy in `regress_mt`.
8. **Exclude bad samples** — Removes RESV_NC and ETBR_NC1 (quality control failures noted in the paper).

**Outputs:** 3 RDS files
- `Experiment_Samples.rds` — Sample metadata data.frame (126 rows x ~20 columns: Dose, Chemical, Timepoint, deltaCT_Avg, Med_log2FC_MT, etc.)
- `RNA_TMM.rds` — TMM-normalized log-CPM expression matrix (~11,800 genes x 126 samples)
- `mtDNA.rds` — qPCR deltaCT measurements linked to sample IDs
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

Run all 6 gene-wise regression models. The 5 dosage regressions each run with `scatter: "chunks:4"` — splitting ~11,800 genes across 4 parallel pods for ~4x speedup (20 pods total). The longitudinal regression runs as a single pod.

:::{important}
**6 models, 5 tools x 4 gene chunks + 1 longitudinal = 21 parallel pods:**

1. **`regress_dose`** (scatter: "chunks:4") — log(CPM) ~ Dose (scaled 0-6 from EtBr concentration 0-150 ng/mL)
   - OLS on lane-averaged replicates (simpler, fewer samples)
   - Mixed model with `(1|SampleName)` random effect (accounts for lane replication)
   - Tests: does gene expression change with EtBr dose?
   - Key output columns: `Estimate`, `Std. Error`, `p_holm` for both OLS and mixed

2. **`regress_cn`** (scatter: "chunks:4") — log(CPM) ~ deltaCT_Avg (qPCR-measured mtDNA copy number)
   - Linear: OLS + mixed model (2 result sets)
   - Quadratic: OLS + mixed model with deltaCT + deltaCT² terms (2 more result sets)
   - Tests: does gene expression associate with actual measured mtDNA-CN?
   - The CN model is the paper's primary significance test — Holm p < 0.05 on CN mixed model defines "significant genes"

3. **`regress_cn_spline`** (scatter: "chunks:4") — log(CPM) ~ ns(deltaCT_Avg, df=2) and ns(deltaCT_Avg, df=3)
   - Natural spline fits capture non-linear CN-expression relationships
   - df=2: one internal knot (flexible S-curve)
   - df=3: two internal knots (more complex curvature)
   - Tests: do some genes have threshold or saturation effects in their CN response?

4. **`regress_cn_quadratic`** (scatter: "chunks:4") — AIC model comparison across all CN models
   - Compares 4 models: linear, quadratic, ns(df=2), ns(df=3)
   - Likelihood ratio tests: P12 (linear vs quad), P13 (linear vs ns2), P14 (linear vs ns3)
   - AIC selection identifies best-fitting model per gene
   - Expected: ~500 genes where non-linear models significantly improve over linear

5. **`regress_mt`** (scatter: "chunks:4") — log(CPM) ~ Med_log2FC_MT (median MT gene fold-change proxy)
   - OLS + mixed model
   - Tests: does gene expression associate with the transcriptional state of the mitochondrial genome?
   - High concordance with CN model expected (both measure mitochondrial perturbation, different proxies)

6. **`regress_longitudinal`** — log(CPM) ~ Experiment (Control vs EtBr) at each of 7 timepoints
   - 7 separate contrasts: 0h, 24h, 48h, 72h, 144h, 168h, 192h
   - Both Treatment arm (continuous EtBr) and Recovery arm (EtBr removed at 72h)
   - Both OLS and mixed model per timepoint
   - Tests: when do genes respond? When do they recover? Do they recover at all?
   - The temporal dimension is what enables the Linear/Switch/Delayed pattern classification
:::

:::{prompt}
Run all 6 regressions: regress_dose, regress_cn, regress_cn_spline, regress_cn_quadratic, and regress_mt with scatter: "chunks:4", plus regress_longitudinal.
:::

---

# Step 3 — Merge Chunks

After all scattered regressions complete, merge the chunked RDS outputs into single files.

:::{note}
**What this step does:**
- Finds all `*_chunk*.rds` files from the 20 scattered regression pods
- Groups by base name (e.g. `ETBR_CN_lmer_chunk00.rds` through `..._chunk03.rds`)
- `rbind()`s all chunks into a single merged RDS per output (e.g. `ETBR_CN_lmer.rds`)
- Produces the same output files as if regressions had run without chunking
:::

:::{prompt}
Merge the chunked regression outputs using merge_chunks.
:::

---

# Step 4 — Classify Response Patterns

Combine all regression results and classify each significant gene into one of four temporal response trajectory patterns using likelihood-based templates.

:::{note}
**How pattern classification works:**

1. **Select significant genes** — Holm p < 0.05 on the CN mixed model (the paper's primary test). Expected: ~1,500 genes.

2. **Compute scaled fold-change** — For each gene, compute the mean log2 fold-change at 4 dose-binned timepoints representing different stages of the mtDNA depletion/recovery arc. These 4 values form the gene's "trajectory signature."

3. **Template matching via dnorm product likelihoods** — Each gene's 4-value trajectory is compared to 3 pattern templates:
   - **Linear** `(0, 0.5, 1, 0.5)` — Expression change is proportional to CN depletion throughout. Response scales linearly with dose, partially recovers when CN recovers. Canonical example: OXPHOS genes.
   - **Switch** `(0, 1, 1, 0)` — Expression changes rapidly once a CN threshold is crossed, maintains the change during treatment, but recovers quickly once CN begins to recover. Mirrors MT gene expression. Canonical example: MT-encoded genes themselves.
   - **Delayed** `(0, 0.25, 1, 1)` — Minimal initial response, increases over time even as CN stabilizes, and persists during recovery. Suggests secondary/stress signaling rather than direct CN sensing. Canonical example: PSAT1 (serine biosynthesis).
   - **None** — Gene is significant but trajectory doesn't clearly match any template (top two likelihoods lack >2x separation ratio).

4. **MitoCarta 3.0 annotation** — Each gene is labeled as:
   - **MT-encoded** (13 protein-coding genes on mitochondrial genome)
   - **MitoCarta** (nuclear-encoded mitochondrial proteins, ~1,100 genes)
   - **Nuclear** (all other nuclear genes)
   - MitoCarta pathway annotations include central dogma categories (mtRNA metabolism, Translation, Protein import, mtDNA maintenance) and OXPHOS complex assignments (I-V).

**Expected pattern breakdown from the paper:**
- ~30% Linear, ~25% Switch, ~20% Delayed, ~25% None
- Central dogma genes (POLG, TWINKLE, TFAM, MRPS/MRPL family) tend to be Linear and upregulated — compensatory response
- Glycolytic enzymes tend to be Delayed — metabolic reprogramming is a secondary response
:::

:::{prompt}
Classify response patterns using /pfs/promptable-public-test-data/arking-mtdna/Human.MitoCarta3.0.xlsx and /pfs/promptable-public-test-data/arking-mtdna/mart_export.txt with the regression results from the previous step.
:::

---

# Results — Statistical Summary

After the pipeline completes, ask for the core statistical findings. The specialist will read the output RDS files and report actual computed values.

:::{tip}
**Expected results from the paper:**
- ~11,800 genes tested after filtering (median > 49 in controls, ENSG prefix)
- ~1,500 significant genes (Holm p < 0.05 on CN mixed model)
- Pattern breakdown: roughly 30% Linear, 25% Switch, 20% Delayed, 25% None
- Quadratic/spline models significantly improve fit for ~500 genes (LRT P12/P13/P14 < 0.05)
- 0 PCA outliers expected (the data is clean), 2 samples excluded (RESV_NC, ETBR_NC1)
- High concordance between CN and MT models — both measure mitochondrial perturbation via different proxies
:::

:::{prompt}
Summarize the complete statistical results:
1. How many genes tested vs significant (Holm p<0.05) for each model (dose, CN linear, CN quadratic, CN spline, MT, longitudinal per timepoint)?
2. Pattern classification breakdown (Linear/Switch/Delayed/None) — counts and percentages
3. How many genes where quadratic or spline models significantly improved over linear (P12, P13, P14)?
4. How does significance compare between CN association and MT association (concordance)?
:::

---

# Results — Biological Interpretation

Dig into the biology. Even without figures and enrichment, the regression results and pattern classification reveal the paper's key findings about compensatory mitochondrial gene regulation and metabolic reprogramming.

:::{note}
**Key biological questions to explore:**
- Do MitoCarta central dogma genes (POLG, POLG2, TWINKLE/TWNK, TFAM, MRPS/MRPL ribosomal proteins, mtRNA metabolism enzymes) show compensatory upregulation as mtDNA-CN drops?
- Are nuclear-encoded OXPHOS genes (NDUFS/NDUFV for Complex I, SDHA/B for Complex II, UQCRC for Complex III, COX for Complex IV, ATP5 for Complex V) significantly downregulated?
- Which pattern do glycolytic enzymes fall into — Linear, Switch, or Delayed? A Delayed pattern would suggest metabolic reprogramming is a secondary response.
- Do the 13 MT-encoded genes show the strongest effects in the dose model? Does fold-change magnitude correlate with genomic position (HSP transcription gradient)?
- Which regression model (dose, CN, MT) identifies the most significant genes? Why might they differ?
:::

:::{prompt}
Provide biological interpretation of the regression and classification results:
1. MitoCarta pathway analysis: How do central dogma genes (mtRNA metabolism, Translation, mtDNA maintenance) respond? Are they compensatorily upregulated? List specific genes and their patterns.
2. OXPHOS genes: Are nuclear-encoded subunits significantly affected? Which complexes (I-V)? What patterns?
3. Glycolysis: What pattern classification did the 11 pathway enzymes (SLC2A1, HK1, GPI, PFKL, ALDOA, TPI1, GAPDH, PGK1, PGM1, ENO1, PKM) receive?
4. MT-encoded genes: Which show the strongest dose-response? Is there a positional gradient on the mitochondrial genome?
5. Model comparison: For genes where non-linear models (quadratic/spline) beat linear, what biological categories are enriched?
:::
