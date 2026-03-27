---
title: Xie et al. (2026) — Full Pipeline
description: Complete mtDNA-CN depletion transcriptomics analysis — normalization, 6 regressions, pattern classification, stratified GO enrichment with treemaps, 9 publication figures, and HTML report
---

# Setup — Add Specialist

Add the mtDNA Transcriptomics Specialist to your notebook's workflow team. This agent implements the complete Xie et al. (2026) analysis pipeline from the Arking Lab at Johns Hopkins — "Transcriptomic analysis of cells following decreased mitochondrial DNA-copy number reveals distinct transcriptional response patterns."

:::{prompt}
Add the mtdna-transcriptomics-specialist to the team.
:::

---

# Full Pipeline

Run the complete analysis in a single execution. The pipeline reproduces four Rmd notebooks (`data_1_read_data`, `stat_1_regression`, `stat_2_lon_track`, `stat_X_plots`) as 12 containerized tools with automatic dependency wiring.

:::{note}
**Pipeline overview (12 tools, 21 parallel pods, ~15 min total):**

**Phase 1 — Data Preparation**
- `prepare_expression_data` — Read raw counts (genes x 112 samples), parse SampleName to extract Chemical/Dose/Curve/Time/Rep, filter to EtBr-only (126 samples across Dosage/Treatment/Recovery arms), TMM-normalize with edgeR (median > 49 in controls, ENSG prefix), PCA outlier removal at 4 SD on first 10 PCs, compute Med_log2FC_MT (median log2 fold-change of 13 MT protein-coding genes vs condition-matched controls)

**Phase 2 — Gene-wise Regressions (5 tools x 4 gene chunks = 20 parallel pods)**
- `regress_dose` — log(CPM) ~ Dose (scaled 0-6) for all ~11,800 genes; both OLS (lane-averaged) and mixed model with (1|SampleName) random effect
- `regress_cn` — log(CPM) ~ deltaCT_Avg (mtDNA copy number); linear + quadratic, both OLS and mixed — 4 result sets
- `regress_cn_spline` — log(CPM) ~ ns(deltaCT_Avg, df=2) and ns(deltaCT_Avg, df=3) natural spline models
- `regress_cn_quadratic` — AIC model comparison: linear vs quadratic vs ns(2) vs ns(3); likelihood ratio tests (P12, P13, P14)
- `regress_mt` — log(CPM) ~ Med_log2FC_MT (MT gene expression proxy); OLS + mixed model
- Each regression tool runs with `scatter: "chunks:4"` — genes are split across 4 parallel pods for ~4x speedup

**Phase 3 — Longitudinal Analysis**
- `regress_longitudinal` — log(CPM) ~ Experiment (Control vs EtBr) at each of 7 timepoints (0, 24, 48, 72, 144, 168, 192h) in Treatment and Recovery arms; both OLS and mixed model per timepoint

**Phase 4 — Merge Chunks**
- `merge_chunks` — Finds all `*_chunk*.rds` files from the scattered regression outputs, groups by base name (e.g. `ETBR_CN_lmer_chunk00.rds` through `..._chunk03.rds` → `ETBR_CN_lmer.rds`), rbinds into single merged RDS files

**Phase 5 — Pattern Classification**
- `classify_response_patterns` — Combine all regression results, annotate genes with MitoCarta 3.0 (MT-encoded/MitoCarta/Nuclear), compute scaled fold-change at 4 dose-binned timepoints, classify into Linear/Switch/Delayed/None via dnorm product likelihoods at timepoints 2 and 4 with >2x ratio threshold

**Phase 6 — GO Enrichment**
- `run_go_enrichment` — 15 enrichGO calls (all-sig + up/down, then 4 patterns x 2 directions), 15-way goana table comparing all gene lists simultaneously, rrvgo semantic similarity treemaps for 4 key enrichment results

**Phase 7 — Figures**
- `generate_publication_figures` — 9 publication-quality PNG figures: CN dose-response (Fig 1A), MT gene dose grid (Fig 1B), MT position vs logFC by dose (Fig 2A), volcano with MitoCarta central dogma labels (Fig X), glycolysis enzyme curves (Fig 4, 11 genes: SLC2A1 through PKM), longitudinal CN time course (Fig 6A), MT gene longitudinal grid (Fig 6B), scaled lag comparison with PSAT1 (Fig 7)

**Phase 8 — Report**
- `generate_html_report` — Self-contained HTML with embedded base64 figures, pattern classification tables, model significance comparisons, MitoCarta pathway breakdown (central dogma, OXPHOS), and GO enrichment highlights

**Input files** (all in vault):
- `merged.csv` — RNA-seq count matrix (genes x samples)
- `AllSamples.Randomized.RNAseq.csv` — sample metadata bridge file (112 rows x 18 columns)
- `ChemicalProjectMatchUpDeltaCT.xlsx` — qPCR deltaCT measurements per sample
- `Human.MitoCarta3.0.xlsx` — MitoCarta 3.0 database (Sheet A: gene mappings, Sheet C: pathways)
- `mart_export.txt` — Ensembl gene annotations from BioMart (gene_id, symbol, chr, position, biotype)
:::

:::{prompt}
Run the full Xie et al. mtDNA depletion pipeline using:
- countFile: /pfs/promptable-public-test-data/arking-mtdna/merged.csv
- bridgeFile: /pfs/promptable-public-test-data/arking-mtdna/AllSamples.Randomized.RNAseq.csv
- deltaCTFile: /pfs/promptable-public-test-data/arking-mtdna/ChemicalProjectMatchUpDeltaCT.xlsx
- mitoCartaFile: /pfs/promptable-public-test-data/arking-mtdna/Human.MitoCarta3.0.xlsx
- annotationFile: /pfs/promptable-public-test-data/arking-mtdna/mart_export.txt

Run all steps: prepare_expression_data, all regressions with scatter: "chunks:4" (dose, CN, CN spline, CN quadratic, MT), regress_longitudinal, merge_chunks, classify_response_patterns, run_go_enrichment, generate_publication_figures, and generate_html_report.
:::

---

# Results — Statistical Summary

After the pipeline completes (~15 min), ask for the core statistical findings. The specialist will read the output RDS files and report actual computed values.

:::{tip}
**Expected results from the paper:**
- ~11,800 genes tested after filtering (median > 49 in controls, ENSG prefix)
- ~1,500 significant genes (Holm p < 0.05 on CN mixed model)
- Pattern breakdown: roughly 30% Linear, 25% Switch, 20% Delayed, 25% None
- Quadratic/spline models significantly improve fit for ~500 genes (LRT P12/P13/P14 < 0.05)
- 0 PCA outliers expected (the data is clean), 2 samples excluded (RESV_NC, ETBR_NC1)
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

Dig into the biology. The paper's key findings involve compensatory mitochondrial gene regulation, a glycolytic metabolic shift, and distinct temporal response patterns in different pathway families.

:::{note}
**Key biological questions from the paper:**
- Do MitoCarta central dogma genes (mtRNA metabolism, Translation, mtDNA maintenance) show compensatory upregulation?
- Are nuclear-encoded OXPHOS genes downregulated while MT-encoded OXPHOS genes show different patterns?
- Is there a coordinated glycolytic shift across all 11 pathway enzymes?
- Does MT gene position on the mitochondrial genome correlate with fold-change magnitude (HSP transcription gradient)?
- Which GO terms distinguish Linear vs Switch vs Delayed genes?
:::

:::{prompt}
Provide biological interpretation:
1. MitoCarta pathway analysis: How do central dogma genes (mtRNA metabolism, Translation, mtDNA maintenance) respond? Are they compensatorily upregulated? List specific genes.
2. OXPHOS genes: Are nuclear-encoded subunits significantly affected? Which complexes (I-V)?
3. Glycolysis: Do all 11 enzymes (SLC2A1, HK1, GPI, PFKL, ALDOA, TPI1, GAPDH, PGK1, PGM1, ENO1, PKM) show coordinated upregulation?
4. Enrichment highlights: What are the top GO terms for each direction (up/down) and pattern (Linear/Switch/Delayed)?
5. What do the rrvgo treemaps show — which GO clusters emerge for delayed-response genes?
:::

---

# Results — Figures

Review the 9 publication figures. Each figure tells a specific part of the biological story.

:::{note}
**Figure guide:**
- **Fig 1A** — CN dose-response: deltaCT should decrease monotonically with EtBr dose (0-150 ng/mL), confirming mtDNA depletion
- **Fig 1B** — MT gene grid: all 13 protein-coding genes should show dose-dependent downregulation, but magnitude varies by genomic position
- **Fig 2A** — MT position vs logFC: genes closer to the heavy-strand promoter (HSP, position ~550) have larger fold-changes at high doses — evidence for a transcription gradient
- **Fig X** — Volcano: MitoCarta central dogma genes labeled, showing compensatory upregulation (upper-left quadrant)
- **Fig 4** — Glycolysis: 11 enzymes showing coordinated upregulation — metabolic shift from OXPHOS to glycolysis
- **Fig 6A** — Longitudinal CN: deltaCT recovery after EtBr removal (Treatment arm depleted at 48h, Recovery arm shows partial return by 192h)
- **Fig 6B** — MT gene longitudinal: expression recovery lags behind CN recovery
- **Fig 7** — The three response patterns overlaid with PSAT1 as the canonical Delayed gene — demonstrates temporal separation
:::

:::{prompt}
Describe each of the 9 figures generated. For each: what does it show biologically, and does it match the expected pattern from the paper? Are there any surprising features in the dose-response curves or longitudinal trajectories?
:::

---

# Results — GO Enrichment Detail

The enrichment analysis runs 15 separate enrichGO analyses plus a 15-way goana comparison. This is the richest part of the analysis.

:::{tip}
**What to look for in the enrichment results:**
- **Upregulated genes (Estimate < 0, higher expression with lower CN):** amino acid metabolism, tRNA aminoacylation, nutrient sensing — consistent with compensatory/stress response
- **Downregulated genes (Estimate > 0, lower expression with lower CN):** cholesterol metabolism, lipid processing — metabolic shutdown
- **Delayed-up genes specifically:** "response to nutrient levels", "amino acid import across plasma membrane" — lagged stress signaling
- **Delayed-down genes:** "ribonucleoprotein complex biogenesis" — structural/growth programs suppressed late
- The **goana 15-way table** lets you compare any two pattern categories head-to-head for the same GO term
:::

:::{prompt}
Summarize the GO enrichment results in detail:
1. How many significant GO terms (BP, BH p<0.05) for each of the 15 enrichGO analyses?
2. What are the top 5 enriched terms for upregulated vs downregulated genes overall?
3. Which terms distinguish Delayed genes from Linear/Switch genes?
4. What do the rrvgo treemaps reveal — which semantic clusters emerge?
5. From the 15-way goana table, which GO terms show the most dramatic pattern-specific enrichment?
:::
