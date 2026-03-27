---
title: Xie et al. (2026) — Full Pipeline + Custom R Exploration
description: Complete mtDNA-CN depletion pipeline (11 tools) followed by interactive R exploration — write and execute custom R code against the pipeline outputs to answer biological questions beyond what the fixed tools compute
---

# Setup — Add Specialists

Add both the mtDNA Transcriptomics Specialist (runs the fixed pipeline) and the Coding Agent (writes and executes custom R code) to your notebook's workflow team.

:::{prompt}
Add the mtdna-transcriptomics-specialist to the team.
:::

:::{prompt}
Add the coding-agent to the team.
:::

---

# Full Pipeline

Run the complete Xie et al. (2026) analysis pipeline — all 12 tools from data preparation through HTML report generation.

:::{note}
**Pipeline overview (11 tools, 28 parallel pods, ~12 min total):**

**Phase 1 — Data Preparation**
- `prepare_expression_data` — Read raw counts (genes x 112 samples), parse SampleName to extract Chemical/Dose/Curve/Time/Rep, filter to EtBr-only (126 samples across Dosage/Treatment/Recovery arms), TMM-normalize with edgeR (median > 49 in controls, ENSG prefix), PCA outlier removal at 4 SD on first 10 PCs, compute Med_log2FC_MT (median log2 fold-change of 13 MT protein-coding genes vs condition-matched controls)

**Phase 2 — Gene-wise Regressions (6 tools x 4 gene chunks = 24 parallel pods + auto-gather)**
- `regress_dose` — log(CPM) ~ Dose (scaled 0-6); OLS + mixed model
- `regress_cn` — log(CPM) ~ deltaCT_Avg; linear + quadratic, OLS + mixed (4 result sets)
- `regress_cn_spline` — log(CPM) ~ ns(deltaCT_Avg, df=2) and df=3
- `regress_cn_quadratic` — AIC model comparison + likelihood ratio tests (P12, P13, P14)
- `regress_mt` — log(CPM) ~ Med_log2FC_MT; OLS + mixed model
- `regress_longitudinal` — log(CPM) ~ Experiment at each of 7 timepoints (0-192h)
- Each regression runs with `scatter: "chunks:4"` — genes split across 4 parallel pods
- `merge_chunks` auto-injected after each scattered tool (no need to plan it)

**Phase 3 — Pattern Classification**
- `classify_response_patterns` — Classify genes into Linear/Switch/Delayed/None via dnorm likelihoods

**Phase 4 — GO Enrichment**
- `run_go_enrichment` — 15 enrichGO analyses + 15-way goana + rrvgo treemaps

**Phase 5 — Figures**
- `generate_publication_figures` — 9 publication-quality PNG figures

**Phase 6 — Report**
- `generate_html_report` — Self-contained HTML with embedded figures and tables

**Input files** (all in vault):
- `merged.csv` — RNA-seq count matrix
- `AllSamples.Randomized.RNAseq.csv` — sample metadata bridge file
- `ChemicalProjectMatchUpDeltaCT.xlsx` — qPCR deltaCT measurements
- `Human.MitoCarta3.0.xlsx` — MitoCarta 3.0 database
- `mart_export.txt` — Ensembl gene annotations
:::

:::{prompt}
Run the full Xie et al. mtDNA depletion pipeline using:
- countFile: /pfs/promptable-public-test-data/arking-mtdna/merged.csv
- bridgeFile: /pfs/promptable-public-test-data/arking-mtdna/AllSamples.Randomized.RNAseq.csv
- deltaCTFile: /pfs/promptable-public-test-data/arking-mtdna/ChemicalProjectMatchUpDeltaCT.xlsx
- mitoCartaFile: /pfs/promptable-public-test-data/arking-mtdna/Human.MitoCarta3.0.xlsx
- annotationFile: /pfs/promptable-public-test-data/arking-mtdna/mart_export.txt

Run all steps: prepare_expression_data, all regressions with scatter: "chunks:4" (dose, CN, CN spline, CN quadratic, MT, longitudinal), classify_response_patterns, run_go_enrichment, generate_publication_figures, and generate_html_report. (merge_chunks is auto-injected after each scattered regression — do NOT include it in the plan.)
:::

---

# Explore — OXPHOS Complex Breakdown

The fixed pipeline reports overall OXPHOS significance, but doesn't break it down by electron transport chain complex. Write custom R code to compute per-complex statistics.

:::{note}
**What this custom analysis does:**
- Loads the CN mixed-model regression results (RDS) and MitoCarta 3.0 annotations
- Filters to genes in MitoCarta's "OXPHOS" pathway
- Groups by ETC complex (I: NDUFS/NDUFV/NDUFA/NDUFB; II: SDH; III: UQCR; IV: COX; V: ATP5)
- For each complex: counts significant genes, computes median effect size, identifies most/least affected subunits
- Outputs a summary table comparing complex-level responses

**Biological question:** Are all 5 complexes equally affected by mtDNA-CN depletion, or are some (e.g., Complex I with its many nuclear-encoded subunits) more sensitive than others?
:::

:::{prompt}
Using the coding agent, write and run R code that:
1. Reads the CN mixed-model regression RDS from the pipeline outputs
2. Reads MitoCarta 3.0 from /pfs/promptable-public-test-data/arking-mtdna/Human.MitoCarta3.0.xlsx
3. Filters to OXPHOS genes, groups by ETC complex (I-V)
4. For each complex: count significant genes (Holm p<0.05), median estimate, min/max genes
5. Print a summary table

Use image gcr.io/promptable-mvp/mtdna-transcriptomics-toolbox:v1.0.0 which has all the R packages needed (edgeR, lme4, openxlsx, etc.).
:::

---

# Explore — Amino Acid Metabolism Deep Dive

The GO enrichment flags "amino acid metabolic process" as enriched in Delayed-up genes. Write custom R to identify exactly which amino acid biosynthesis pathways are driving this signal.

:::{note}
**What this custom analysis does:**
- Loads the pattern classification results + GO enrichment results
- Identifies all Delayed-up genes annotated with amino acid metabolism GO terms
- Cross-references with KEGG pathways for specific amino acid biosynthesis routes
- Checks if serine/glycine biosynthesis genes (PHGDH, PSAT1, PSPH, SHMT1/2) are ALL Delayed — this would indicate a coordinated pathway-level response
- Computes the average lag score for amino acid vs non-amino acid Delayed genes

**Biological question:** Is the amino acid metabolic response concentrated in specific biosynthesis pathways (serine/glycine one-carbon metabolism), or spread across many amino acid pathways?
:::

:::{prompt}
Using the coding agent, write and run R code that:
1. Reads the pattern classification RDS from pipeline outputs
2. Filters to Delayed-up genes
3. Checks which amino acid biosynthesis genes are in the Delayed category — specifically look for PHGDH, PSAT1, PSPH, SHMT1, SHMT2, CBS, CTH, ASNS, ATF4 (ISR pathway genes)
4. For each gene found: print its pattern, estimate, p-value, and which pathway it belongs to
5. Test whether serine/glycine biosynthesis is specifically enriched in Delayed vs Linear/Switch

Use image gcr.io/promptable-mvp/mtdna-transcriptomics-toolbox:v1.0.0.
:::

---

# Explore — Ribosomal Protein Compensation

The paper notes compensatory upregulation of mitochondrial translation machinery. Write custom R to compare mitochondrial ribosomal proteins (MRPS/MRPL) with cytoplasmic ribosomal proteins (RPS/RPL) — are they regulated in opposite directions?

:::{note}
**What this custom analysis does:**
- Loads regression results
- Separates mitochondrial ribosomal proteins (MRPS1-MRPS39, MRPL1-MRPL58) from cytoplasmic ribosomal proteins (RPS2-RPS30, RPL3-RPL41)
- Compares direction and magnitude of CN association between the two groups
- Tests whether mitoribosomal proteins are significantly more upregulated than cytoribosmal proteins
- Generates a paired comparison (box plot or summary table)

**Biological question:** If MRPS/MRPL genes are upregulated (compensatory mitochondrial translation) while RPS/RPL genes are unaffected or downregulated (general translation suppression), this would demonstrate a targeted compensatory response rather than a general translation change.
:::

:::{prompt}
Using the coding agent, write and run R code that:
1. Reads the CN mixed-model regression results from pipeline outputs
2. Identifies mitochondrial ribosomal proteins (grep "^MRPS|^MRPL" on gene symbols) and cytoplasmic ribosomal proteins (grep "^RPS[0-9]|^RPL[0-9]" on gene symbols)
3. For each group: count significant genes, median estimate (direction), mean -log10(p)
4. Run a Wilcoxon test comparing estimates between the two groups
5. Print summary: are mitoribosomal proteins significantly more upregulated than cytoribosomal proteins?

Use image gcr.io/promptable-mvp/mtdna-transcriptomics-toolbox:v1.0.0.
:::

---

# Explore — Temporal Response Heatmap

Create a custom heatmap showing how the top 50 genes (by CN p-value) change across all 7 longitudinal timepoints — visualize the temporal dynamics that the pattern classification summarizes.

:::{note}
**What this custom analysis does:**
- Loads longitudinal regression results (7 timepoints x OLS estimates per gene)
- Selects top 50 genes by CN mixed-model significance
- Builds a matrix: genes x timepoints, values = estimate (fold-change direction)
- Clusters genes by temporal trajectory (hierarchical clustering)
- Generates a heatmap with gene symbols, timepoint labels, and pattern annotation color bar
- Saves as PNG

**Biological question:** Do genes with the same pattern classification (Linear/Switch/Delayed) cluster together in the heatmap? Are there sub-patterns within each class that the 3-template classification misses?
:::

:::{prompt}
Using the coding agent, write and run R code that:
1. Reads the longitudinal regression RDS files and the CN mixed-model results from pipeline outputs
2. Selects the top 50 genes by CN significance (smallest Holm p-values)
3. Builds a timepoint x gene matrix of fold-change estimates at 0h, 24h, 48h, 72h, 144h, 168h, 192h
4. Generates a heatmap (use pheatmap or ComplexHeatmap if available, otherwise base R heatmap) with:
   - Row clustering (genes grouped by temporal trajectory)
   - Column order = chronological timepoints
   - Color annotation sidebar showing pattern classification (Linear=blue, Switch=green, Delayed=red, None=gray)
5. Save as /output/temporal_heatmap_top50.png

Use image gcr.io/promptable-mvp/mtdna-transcriptomics-toolbox:v1.0.0.
:::

---

# Explore — Custom Gene Set Analysis

Run a custom analysis on any gene set of interest. Replace the gene list below with your own.

:::{tip}
**Ideas for custom gene sets:**
- **ISR (Integrated Stress Response):** ATF4, DDIT3/CHOP, ASNS, TRIB3, SLC7A11, MTHFD2, PSAT1 — does the ISR activate during mtDNA depletion?
- **Autophagy/mitophagy:** PINK1, PRKN/Parkin, BNIP3, BNIP3L, FUNDC1, ATG5, ATG7, BECN1 — is damaged mitochondria being cleared?
- **TCA cycle:** CS, ACO2, IDH2, IDH3A, OGDH, SUCLA2, SUCLG1, FH, MDH2 — is TCA cycle downregulated alongside OXPHOS?
- **Ferroptosis:** GPX4, SLC7A11, ACSL4, LPCAT3, HMOX1 — does mtDNA depletion sensitize to ferroptosis?
- **mtDNA maintenance:** POLG, POLG2, TWNK, SSBP1, TFAM, MGME1, DNA2, RNASEH1 — which maintenance genes compensate most?
:::

:::{prompt}
Using the coding agent, write and run R code that:
1. Reads the CN mixed-model regression and pattern classification results from pipeline outputs
2. For the following gene set (Integrated Stress Response): ATF4, DDIT3, ASNS, TRIB3, SLC7A11, MTHFD2, PSAT1, SLC7A5, EIF2AK4
3. For each gene: print symbol, CN estimate, Holm p-value, pattern classification, and whether it's in MitoCarta
4. Summarize: are these genes coordinately regulated? What pattern do they share? Is the ISR activated during mtDNA depletion?

Use image gcr.io/promptable-mvp/mtdna-transcriptomics-toolbox:v1.0.0.
:::
