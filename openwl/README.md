<p align="center">
  <img src="icon.png" width="128" alt="mtDNA Transcriptomics Specialist">
</p>

<h3 align="center">mtDNA Transcriptomics Specialist</h3>

<p align="center">
  An <a href="https://openwl.org">OpenWL</a> agent for the Xie et al. (2026) EtBr-induced mtDNA copy number depletion study.<br>
  TMM-normalized RNA-seq analysis, gene-wise regressions, response pattern classification, GO enrichment, and publication figures.
</p>

---

## Directory Structure

```
openwl/
├── agent.yaml          # Agent definition: metadata, intelligence instructions, tool registry
├── compile.py          # Compiler: agent.yaml + tools/ + shared/ → single JSON spec for deployment
├── icon.png            # Agent avatar
├── compiled/           # Build output (gitignored) — compiled JSON specs ready for deployment
├── shared/             # R helper libraries injected into every tool at compile time
│   ├── data_helpers.R          # Sample metadata parsing, TMM normalization utilities
│   └── regression_helpers.R    # Gene chunking (get_gene_chunk), output suffixing (chunk_suffix)
├── templates/          # Promptable notebook templates — importable multi-step workflows
│   ├── scatter-test-minimal.md           # Quick demo: 2 regressions × 2 chunks + figures (~5 min)
│   ├── regression-only.md                # Statistical core: all 6 regressions + classification
│   ├── full-pipeline.md                  # Complete pipeline: 12 tools from prep to HTML report
│   └── full-pipeline-with-exploration.md # Full pipeline + custom R exploration via coding agent
└── tools/              # One directory per tool — each has tool.yaml (params/outputs) + script.R
    ├── prepare_expression_data/    # Read counts, parse metadata, TMM-normalize, PCA outliers
    ├── regress_dose/               # log(CPM) ~ Dose; OLS + mixed model
    ├── regress_cn/                 # log(CPM) ~ deltaCT_Avg; linear + quadratic
    ├── regress_cn_spline/          # log(CPM) ~ ns(deltaCT_Avg, df=2/3)
    ├── regress_cn_quadratic/       # AIC model comparison across CN models
    ├── regress_mt/                 # log(CPM) ~ Med_log2FC_MT
    ├── regress_longitudinal/       # log(CPM) ~ Experiment per timepoint
    ├── merge_chunks/               # Gather: rbind scattered *_chunk*.rds → merged RDS
    ├── classify_response_patterns/ # Linear/Switch/Delayed/None via likelihood templates
    ├── run_go_enrichment/          # Stratified GO enrichment + rrvgo treemaps
    ├── generate_publication_figures/ # 9 publication-quality PNG figures
    └── generate_html_report/       # Self-contained HTML with embedded figures + tables
```

## Pipeline

```
prepare_expression_data
         │
         ├──→ regress_dose ──────┐
         ├──→ regress_cn ────────┤
         ├──→ regress_cn_spline ─┤  (each scattered across 4 parallel pods)
         ├──→ regress_cn_quadratic┤
         ├──→ regress_mt ────────┤
         └──→ regress_longitudinal┘
                     │
              merge_chunks
                     │
        classify_response_patterns
                     │
            ┌────────┴────────┐
    run_go_enrichment   generate_publication_figures
            └────────┬────────┘
         generate_html_report
```

## Build & Deploy

```bash
# Compile agent.yaml + tools/ + shared/ into a single JSON spec
python3 openwl/compile.py

# Deploy to Promptable
ptbl-admin agent update mtdna-transcriptomics-specialist openwl/compiled/mtdna-transcriptomics-specialist.owl.json
```

The compiler embeds shared R helpers into each tool's bash wrapper, reads the git SHA as the version, and links back to the exact commit on GitHub. Always commit before compiling so the source URL is accurate.
