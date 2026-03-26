#!/usr/bin/env Rscript
# regress_mt.R — GE ~ Med_log2FC_MT (median MT gene fold change)
# From: chemcells_rnaseq_stat_1_regression.Rmd §2.3

suppressPackageStartupMessages(library(tidyverse))

shared_dir <- Sys.getenv("SHARED_DIR", "../../shared")
source(file.path(shared_dir, "regression_helpers.R"))

pheno_file <- Sys.getenv("PHENO_FILE")
rna_file   <- Sys.getenv("RNA_FILE")
output_dir <- Sys.getenv("OUTPUT_DIR", "/output")

cat("=== regress_mt ===\n")

dat_pheno <- readRDS(pheno_file)
dat_rna   <- readRDS(rna_file)

# Verify Med_log2FC_MT exists
stopifnot("Med_log2FC_MT" %in% names(dat_pheno))

pheno_sub <- dat_pheno %>% subset(Chemical == "ETBR" & Curve == "Dosage")
genes <- rownames(dat_rna)

cat("Genes:", length(genes), "| Samples:", nrow(pheno_sub), "\n")
cat("Med_log2FC_MT range:",
    round(range(pheno_sub$Med_log2FC_MT, na.rm = TRUE), 3), "\n")

# OLS with lane averaging
cat("Running TWAS.lm (Med_log2FC_MT)...\n")
ETBR_MT_lm <- TWAS.lm(genes, dat_rna, "Med_log2FC_MT", pheno_sub)

# Mixed model
cat("Running TWAS.lmer (Med_log2FC_MT)...\n")
ETBR_MT_lmer <- TWAS.lmer(genes, dat_rna, "Med_log2FC_MT", pheno_sub)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(ETBR_MT_lm, file.path(output_dir, "ETBR_MT_lm.rds"))
saveRDS(ETBR_MT_lmer, file.path(output_dir, "ETBR_MT_lmer.rds"))

n_sig <- sum(p.adjust(ETBR_MT_lmer$P, "holm") < 0.05, na.rm = TRUE)
cat("Significant genes (MT lmer, Holm 0.05):", n_sig, "\n")
cat("=== DONE ===\n")
