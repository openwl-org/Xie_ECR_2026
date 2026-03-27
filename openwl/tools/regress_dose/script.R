#!/usr/bin/env Rscript
# regress_dose.R — GE ~ Dose linear regression
# From: chemcells_rnaseq_stat_1_regression.Rmd §2.1

suppressPackageStartupMessages(library(tidyverse))

shared_dir <- Sys.getenv("SHARED_DIR", "../../shared")
source(file.path(shared_dir, "regression_helpers.R"))

pheno_file <- Sys.getenv("PHENO_FILE")
rna_file   <- Sys.getenv("RNA_FILE")
output_dir <- Sys.getenv("OUTPUT_DIR", "/output")

cat("=== regress_dose ===\n")

dat_pheno <- readRDS(pheno_file)
dat_rna   <- readRDS(rna_file)

# Subset to ETBR dosage samples
pheno_sub <- dat_pheno %>% subset(Chemical == "ETBR" & Curve == "Dosage")
genes <- get_gene_chunk(rownames(dat_rna))
sfx <- chunk_suffix()

cat("Genes:", length(genes), "| Samples:", nrow(pheno_sub), "\n")

# OLS with lane averaging
cat("Running TWAS.lm (Dose)...\n")
ETBR_dose_lm <- TWAS.lm(genes, dat_rna, "Dose", pheno_sub)

# Mixed model
cat("Running TWAS.lmer (Dose)...\n")
ETBR_dose_lmer <- TWAS.lmer(genes, dat_rna, "Dose", pheno_sub)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(ETBR_dose_lm, file.path(output_dir, paste0("ETBR_dose_lm", sfx, ".rds")))
saveRDS(ETBR_dose_lmer, file.path(output_dir, paste0("ETBR_dose_lmer", sfx, ".rds")))

n_sig_lm   <- sum(p.adjust(ETBR_dose_lm$P, "holm") < 0.05, na.rm = TRUE)
n_sig_lmer <- sum(p.adjust(ETBR_dose_lmer$P, "holm") < 0.05, na.rm = TRUE)
cat("Significant genes (Holm 0.05): lm =", n_sig_lm, "| lmer =", n_sig_lmer, "\n")
cat("=== DONE ===\n")
