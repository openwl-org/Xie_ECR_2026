#!/usr/bin/env Rscript
# regress_cn.R — GE ~ deltaCT_Avg (linear + quadratic)
# From: chemcells_rnaseq_stat_1_regression.Rmd §2.2

suppressPackageStartupMessages(library(tidyverse))

shared_dir <- Sys.getenv("SHARED_DIR", "../../shared")
source(file.path(shared_dir, "regression_helpers.R"))

pheno_file <- Sys.getenv("PHENO_FILE")
rna_file   <- Sys.getenv("RNA_FILE")
mtdna_file <- Sys.getenv("MTDNA_FILE")
output_dir <- Sys.getenv("OUTPUT_DIR", "/output")

cat("=== regress_cn ===\n")

dat_pheno <- readRDS(pheno_file)
dat_rna   <- readRDS(rna_file)
dat_CN    <- readRDS(mtdna_file)

# Join deltaCT_Avg to pheno
dat_CN$PO <- factor(dat_CN$PO)
dat_pheno$deltaCT_Avg <- sapply(dat_pheno$SampleName, function(x) {
  vals <- dat_CN$deltaCT_Avg[which(dat_CN$Sample_Name == x & !is.na(dat_CN$deltaCT_Avg))]
  if (length(vals) == 0) return(NA)
  return(vals[1])
})

pheno_sub <- dat_pheno %>% subset(Chemical == "ETBR" & Curve == "Dosage")
genes <- rownames(dat_rna)

cat("Genes:", length(genes), "| Samples:", nrow(pheno_sub), "\n")

# Linear OLS
cat("Running TWAS.lm (deltaCT_Avg)...\n")
ETBR_CN_lm <- TWAS.lm(genes, dat_rna, "deltaCT_Avg", pheno_sub)

# Quadratic OLS
cat("Running TWAS.lm2 (deltaCT_Avg)...\n")
ETBR_CN_lm2 <- TWAS.lm2(genes, dat_rna, "deltaCT_Avg", pheno_sub)

# Linear mixed model
cat("Running TWAS.lmer (deltaCT_Avg)...\n")
ETBR_CN_lmer <- TWAS.lmer(genes, dat_rna, "deltaCT_Avg", pheno_sub)

# Quadratic mixed model
cat("Running TWAS.lmer2 (deltaCT_Avg)...\n")
ETBR_CN_lmer2 <- TWAS.lmer2(genes, dat_rna, "deltaCT_Avg", pheno_sub)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(ETBR_CN_lm, file.path(output_dir, "ETBR_CN_lm.rds"))
saveRDS(ETBR_CN_lm2, file.path(output_dir, "ETBR_CN_lm2.rds"))
saveRDS(ETBR_CN_lmer, file.path(output_dir, "ETBR_CN_lmer.rds"))
saveRDS(ETBR_CN_lmer2, file.path(output_dir, "ETBR_CN_lmer2.rds"))

n_sig <- sum(p.adjust(ETBR_CN_lmer$P, "holm") < 0.05, na.rm = TRUE)
cat("Significant genes (CN lmer, Holm 0.05):", n_sig, "\n")
cat("=== DONE ===\n")
