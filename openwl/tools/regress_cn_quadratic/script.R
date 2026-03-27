#!/usr/bin/env Rscript
# regress_cn_quadratic.R — AIC model comparison (linear vs quadratic vs spline)
# From: chemcells_rnaseq_stat_1_regression.Rmd §2.2

suppressPackageStartupMessages(library(tidyverse))

shared_dir <- Sys.getenv("SHARED_DIR", "../../shared")
source(file.path(shared_dir, "regression_helpers.R"))

pheno_file <- Sys.getenv("PHENO_FILE")
rna_file   <- Sys.getenv("RNA_FILE")
mtdna_file <- Sys.getenv("MTDNA_FILE")
output_dir <- Sys.getenv("OUTPUT_DIR", "/output")

cat("=== regress_cn_quadratic (ANOVA model comparison) ===\n")

dat_pheno <- readRDS(pheno_file)
dat_rna   <- readRDS(rna_file)
dat_CN    <- readRDS(mtdna_file)

# Join deltaCT_Avg
dat_CN$PO <- factor(dat_CN$PO)
dat_pheno$deltaCT_Avg <- sapply(dat_pheno$SampleName, function(x) {
  vals <- dat_CN$deltaCT_Avg[which(dat_CN$Sample_Name == x & !is.na(dat_CN$deltaCT_Avg))]
  if (length(vals) == 0) return(NA)
  return(vals[1])
})

pheno_sub <- dat_pheno %>% subset(Chemical == "ETBR" & Curve == "Dosage")
genes <- get_gene_chunk(rownames(dat_rna))
sfx <- chunk_suffix()

cat("Genes:", length(genes), "| Samples:", nrow(pheno_sub), "\n")

# ANOVA model comparison
cat("Running TWAS.anova (deltaCT_Avg)...\n")
ETBR_CN_anova <- TWAS.anova(genes, dat_rna, "deltaCT_Avg", pheno_sub)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(ETBR_CN_anova, file.path(output_dir, paste0("ETBR_CN_anova", sfx, ".rds")))

n_sig_quad <- sum(p.adjust(ETBR_CN_anova$P12, "holm") < 0.05, na.rm = TRUE)
n_sig_sp2  <- sum(p.adjust(ETBR_CN_anova$P13, "holm") < 0.05, na.rm = TRUE)
n_sig_sp3  <- sum(p.adjust(ETBR_CN_anova$P14, "holm") < 0.05, na.rm = TRUE)
cat("Significant (Holm 0.05): quadratic =", n_sig_quad,
    "| spline(2) =", n_sig_sp2, "| spline(3) =", n_sig_sp3, "\n")
cat("=== DONE ===\n")
