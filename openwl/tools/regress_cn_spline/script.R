#!/usr/bin/env Rscript
# regress_cn_spline.R — GE ~ ns(deltaCT_Avg, df=2) and ns(deltaCT_Avg, df=3)
# From: chemcells_rnaseq_stat_1_regression.Rmd §2.2

suppressPackageStartupMessages(library(tidyverse))

shared_dir <- Sys.getenv("SHARED_DIR", "../../shared")
source(file.path(shared_dir, "regression_helpers.R"))

pheno_file <- Sys.getenv("PHENO_FILE")
rna_file   <- Sys.getenv("RNA_FILE")
mtdna_file <- Sys.getenv("MTDNA_FILE")
output_dir <- Sys.getenv("OUTPUT_DIR", "/output")

cat("=== regress_cn_spline ===\n")

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

# Natural spline df=2
cat("Running TWAS.lmer3 (deltaCT_Avg, df=2)...\n")
ETBR_CN_ns2 <- TWAS.lmer3(genes, dat_rna, "deltaCT_Avg", pheno_sub, 2)

# Natural spline df=3
cat("Running TWAS.lmer3 (deltaCT_Avg, df=3)...\n")
ETBR_CN_ns3 <- TWAS.lmer3(genes, dat_rna, "deltaCT_Avg", pheno_sub, 3)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(ETBR_CN_ns2, file.path(output_dir, paste0("ETBR_CN_ns2", sfx, ".rds")))
saveRDS(ETBR_CN_ns3, file.path(output_dir, paste0("ETBR_CN_ns3", sfx, ".rds")))

cat("=== DONE ===\n")
