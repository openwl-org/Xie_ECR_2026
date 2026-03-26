#!/usr/bin/env Rscript
# regress_longitudinal.R — GE ~ Experiment (Control vs EtBr) per timepoint
# From: chemcells_rnaseq_stat_1_regression.Rmd §2.4

suppressPackageStartupMessages(library(tidyverse))

shared_dir <- Sys.getenv("SHARED_DIR", "../../shared")
source(file.path(shared_dir, "regression_helpers.R"))

pheno_file <- Sys.getenv("PHENO_FILE")
rna_file   <- Sys.getenv("RNA_FILE")
output_dir <- Sys.getenv("OUTPUT_DIR", "/output")

cat("=== regress_longitudinal ===\n")

dat_pheno <- readRDS(pheno_file)
dat_rna   <- readRDS(rna_file)

# Longitudinal samples only (Treatment + Recovery, not Dosage)
dat_Lon <- dat_pheno %>% subset(Curve != "Dosage")

# Create binary Experiment factor: Control (Dose=0) vs EtBr (Dose=6)
dat_Lon$Experiment <- factor(
  c("Control", "EtBr")[as.numeric(dat_Lon$Dose == 6) + 1],
  levels = c("Control", "EtBr")
)

genes <- rownames(dat_rna)
timepoints <- unique(dat_Lon$Time)
cat("Genes:", length(genes), "| Timepoints:", paste(timepoints, collapse = ", "), "\n")

# Initialize output data frames
ETBR_LON_lm <- data.frame(Gene = character(), Estimate = numeric(),
                            SE = numeric(), P = numeric(), Time = numeric())
ETBR_LON_lmer <- data.frame(Gene = character(), Estimate = numeric(),
                              SE = numeric(), P = numeric(), Time = numeric())

for (tp in timepoints) {
  cat("\nTimepoint:", tp, "\n")
  dat_tp <- dat_Lon %>% subset(Time == tp)

  # OLS
  cat("  Running TWAS.lm (Experiment)...\n")
  short_lm <- TWAS.lm(genes, dat_rna, "Experiment", dat_tp)
  short_lm$Time <- tp
  ETBR_LON_lm <- rbind(ETBR_LON_lm, short_lm)

  # Mixed model
  cat("  Running TWAS.lmer (Experiment)...\n")
  short_lmer <- TWAS.lmer(genes, dat_rna, "Experiment", dat_tp)
  short_lmer$Time <- tp
  ETBR_LON_lmer <- rbind(ETBR_LON_lmer, short_lmer)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(ETBR_LON_lm, file.path(output_dir, "ETBR_LON_lm.rds"))
saveRDS(ETBR_LON_lmer, file.path(output_dir, "ETBR_LON_lmer.rds"))

cat("\n=== Summary ===\n")
for (tp in timepoints) {
  n_sig <- sum(p.adjust(ETBR_LON_lmer$P[ETBR_LON_lmer$Time == tp], "holm") < 0.05,
               na.rm = TRUE)
  cat("  Time", tp, ": sig genes =", n_sig, "\n")
}
cat("=== DONE ===\n")
