#!/usr/bin/env Rscript
# classify_response_patterns.R
# From: chemcells_rnaseq_stat_2_lon_track.Rmd
# Likelihood-based classification: Linear / Switch / Delayed / None

suppressPackageStartupMessages({
  library(readxl)
  library(lme4)
  library(tidyverse)
})

shared_dir <- Sys.getenv("SHARED_DIR", "../../shared")
source(file.path(shared_dir, "data_helpers.R"))
source(file.path(shared_dir, "regression_helpers.R"))

# Input files
pheno_file      <- Sys.getenv("PHENO_FILE")
rna_file        <- Sys.getenv("RNA_FILE")
mtdna_file      <- Sys.getenv("MTDNA_FILE")
annotation_file <- Sys.getenv("ANNOTATION_FILE")
mitocarta_file  <- Sys.getenv("MITO_CARTA_FILE")
dose_lm_file    <- Sys.getenv("DOSE_LM_FILE")
dose_lmer_file  <- Sys.getenv("DOSE_LMER_FILE")
cn_lm_file      <- Sys.getenv("CN_LM_FILE")
cn_lmer_file    <- Sys.getenv("CN_LMER_FILE")
mt_lm_file      <- Sys.getenv("MT_LM_FILE")
mt_lmer_file    <- Sys.getenv("MT_LMER_FILE")
output_dir      <- Sys.getenv("OUTPUT_DIR", "/output")

cat("=== classify_response_patterns ===\n")

# ============================================================
# 1. Load data
# ============================================================

dat_pheno <- readRDS(pheno_file)
dat_rna   <- readRDS(rna_file)
dat_CN    <- readRDS(mtdna_file)
G_list    <- load_gene_list(annotation_file, filter_genes = rownames(dat_rna))
mc        <- load_mitocarta(mitocarta_file)

ETBR_dose_lm   <- readRDS(dose_lm_file)
ETBR_dose_lmer <- readRDS(dose_lmer_file)
ETBR_CN_lm     <- readRDS(cn_lm_file)
ETBR_CN_lmer   <- readRDS(cn_lmer_file)
ETBR_MT_lm     <- readRDS(mt_lm_file)
ETBR_MT_lmer   <- readRDS(mt_lmer_file)

# Join deltaCT_Avg
dat_CN$PO <- factor(dat_CN$PO)
dat_pheno$deltaCT_Avg <- sapply(dat_pheno$SampleName, function(x) {
  vals <- dat_CN$deltaCT_Avg[which(dat_CN$Sample_Name == x & !is.na(dat_CN$deltaCT_Avg))]
  if (length(vals) == 0) return(NA)
  return(vals[1])
})

# ============================================================
# 2. Combine all regression results
# ============================================================

dat_all <- ETBR_dose_lm %>% dplyr::select(Gene) %>%
  left_join(G_list, by = c("Gene" = "ensembl_gene_id")) %>%
  mutate(
    GeneType = case_when(
      chromosome_name == "MT" ~ "3 MT encoded",
      Gene %in% mc$genes$EnsemblGeneID_mapping_version_20200130 ~ "2 MitoCarta",
      TRUE ~ "1 Nuclear"
    )
  ) %>%
  subset(gene_type %in% c("protein_coding", "Mt_rRNA", "Mt_tRNA",
                           "lncRNA", "miRNA", "misc_RNA", "snRNA")) %>%
  left_join(ETBR_dose_lm, by = "Gene") %>%
  rename(dose_lm_b = Estimate, dose_lm_SE = SE, dose_lm_P = P) %>%
  left_join(ETBR_dose_lmer, by = "Gene") %>%
  rename(dose_lmer_b = Estimate, dose_lmer_SE = SE, dose_lmer_P = P) %>%
  left_join(ETBR_CN_lm, by = "Gene") %>%
  rename(CN_lm_b = Estimate, CN_lm_SE = SE, CN_lm_P = P) %>%
  left_join(ETBR_CN_lmer, by = "Gene") %>%
  rename(CN_lmer_b = Estimate, CN_lmer_SE = SE, CN_lmer_P = P) %>%
  left_join(ETBR_MT_lm, by = "Gene") %>%
  rename(MT_lm_b = Estimate, MT_lm_SE = SE, MT_lm_P = P) %>%
  left_join(ETBR_MT_lmer, by = "Gene") %>%
  rename(MT_lmer_b = Estimate, MT_lmer_SE = SE, MT_lmer_P = P) %>%
  arrange(GeneType)

cat("Combined results:", nrow(dat_all), "genes\n")

# ============================================================
# 3. Prepare for path classification
# ============================================================

dat_etbr <- dat_all %>%
  mutate(Estimate = CN_lmer_b, SE = CN_lmer_SE, P = CN_lmer_P) %>%
  select(Gene, entrezgene_id, hgnc_symbol, gene_type, chromosome_name,
         start_position, end_position, strand, gene_description, GeneType,
         Estimate, SE, P)

# ============================================================
# 4. data.Scale.Lag2() — compute scaled fold-change at 4 dose-binned timepoints
# ============================================================
# This function is from plot_functions.R (not in repo), reconstructed from usage:
# The 4 timepoints correspond to dose bins: {0}, {1-2}, {3-4}, {5-6}
# (Dose levels 0,1,2,3,4,5,6 → binned into 4 groups)
# Returns mean and SD of log2FC at each bin, scaled relative to bin 3 (max effect)

data.Scale.Lag2 <- function(gene) {
  pheno_sub <- dat_pheno %>% subset(Chemical == "ETBR" & Curve == "Dosage")
  dose_numeric <- as.numeric(as.character(pheno_sub$Dose))

  # Bin doses: 0 → bin 1, 1-2 → bin 2, 3-4 → bin 3, 5-6 → bin 4
  bins <- cut(dose_numeric, breaks = c(-0.5, 0.5, 2.5, 4.5, 6.5),
              labels = 1:4)

  # Get expression values
  expr <- dat_rna[gene, pheno_sub$RNASeq_label]

  # Average across lanes per SampleName first
  df <- data.frame(
    SampleName = pheno_sub$SampleName,
    Y = as.numeric(expr),
    bin = bins
  ) %>%
    group_by(SampleName, bin) %>%
    summarise(Y = mean(Y), .groups = "drop")

  # Control mean (bin 1)
  ctrl_mean <- mean(df$Y[df$bin == "1"])

  # log2 fold change relative to control
  df$log2FC <- log2(df$Y / ctrl_mean)

  # Mean and SD per bin
  bin_stats <- df %>%
    group_by(bin) %>%
    summarise(Y = mean(log2FC), Ysd = sd(log2FC), .groups = "drop")

  return(data.frame(Y = bin_stats$Y, Ysd = bin_stats$Ysd))
}

# ============================================================
# 5. Compute scaled fold-change matrices
# ============================================================

cat("Computing scaled fold-change time series...\n")

GE_mean <- sapply(dat_etbr$Gene, function(x) {
  a <- tryCatch(data.Scale.Lag2(x), error = function(e) {
    data.frame(Y = rep(0, 4), Ysd = rep(1, 4))
  })
  if (a$Y[3] < 0) a$Y <- a$Y * -1
  Y_mult <- a$Y[3]
  if (Y_mult == 0) Y_mult <- 1  # prevent division by zero
  return(a$Y / Y_mult)
})

GE_sd <- sapply(dat_etbr$Gene, function(x) {
  a <- tryCatch(data.Scale.Lag2(x), error = function(e) {
    data.frame(Y = rep(0, 4), Ysd = rep(1, 4))
  })
  if (a$Y[3] < 0) a$Y <- a$Y * -1
  Y_mult <- a$Y[3]
  if (Y_mult == 0) Y_mult <- 1
  return(a$Ysd / Y_mult)
})

# ============================================================
# 6. Likelihood-based classification
# ============================================================

cat("Classifying response patterns...\n")

# Template vectors for the 3 response patterns
CN_mean  <- c(0, 0.5, 1, 0.5)   # Linear: rises then falls
MT_mean  <- c(0, 1, 1, 0)       # Switch: immediate response then recovery
lag_mean <- c(0, 0.25, 1, 1)    # Delayed: slow onset, sustained

# Compute likelihoods using timepoints 2 and 4
GE_L_CN  <- c()
GE_L_MT  <- c()
GE_L_lag <- c()

for (i in 1:ncol(GE_mean)) {
  GE_L_CN  <- c(GE_L_CN,
                 dnorm(CN_mean[2], GE_mean[2, i], GE_sd[2, i]) *
                   dnorm(CN_mean[4], GE_mean[4, i], GE_sd[4, i]))
  GE_L_MT  <- c(GE_L_MT,
                 dnorm(MT_mean[2], GE_mean[2, i], GE_sd[2, i]) *
                   dnorm(MT_mean[4], GE_mean[4, i], GE_sd[4, i]))
  GE_L_lag <- c(GE_L_lag,
                 dnorm(lag_mean[2], GE_mean[2, i], GE_sd[2, i]) *
                   dnorm(lag_mean[4], GE_mean[4, i], GE_sd[4, i]))
}

dat_etbr <- dat_etbr %>%
  mutate(L_CN = GE_L_CN, L_MT = GE_L_MT, L_lag = GE_L_lag)

# Classify: >2-fold ratio of top/second-best posterior → assigned, else "None"
dat_etbr$post_path <- apply(
  cbind(GE_L_CN / (GE_L_CN + GE_L_MT + GE_L_lag),
        GE_L_MT / (GE_L_CN + GE_L_MT + GE_L_lag),
        GE_L_lag / (GE_L_CN + GE_L_MT + GE_L_lag)),
  1, function(x) {
    x_i <- order(x, decreasing = TRUE)
    ifelse(x[x_i[1]] / x[x_i[2]] > 2, x_i[1], 4)
  })

dat_etbr <- dat_etbr %>%
  mutate(post_path = c("Linear", "Switch", "Delayed", "None")[post_path])

dat_etbr$chromosome_name <- factor(dat_etbr$chromosome_name,
                                    levels = c(1:22, "X", "MT"))

# ============================================================
# 7. Export
# ============================================================

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(dat_etbr, file.path(output_dir, "ETBR_paths.rds"))

cat("\n=== Classification Summary ===\n")
print(table(dat_etbr$post_path))
cat("\nBy GeneType:\n")
print(table(dat_etbr$GeneType, dat_etbr$post_path))
cat("\n=== DONE ===\n")
