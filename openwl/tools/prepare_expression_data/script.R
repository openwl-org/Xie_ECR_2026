#!/usr/bin/env Rscript
# prepare_expression_data.R
# From: chemcells_rnaseq_data_1_read_data.Rmd
# Xie et al. (2026) — EtBr mtDNA-CN depletion study

suppressPackageStartupMessages({
  library(readxl)
  library(edgeR)
  library(ggplot2)
  library(ggfortify)
  library(RColorBrewer)
  library(tidyverse)
})

# Source shared helpers
shared_dir <- Sys.getenv("SHARED_DIR", "../../shared")
source(file.path(shared_dir, "data_helpers.R"))

# Input paths from environment
count_file      <- Sys.getenv("COUNT_FILE")
bridge_file     <- Sys.getenv("BRIDGE_FILE")
deltact_file    <- Sys.getenv("DELTA_CT_FILE")
annotation_file <- Sys.getenv("ANNOTATION_FILE")
output_dir      <- Sys.getenv("OUTPUT_DIR", "/output")

cat("=== prepare_expression_data ===\n")
cat("Count file:", count_file, "\n")
cat("Bridge file:", bridge_file, "\n")
cat("DeltaCT file:", deltact_file, "\n")
cat("Annotation file:", annotation_file, "\n")

# ============================================================
# 1. Read raw data
# ============================================================

mrna.count <- read.csv(count_file)
chem.bridge <- read.csv(bridge_file)[1:112, 1:18]
mtdna.data <- as.data.frame(read_xlsx(deltact_file)[, -1])

# ============================================================
# 2. Parse sample metadata from column names
# ============================================================

sample_no <- sapply(names(mrna.count)[-1], function(x) {
  as.numeric(substring(strsplit(x, "_")[[1]][3], 7))
})
lane_no <- sapply(names(mrna.count)[-1], function(x) {
  as.numeric(substring(strsplit(x, "_")[[1]][5], 4))
})

mrna.chara <- data.frame(
  RandomNumber = sapply(sample_no, function(x) {
    chem.bridge$SampleNumber[which(chem.bridge$RandomNumber == x)]
  }),
  SampleNumber = sample_no,
  Lane = lane_no,
  SampleName = sapply(sample_no, function(x) {
    chem.bridge$SampleName[which(chem.bridge$RandomNumber == x)]
  }),
  RNASeq_label = names(mrna.count)[-1],
  Chemical = NA,
  Dose = NA,
  Experiment = "Dosage",
  Curve = NA,
  Time = 48,
  DMSO = FALSE,
  row.names = NULL
)

# ============================================================
# 3. Annotate experiment type, chemical, dose, time, curve
# ============================================================

mrna.chara$DMSO <- grepl("DMSO", mrna.chara$SampleName)

# Longitudinal samples
mrna.chara$Experiment[grepl("BaseLine", mrna.chara$SampleName) |
                        grepl("hr", mrna.chara$SampleName)] <- "Longitudinal"
mrna.chara$Chemical[mrna.chara$Experiment == "Longitudinal"] <- "ETBR"

mrna.chara$Curve[mrna.chara$Experiment == "Longitudinal"] <- sapply(
  mrna.chara$SampleName[mrna.chara$Experiment == "Longitudinal"],
  function(x) strsplit(strsplit(x, "_")[[1]][1], "-")[[1]][3]
)
mrna.chara$Curve[mrna.chara$Experiment == "Longitudinal" &
                   grepl("BaseLine", mrna.chara$SampleName)] <- "Treatment"
mrna.chara$Curve[mrna.chara$Experiment == "Dosage"] <- "Dosage"

# Time parsing
mrna.chara$Time[mrna.chara$Experiment == "Longitudinal" &
                  grepl("BaseLine", mrna.chara$SampleName)] <- 0
mrna.chara$Time[mrna.chara$Experiment == "Longitudinal" &
                  !grepl("BaseLine", mrna.chara$SampleName)] <-
  sapply(mrna.chara$SampleName[mrna.chara$Experiment == "Longitudinal" &
                                  !grepl("BaseLine", mrna.chara$SampleName)],
         function(x) as.numeric(gsub("([0-9]+).*$", "\\1", substring(x, 1, 3))))

mrna.chara$Time[mrna.chara$Curve == "Recovery"] <-
  mrna.chara$Time[mrna.chara$Curve == "Recovery"] + 48

# Dose parsing
mrna.chara$Dose[grepl("NC", mrna.chara$SampleName) |
                  grepl("DMSO", mrna.chara$SampleName)] <- 0
mrna.chara$Dose[mrna.chara$Experiment == "Longitudinal" &
                  is.na(mrna.chara$Dose) &
                  mrna.chara$Chemical == "Acet"] <- 7
mrna.chara$Dose[mrna.chara$Experiment == "Longitudinal" &
                  is.na(mrna.chara$Dose)] <- 6

# Chemical and dose from SampleName (parse_number / unit)
for (i in 1:3) {
  chem <- c("Acet", "ETBR", "RESV")[i]
  unit <- c(0.25, 25, 4)[i]
  mrna.chara$Chemical[grep(chem, mrna.chara$SampleName)] <- chem
  n <- which(grepl(chem, mrna.chara$SampleName) & is.na(mrna.chara$Dose))
  mrna.chara$Dose[n] <- sapply(mrna.chara$SampleName[n], function(x) {
    parse_number(strsplit(x, "-")[[1]][1]) / unit
  })
}

# ============================================================
# 4. Filter to ETBR only
# ============================================================

mrna.chara <- mrna.chara %>% subset(Chemical == "ETBR")
mrna.count <- mrna.count[, c("Gene", mrna.chara$RNASeq_label)]

# Parse replicate number
mrna.chara$Rep <- sapply(mrna.chara$SampleName, function(x) {
  substr(x, nchar(x), nchar(x))
})
mrna.chara$Rep[which(mrna.chara$Rep == "o")] <- 1

pheno.out <- mrna.chara[, c("SampleName", "Rep", "Lane", "RNASeq_label",
                             "Chemical", "Dose", "Curve", "Time")]

# Factor levels
pheno.out$Rep <- factor(pheno.out$Rep)
pheno.out$Lane <- factor(pheno.out$Lane)
pheno.out$Dose <- factor(pheno.out$Dose)
pheno.out$Time <- factor(pheno.out$Time, levels = seq(0, 192, 24))
pheno.out$Curve <- factor(pheno.out$Curve,
                           levels = c("Dosage", "Treatment", "Recovery"))
pheno.out <- pheno.out %>% arrange(Curve, Time, Dose, Rep, Lane)

# ============================================================
# 5. TMM normalization
# ============================================================

dat.cnt <- as.matrix(mrna.count[, -1])
row.names(dat.cnt) <- mrna.count$Gene

# Gene filtering: median count > 49 in controls, ENSG prefix
gene_median <- apply(dat.cnt[, mrna.chara$RNASeq_label[mrna.chara$Dose == 0]],
                     1, median)
valid_gene <- which(gene_median > 49 & grepl("ENSG", mrna.count$Gene))

cat("Valid genes (pre-outlier):", length(valid_gene), "\n")

# TMM + CPM
y <- DGEList(dat.cnt[valid_gene, ])
y <- calcNormFactors(y, method = "TMM")
norm.tmm <- cpm(y)
row.names(norm.tmm) <- mrna.count$Gene[valid_gene]

# ============================================================
# 6. PCA outlier detection
# ============================================================

log.pca <- prcomp(t(log(norm.tmm + 1)), scale. = TRUE)
outliers <- find_PCA_outliers(log.pca, sd_num = 4)
cat("PCA outliers at 4 SD:", length(outliers), "\n")

mrna.chara$Exclude <- FALSE
mrna.chara$Exclude[mrna.chara$RNASeq_label %in% outliers] <- TRUE
mrna.chara$Exclude[grep("RESV_NC", mrna.chara$SampleName)] <- TRUE
mrna.chara$Exclude[which("ETBR_NC1" == mrna.chara$SampleName)] <- TRUE

cat("Samples excluded:", sum(mrna.chara$Exclude), "\n")

# Rerun without outliers
mrna.chara <- mrna.chara[!mrna.chara$Exclude, ]
pheno.out <- pheno.out %>% subset(RNASeq_label %in% mrna.chara$RNASeq_label)

dat.cnt <- as.matrix(mrna.count[, mrna.chara$RNASeq_label])
gene_median <- apply(dat.cnt[, mrna.chara$RNASeq_label[mrna.chara$Dose == 0 &
                                                          !mrna.chara$Exclude]],
                     1, median)
valid_gene <- which(gene_median > 49 & grepl("ENSG", mrna.count$Gene))

cat("Valid genes (post-outlier):", length(valid_gene), "\n")

y <- DGEList(dat.cnt[valid_gene, ])
y <- calcNormFactors(y, method = "TMM")
norm.tmm <- cpm(y)
row.names(norm.tmm) <- mrna.count$Gene[valid_gene]

# ============================================================
# 7. Compute Med_log2FC_MT
# ============================================================

G_list <- load_gene_list(annotation_file, filter_genes = rownames(norm.tmm))
pheno.out$Med_log2FC_MT <- compute_med_log2fc_mt(pheno.out, norm.tmm, G_list)

cat("Med_log2FC_MT range:",
    round(range(pheno.out$Med_log2FC_MT, na.rm = TRUE), 3), "\n")

# ============================================================
# 8. Export
# ============================================================

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(pheno.out, file = file.path(output_dir, "Experiment_Samples.rds"))
saveRDS(norm.tmm, file = file.path(output_dir, "RNA_TMM.rds"))
saveRDS(mtdna.data, file = file.path(output_dir, "mtDNA.rds"))

cat("\n=== Output files ===\n")
cat("  Experiment_Samples.rds:", nrow(pheno.out), "samples\n")
cat("  RNA_TMM.rds:", nrow(norm.tmm), "genes x", ncol(norm.tmm), "samples\n")
cat("  mtDNA.rds:", nrow(mtdna.data), "rows\n")

dose_levels <- sort(unique(as.numeric(as.character(pheno.out$Dose))))
cat("  Dose levels:", paste(dose_levels, collapse = ", "), "\n")
cat("  Curves:", paste(levels(pheno.out$Curve), collapse = ", "), "\n")
cat("\n=== DONE ===\n")
