#!/usr/bin/env Rscript
# generate_publication_figures.R
# From: chemcells_rnaseq_stat_X_plots.Rmd (Figures 1, 2, X, 4, 6, 7)

suppressPackageStartupMessages({
  library(readxl)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(ggh4x)
  library(RColorBrewer)
  library(tidyverse)
})

shared_dir <- Sys.getenv("SHARED_DIR", "../../shared")
source(file.path(shared_dir, "data_helpers.R"))
source(file.path(shared_dir, "regression_helpers.R"))

pheno_file      <- Sys.getenv("PHENO_FILE")
rna_file        <- Sys.getenv("RNA_FILE")
mtdna_file      <- Sys.getenv("MTDNA_FILE")
paths_file      <- Sys.getenv("PATHS_FILE")
annotation_file <- Sys.getenv("ANNOTATION_FILE")
mitocarta_file  <- Sys.getenv("MITO_CARTA_FILE")
go_results_file <- Sys.getenv("GO_RESULTS_FILE", "")
output_dir      <- Sys.getenv("OUTPUT_DIR", "/output")

cat("=== generate_publication_figures ===\n")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load data
# ============================================================

dat_pheno <- readRDS(pheno_file)
dat_rna   <- readRDS(rna_file)
dat_CN    <- readRDS(mtdna_file)
dat_etbr  <- readRDS(paths_file)
G_list    <- load_gene_list(annotation_file, filter_genes = rownames(dat_rna))
mc        <- load_mitocarta(mitocarta_file)

# Join deltaCT_Avg
dat_CN$PO <- factor(dat_CN$PO)
dat_pheno$deltaCT_Avg <- sapply(dat_pheno$SampleName, function(x) {
  vals <- dat_CN$deltaCT_Avg[which(dat_CN$Sample_Name == x & !is.na(dat_CN$deltaCT_Avg))]
  if (length(vals) == 0) return(NA)
  return(vals[1])
})

# ============================================================
# Helper: plot gene dose-response (replaces plot.Gene.Dose3)
# ============================================================

plot_gene_dose <- function(gene_symbol, chemical = "ETBR") {
  gene_id <- G_list$ensembl_gene_id[G_list$hgnc_symbol == gene_symbol]
  if (length(gene_id) == 0 || !gene_id[1] %in% rownames(dat_rna)) {
    return(ggplot() + ggtitle(paste(gene_symbol, "- not found")))
  }
  gene_id <- gene_id[1]

  pheno_sub <- dat_pheno %>% subset(Chemical == chemical & Curve == "Dosage")
  df <- data.frame(
    SampleName = pheno_sub$SampleName,
    Dose = as.numeric(as.character(pheno_sub$Dose)) * 25,
    Y = as.numeric(dat_rna[gene_id, pheno_sub$RNASeq_label])
  ) %>%
    group_by(SampleName, Dose) %>%
    summarise(Y = mean(Y), .groups = "drop")

  ggplot(df, aes(x = factor(Dose), y = Y)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 0.8, alpha = 0.6) +
    theme_classic() +
    xlab("EtBr (ng/mL)") + ylab("CPM") +
    ggtitle(gene_symbol)
}

# ============================================================
# Helper: plot gene longitudinal (replaces plot.Gene.Long2)
# ============================================================

plot_gene_longitudinal <- function(gene_symbol) {
  gene_id <- G_list$ensembl_gene_id[G_list$hgnc_symbol == gene_symbol]
  if (length(gene_id) == 0 || !gene_id[1] %in% rownames(dat_rna)) {
    return(ggplot() + ggtitle(paste(gene_symbol, "- not found")))
  }
  gene_id <- gene_id[1]

  pheno_sub <- dat_pheno %>% subset(Chemical == "ETBR" & Curve != "Dosage")
  pheno_sub$EtBr <- ifelse(pheno_sub$Dose == 6, "Experimental", "Control")
  pheno_sub$EtBr <- factor(pheno_sub$EtBr, levels = c("Control", "Experimental"))

  df <- data.frame(
    Time = pheno_sub$Time,
    EtBr = pheno_sub$EtBr,
    SampleName = pheno_sub$SampleName,
    Y = as.numeric(dat_rna[gene_id, pheno_sub$RNASeq_label])
  ) %>%
    group_by(SampleName, Time, EtBr) %>%
    summarise(Y = mean(Y), .groups = "drop")

  df_summary <- df %>%
    group_by(Time, EtBr) %>%
    summarise(Ymean = mean(Y), Ysd = sd(Y), .groups = "drop")

  ggplot(df_summary, aes(x = Time, y = Ymean, group = EtBr, color = EtBr)) +
    geom_line(alpha = 0.5) +
    geom_point(data = df, aes(x = Time, y = Y, color = EtBr),
               position = position_dodge(0.4), size = 0.8) +
    geom_errorbar(aes(ymin = Ymean - Ysd, ymax = Ymean + Ysd),
                  width = 1, position = position_dodge(0.4)) +
    theme_classic() +
    scale_color_brewer(type = "qual", palette = "Dark2") +
    guides(color = "none") +
    xlab("Time (hours)") + ylab("CPM") +
    ggtitle(gene_symbol)
}

# ============================================================
# Figure 1A: CN dose-response
# ============================================================

cat("Fig 1A: CN dose-response...\n")
pheno_dose <- dat_pheno %>% subset(Chemical == "ETBR" & Curve == "Dosage" & Lane == 1)
pheno_dose$Dose_ng <- as.numeric(as.character(pheno_dose$Dose)) * 25

png(file.path(output_dir, "fig1a_cn_dose.png"), width = 900, height = 600, res = 150)
ggplot(pheno_dose, aes(x = factor(Dose_ng), y = deltaCT_Avg)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1) +
  theme_classic() +
  xlab("EtBr (ng/mL)") + ylab("mtDNA-CN (deltaCT)")
dev.off()

# ============================================================
# Figure 1B: MT gene dose grid
# ============================================================

cat("Fig 1B: MT gene dose grid...\n")
mt_symbols <- dat_etbr %>%
  subset(chromosome_name == "MT") %>%
  arrange(start_position) %>%
  pull(hgnc_symbol)

mt_plots <- lapply(mt_symbols, function(g) {
  plot_gene_dose(g, "ETBR") +
    theme(rect = element_rect(fill = "transparent"),
          plot.title = element_text(size = 8))
})

png(file.path(output_dir, "fig1b_mt_gene_dose_grid.png"),
    width = 2400, height = 2400, res = 150)
print(ggarrange(plotlist = mt_plots, ncol = 4, nrow = ceiling(length(mt_plots) / 4)))
dev.off()

# ============================================================
# Figure 2A: MT position vs logFC
# ============================================================

cat("Fig 2A: MT position vs logFC...\n")

# Compute per-dose log2FC for each MT gene
cnfc <- sapply(1:6, function(i) {
  coef(lm(deltaCT_Avg ~ Dose, data = dat_pheno %>%
             subset(Curve == "Dosage" & Dose %in% c(0, i))))[2]
})

logFC1 <- data.frame(Dose = character(), Gene = character(),
                      logFC = numeric(), logSD = numeric(),
                      stringsAsFactors = FALSE)

for (gene_sym in mt_symbols) {
  gene_id <- G_list$ensembl_gene_id[G_list$hgnc_symbol == gene_sym]
  if (length(gene_id) == 0 || !gene_id[1] %in% rownames(dat_rna)) next
  gene_id <- gene_id[1]

  pheno_sub <- dat_pheno %>% subset(Chemical == "ETBR" & Curve == "Dosage")
  df <- data.frame(
    Dose = as.numeric(as.character(pheno_sub$Dose)),
    SampleName = pheno_sub$SampleName,
    Y = as.numeric(dat_rna[gene_id, pheno_sub$RNASeq_label])
  )

  for (d in unique(df$Dose[df$Dose != 0])) {
    tryCatch({
      a <- summary(lmer(log2(Y) ~ I(Dose != 0) + (1 | SampleName),
                         df %>% subset(Dose %in% c(0, d))))$coef
      logFC1 <- rbind(logFC1, data.frame(
        Dose = d, Gene = gene_sym, logFC = a[2, 1], logSD = a[2, 2]
      ))
    }, error = function(e) NULL)
  }
}

if (nrow(logFC1) > 0) {
  logFC1$Gene <- factor(logFC1$Gene, levels = mt_symbols)
  logFC1$Dose_label <- factor(as.numeric(logFC1$Dose) * 25)
  logFC1$X <- sapply(logFC1$Gene, function(x) {
    info <- dat_etbr %>% subset(hgnc_symbol == x)
    if (nrow(info) > 0) (info$start_position[1] + info$end_position[1]) / 2 else NA
  })

  png(file.path(output_dir, "fig2a_mt_position_logfc.png"),
      width = 1200, height = 800, res = 150)
  print(ggplot(logFC1, aes(x = X, y = logFC, group = Dose_label, color = Dose_label)) +
          geom_line() + geom_point() +
          theme_classic() + xlab("Position (bp)") + ylab("log2 FC"))
  dev.off()
}

# ============================================================
# Figure X: Volcano plot
# ============================================================

cat("Fig X: Volcano plot...\n")
dat_plot <- dat_etbr %>%
  mutate(log2FC = Estimate / log(2), log10P = log10(P)) %>%
  subset(chromosome_name != "MT")

p_cut <- max(dat_plot$P[which(p.adjust(dat_plot$P, "holm") < 0.05)], na.rm = TRUE)

# Label MitoCarta central dogma genes
dogma_symbols <- mc$genes %>%
  subset(grepl("dogma", MitoCarta3.0_MitoPathways)) %>%
  pull(Symbol)

dat_plot$delabel <- ifelse(
  dat_plot$log10P <= log10(p_cut) & dat_plot$hgnc_symbol %in% dogma_symbols,
  dat_plot$hgnc_symbol, ""
)

png(file.path(output_dir, "figX_volcano.png"), width = 1200, height = 800, res = 150)
print(ggplot(dat_plot, aes(x = log2FC, y = -log10P, label = delabel)) +
        geom_point(size = 0.5, alpha = 0.5) +
        geom_label_repel(max.overlaps = 100, label.size = NA,
                          fill = alpha("white", 0.1), size = 2) +
        theme_classic() +
        geom_hline(yintercept = -log10(p_cut), col = "red"))
dev.off()

# ============================================================
# Figure 4: Glycolysis genes
# ============================================================

cat("Fig 4: Glycolysis genes...\n")
glycolysis_genes <- c("SLC2A1", "HK1", "GPI", "PFKL", "ALDOA", "TPI1",
                       "GAPDH", "PGK1", "PGM1", "ENO1", "PKM")

glyc_plots <- lapply(glycolysis_genes, function(g) {
  plot_gene_dose(g, "ETBR") +
    theme(plot.title = element_text(size = 8))
})

png(file.path(output_dir, "fig4_glycolysis_genes.png"),
    width = 2000, height = 1200, res = 150)
print(ggarrange(plotlist = glyc_plots, ncol = 4, nrow = 3))
dev.off()

# ============================================================
# Figure 6A: Longitudinal CN time course
# ============================================================

cat("Fig 6A: CN longitudinal...\n")
dat_lon_cn <- dat_pheno %>%
  subset(Chemical == "ETBR" & Curve != "Dosage" & Lane == 1)
dat_lon_cn$EtBr <- ifelse(dat_lon_cn$Dose == 6, "Experimental", "Control")
dat_lon_cn$EtBr <- factor(dat_lon_cn$EtBr, levels = c("Control", "Experimental"))

dat_lon_summary <- dat_lon_cn %>%
  group_by(Time, EtBr) %>%
  summarise(Y = mean(deltaCT_Avg), Y_sd = sd(deltaCT_Avg), .groups = "drop")

png(file.path(output_dir, "fig6a_cn_longitudinal.png"), width = 900, height = 600, res = 150)
print(ggplot(dat_lon_summary, aes(x = Time, y = Y, group = EtBr)) +
        geom_line(aes(color = EtBr), alpha = 0.5) +
        geom_point(data = dat_lon_cn, aes(y = deltaCT_Avg, color = EtBr),
                   position = position_dodge(0.4), size = 1) +
        geom_errorbar(aes(ymin = Y - Y_sd, ymax = Y + Y_sd, color = EtBr),
                      width = 1, position = position_dodge(0.4)) +
        theme_classic() +
        scale_color_brewer(type = "qual", palette = "Dark2") +
        ylab("mtDNA-CN (deltaCT)") + xlab("Time (hours)"))
dev.off()

# ============================================================
# Figure 6B: MT gene longitudinal grid
# ============================================================

cat("Fig 6B: MT gene longitudinal grid...\n")
mt_lon_plots <- lapply(mt_symbols, function(g) {
  plot_gene_longitudinal(g) +
    theme(legend.position = "none", plot.title = element_text(size = 8))
})

png(file.path(output_dir, "fig6b_mt_gene_longitudinal_grid.png"),
    width = 2400, height = 2400, res = 150)
print(ggarrange(plotlist = mt_lon_plots, ncol = 4,
                nrow = ceiling(length(mt_lon_plots) / 4)))
dev.off()

# ============================================================
# Figure 7: Scaled lag comparison (PSAT1)
# ============================================================

cat("Fig 7: Scaled lag comparison...\n")

# data.Scale.Lag2 implementation (same as in classify_response_patterns)
data_scale_lag2 <- function(gene_symbol) {
  gene_id <- G_list$ensembl_gene_id[G_list$hgnc_symbol == gene_symbol]
  if (length(gene_id) == 0 || !gene_id[1] %in% rownames(dat_rna)) return(NULL)
  gene_id <- gene_id[1]

  pheno_sub <- dat_pheno %>% subset(Chemical == "ETBR" & Curve == "Dosage")
  dose_numeric <- as.numeric(as.character(pheno_sub$Dose))
  bins <- cut(dose_numeric, breaks = c(-0.5, 0.5, 2.5, 4.5, 6.5), labels = 1:4)

  df <- data.frame(
    SampleName = pheno_sub$SampleName, bin = bins,
    Y = as.numeric(dat_rna[gene_id, pheno_sub$RNASeq_label])
  ) %>%
    group_by(SampleName, bin) %>%
    summarise(Y = mean(Y), .groups = "drop")

  ctrl_mean <- mean(df$Y[df$bin == "1"])
  df$log2FC <- log2(df$Y / ctrl_mean)

  bin_stats <- df %>%
    group_by(bin) %>%
    summarise(Y = mean(log2FC), Ysd = sd(log2FC), .groups = "drop")

  return(bin_stats)
}

# Compute scaled reference curves
cn_template  <- data.frame(bin = 1:4, Y = c(0, 0.5, 1, 0.5), curve = "Linear: mtDNA-CN")
mt_template  <- data.frame(bin = 1:4, Y = c(0, 1, 1, 0), curve = "Switch: MT Med FC")

psat1_data <- data_scale_lag2("PSAT1")
if (!is.null(psat1_data)) {
  # Scale relative to bin 3
  if (psat1_data$Y[3] < 0) psat1_data$Y <- psat1_data$Y * -1
  scale_factor <- psat1_data$Y[3]
  if (scale_factor == 0) scale_factor <- 1
  psat1_scaled <- data.frame(bin = 1:4, Y = psat1_data$Y / scale_factor,
                              curve = "Delayed: PSAT1")

  plot_data <- rbind(cn_template, mt_template, psat1_scaled)

  png(file.path(output_dir, "fig7_scaled_lag.png"), width = 1200, height = 800, res = 150)
  print(ggplot(plot_data, aes(x = bin, y = Y, color = curve)) +
          geom_line(size = 1) + geom_point(size = 2) +
          theme_classic() +
          scale_color_manual(values = c("black", "red", "blue")) +
          ylim(-0.4, 1.2) + ylab("Scaled FC") + xlab("Dose bin") +
          theme(legend.position = "bottom"))
  dev.off()
}

# ============================================================
# Summary
# ============================================================

figures <- list.files(output_dir, pattern = "\\.png$")
cat("\n=== Generated", length(figures), "figures ===\n")
for (f in figures) {
  size <- file.size(file.path(output_dir, f))
  cat("  ", f, ":", round(size / 1024), "KB\n")
}
cat("=== DONE ===\n")
