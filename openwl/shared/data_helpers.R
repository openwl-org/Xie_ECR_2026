# data_helpers.R
# Shared data loading and gene annotation helpers from Xie et al. (2026)
# Extracted from chemcells_rnaseq_data_1_read_data.Rmd and stat_2_lon_track.Rmd

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
})

#' Load gene annotation list from mart_export.txt
#' @param path Path to mart_export.txt
#' @param filter_genes Optional character vector of gene IDs to subset
load_gene_list <- function(path, filter_genes = NULL) {
  G_list <- read.delim(path)
  names(G_list) <- c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol",
                     "chromosome_name", "start_position", "end_position",
                     "strand", "gene_type", "gene_description", "gene_name")
  G_list <- G_list[!duplicated(G_list$ensembl_gene_id), ]
  if (!is.null(filter_genes)) {
    G_list <- G_list %>% subset(ensembl_gene_id %in% filter_genes)
  }
  G_list$chromosome_name <- factor(G_list$chromosome_name,
                                    levels = c(1:22, "X", "MT"))
  return(G_list)
}

#' Load MitoCarta 3.0 gene and pathway sheets
#' @param path Path to Human.MitoCarta3.0.xlsx
load_mitocarta <- function(path) {
  genes <- read_xlsx(path, sheet = "A Human MitoCarta3.0",
                     col_types = "text") %>% as.data.frame()
  paths <- read_xlsx(path, sheet = "C MitoPathways",
                     col_types = "text") %>% as.data.frame()
  return(list(genes = genes, paths = paths))
}

#' Classify genes by GeneType using chromosome and MitoCarta
#' @param genes Character vector of Ensembl gene IDs
#' @param G_list Gene annotation data.frame
#' @param mitocarta_genes MitoCarta genes data.frame
classify_gene_type <- function(genes, G_list, mitocarta_genes) {
  gene_info <- data.frame(Gene = genes) %>%
    left_join(G_list, by = c("Gene" = "ensembl_gene_id"))

  gene_info$GeneType <- case_when(
    gene_info$chromosome_name == "MT" ~ "3 MT encoded",
    gene_info$Gene %in% mitocarta_genes$EnsemblGeneID_mapping_version_20200130 ~ "2 MitoCarta",
    TRUE ~ "1 Nuclear"
  )
  return(gene_info)
}

#' PCA outlier detection on first N principal components
#' @param pca_result prcomp object
#' @param sd_num Number of standard deviations for outlier threshold
#' @param n_pcs Number of PCs to check (default 10)
find_PCA_outliers <- function(pca_result, sd_num, n_pcs = 10) {
  all_outliers <- c()
  for (i in 1:min(n_pcs, ncol(pca_result$x))) {
    scores <- pca_result$x[, i]
    upper <- rownames(pca_result$x)[scores > (mean(scores) + sd_num * sd(scores))]
    lower <- rownames(pca_result$x)[scores < (mean(scores) - sd_num * sd(scores))]
    all_outliers <- c(all_outliers, upper, lower)
  }
  return(unique(all_outliers))
}

#' Compute Med_log2FC_MT — median log2 fold change of MT protein-coding genes
#' vs condition-matched control means.
#' @param pheno Phenotype data.frame with Chemical, Curve, Time, Dose, RNASeq_label
#' @param norm_tmm TMM-normalized CPM matrix (genes x samples)
#' @param G_list Gene annotation data.frame
compute_med_log2fc_mt <- function(pheno, norm_tmm, G_list) {
  mt_genes <- G_list$ensembl_gene_id[G_list$chromosome_name == "MT" &
                                       G_list$gene_type == "protein_coding"]
  mt_genes <- mt_genes[mt_genes %in% rownames(norm_tmm)]

  # Compute condition-matched NC means for each MT gene
  conditions <- pheno %>% select(Chemical, Curve, Time) %>% distinct()
  MT_NC <- conditions

  for (gene in mt_genes) {
    NC_means <- c()
    for (j in 1:nrow(conditions)) {
      sample_j <- pheno %>%
        subset(Chemical == conditions$Chemical[j] &
                 Curve == conditions$Curve[j] &
                 Time == conditions$Time[j] &
                 Dose == 0) %>%
        select(RNASeq_label) %>% unlist()
      NC_means <- c(NC_means, mean(unlist(norm_tmm[gene, sample_j])))
    }
    MT_NC[[gene]] <- NC_means
  }

  # Compute per-sample median log2 FC
  Median_log2FC_MT <- c()
  for (i in 1:nrow(pheno)) {
    Gene_log2FC_MT <- c()
    RNA_sample <- pheno$RNASeq_label[i]
    NC_i <- which(pheno$Chemical[i] == MT_NC$Chemical &
                    pheno$Curve[i] == MT_NC$Curve &
                    pheno$Time[i] == MT_NC$Time)
    for (gene in mt_genes) {
      Gene_log2FC_MT <- c(Gene_log2FC_MT,
                          log2(unlist(norm_tmm[gene, RNA_sample]) /
                                 unlist(MT_NC[NC_i, gene])))
    }
    Median_log2FC_MT <- c(Median_log2FC_MT, median(Gene_log2FC_MT))
  }

  return(Median_log2FC_MT)
}
