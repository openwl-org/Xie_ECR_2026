#!/usr/bin/env Rscript
# run_go_enrichment.R — Stratified GO enrichment
# From: chemcells_rnaseq_stat_X_plots.Rmd (Figure 3, Table 3, Figure S2, Sup Figure 1)

suppressPackageStartupMessages({
  library(limma)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  library(rrvgo)
  library(tidyverse)
})

paths_file    <- Sys.getenv("PATHS_FILE")
cn_anova_file <- Sys.getenv("CN_ANOVA_FILE")
mt_lmer_file  <- Sys.getenv("MT_LMER_FILE")
output_dir    <- Sys.getenv("OUTPUT_DIR", "/output")

cat("=== run_go_enrichment ===\n")

# ============================================================
# 1. Load and prepare data
# ============================================================

dat_etbr <- readRDS(paths_file) %>% subset(gene_type != "Mt_tRNA")

# Load additional regression results for significance columns
ETBR_CN_anova <- readRDS(cn_anova_file)
ETBR_MT_lmer  <- readRDS(mt_lmer_file)

# Join additional columns
dat_etbr <- dat_etbr %>%
  left_join(ETBR_MT_lmer, by = "Gene", suffix = c("", "_mt")) %>%
  left_join(ETBR_CN_anova, by = "Gene") %>%
  mutate(
    Sig_lin = p.adjust(P, "holm") < 0.05,
    Sig_mt  = p.adjust(P_mt, "holm") < 0.05,
    Sig_sq  = p.adjust(P12, "holm") < 0.05,
    Sig_sp2 = p.adjust(P13, "holm") < 0.05,
    Sig_sp3 = p.adjust(P14, "holm") < 0.05
  )

dat_etbr$post_path <- factor(dat_etbr$post_path,
                               levels = c("Switch", "Linear", "Delayed", "None"))

# Universe: all non-MT genes with entrezgene_id
uni <- unique(dat_etbr$entrezgene_id[!is.na(dat_etbr$entrezgene_id)])
uni <- uni[!uni %in% dat_etbr$entrezgene_id[dat_etbr$chromosome_name == "MT"]]
uni <- as.character(uni)

# Significant genes
top_genes <- dat_etbr %>%
  subset(chromosome_name != "MT" & !is.na(entrezgene_id)) %>%
  subset(Sig_lin | Sig_mt | Sig_sq | Sig_sp2 | Sig_sp3)

cat("Universe:", length(uni), "genes | Significant:", nrow(top_genes), "genes\n")

# ============================================================
# 2. enrichGO: stratified by direction
# ============================================================

run_enrichGO <- function(gene_ids, label) {
  cat("  enrichGO:", label, "(", length(gene_ids), "genes )\n")
  tryCatch({
    enrichGO(gene = gene_ids, universe = uni, OrgDb = org.Hs.eg.db,
             ont = "BP", pAdjustMethod = "BH",
             pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  }, error = function(e) {
    cat("    Error:", e$message, "\n")
    NULL
  })
}

ego_results <- list()

# All significant
ids_all <- top_genes$entrezgene_id[!is.na(top_genes$entrezgene_id)]
ego_results$all <- run_enrichGO(ids_all, "all")

# Up (Estimate < 0 means increased expression with CN depletion)
ids_up <- top_genes$entrezgene_id[!is.na(top_genes$entrezgene_id) &
                                     top_genes$Estimate < 0 & top_genes$Sig_lin]
ego_results$all_up <- run_enrichGO(ids_up, "all_up")

# Down
ids_down <- top_genes$entrezgene_id[!is.na(top_genes$entrezgene_id) &
                                       top_genes$Estimate > 0 & top_genes$Sig_lin]
ego_results$all_down <- run_enrichGO(ids_down, "all_down")

# By pattern + direction
for (pattern in c("Delayed", "Linear", "Switch", "None")) {
  for (direction in c("up", "down")) {
    if (direction == "up") {
      ids <- top_genes$entrezgene_id[!is.na(top_genes$entrezgene_id) &
                                        top_genes$Sig_lin &
                                        top_genes$post_path == pattern &
                                        top_genes$Estimate < 0]
    } else {
      ids <- top_genes$entrezgene_id[!is.na(top_genes$entrezgene_id) &
                                        top_genes$Sig_lin &
                                        top_genes$post_path == pattern &
                                        top_genes$Estimate > 0]
    }
    key <- paste0(tolower(pattern), "_", direction)
    ego_results[[key]] <- run_enrichGO(ids, key)
  }
}

# ============================================================
# 3. goana: 15-way comparison table
# ============================================================

cat("\nRunning goana (15-way)...\n")

top_linear <- top_genes %>% subset(Sig_lin)

go_etbr <- tryCatch({
  goana(
    de = list(
      all     = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id)],
      all_u   = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id) & top_linear$Estimate < 0],
      all_d   = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id) & top_linear$Estimate > 0],
      mt      = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id) & top_linear$post_path == "Switch"],
      mt_u    = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id) & top_linear$post_path == "Switch" & top_linear$Estimate < 0],
      mt_d    = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id) & top_linear$post_path == "Switch" & top_linear$Estimate > 0],
      cn      = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id) & top_linear$post_path == "Linear"],
      cn_u    = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id) & top_linear$post_path == "Linear" & top_linear$Estimate < 0],
      cn_d    = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id) & top_linear$post_path == "Linear" & top_linear$Estimate > 0],
      delay   = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id) & top_linear$post_path == "Delayed"],
      delay_u = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id) & top_linear$post_path == "Delayed" & top_linear$Estimate < 0],
      delay_d = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id) & top_linear$post_path == "Delayed" & top_linear$Estimate > 0],
      none    = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id) & top_linear$post_path == "None"],
      none_u  = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id) & top_linear$post_path == "None" & top_linear$Estimate < 0],
      none_d  = top_linear$entrezgene_id[!is.na(top_linear$entrezgene_id) & top_linear$post_path == "None" & top_linear$Estimate > 0]
    ),
    universe = uni, species = "Hs", FDR = 0.05
  )
}, error = function(e) {
  cat("goana error:", e$message, "\n")
  NULL
})

# ============================================================
# 4. rrvgo treemap simplification
# ============================================================

cat("\nComputing rrvgo treemap simplifications...\n")

treemap_results <- list()
for (name in c("all_up", "all_down", "delayed_up", "delayed_down")) {
  ego <- ego_results[[name]]
  if (is.null(ego)) next
  go_result <- ego@result %>% subset(p.adjust < 0.05)
  if (nrow(go_result) == 0) next
  tryCatch({
    sim <- calculateSimMatrix(go_result$ID, orgdb = "org.Hs.eg.db",
                               ont = "BP", method = "Rel")
    scores <- setNames(-log10(go_result$p.adjust), go_result$ID)
    reduced <- reduceSimMatrix(sim, scores, threshold = 0.7,
                                orgdb = "org.Hs.eg.db")
    treemap_results[[name]] <- list(simMatrix = sim, reduced = reduced)
    cat("  ", name, ":", nrow(reduced), "clusters\n")
  }, error = function(e) {
    cat("  ", name, ": rrvgo error -", e$message, "\n")
  })
}

# ============================================================
# 5. Export
# ============================================================

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(ego_results, file.path(output_dir, "GO_enrichGO_results.rds"))
saveRDS(go_etbr, file.path(output_dir, "GO_goana_table.rds"))
saveRDS(treemap_results, file.path(output_dir, "GO_treemap_results.rds"))

if (!is.null(go_etbr)) {
  write.table(go_etbr %>% tibble::rownames_to_column("GOID"),
              file.path(output_dir, "GO_goana_table.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("\n=== enrichGO results summary ===\n")
for (name in names(ego_results)) {
  ego <- ego_results[[name]]
  n <- if (!is.null(ego)) nrow(ego@result %>% subset(p.adjust < 0.05)) else 0
  cat("  ", name, ":", n, "significant terms\n")
}

cat("\n=== DONE ===\n")
