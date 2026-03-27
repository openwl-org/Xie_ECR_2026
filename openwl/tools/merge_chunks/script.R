#!/usr/bin/env Rscript
# merge_chunks.R — Merge chunked regression RDS outputs via rbind
# Finds *_chunk*.rds files, groups by base name, merges into single RDS

input_dirs <- Sys.getenv("INPUT_DIRS", "")
output_dir <- Sys.getenv("OUTPUT_DIR", "/output")

cat("=== merge_chunks ===\n")

# Parse comma-separated input directories
dirs <- trimws(unlist(strsplit(input_dirs, ",")))
dirs <- dirs[dirs != ""]

if (length(dirs) == 0) {
  stop("INPUT_DIRS is empty — no directories to scan for chunk files")
}

# Find all chunk RDS files across all input directories
chunk_files <- character(0)
for (d in dirs) {
  if (!dir.exists(d)) {
    cat("Warning: directory does not exist:", d, "\n")
    next
  }
  found <- list.files(d, pattern = "_chunk[0-9]+\\.rds$", full.names = TRUE)
  chunk_files <- c(chunk_files, found)
  cat("Found", length(found), "chunk files in", d, "\n")
}

if (length(chunk_files) == 0) {
  stop("No *_chunk*.rds files found in any input directory")
}

# Group by base name (strip _chunkNN suffix)
base_names <- sub("_chunk[0-9]+\\.rds$", ".rds", basename(chunk_files))
groups <- split(chunk_files, base_names)

cat("Merging", length(groups), "output groups from", length(chunk_files), "chunk files\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (out_name in names(groups)) {
  files <- sort(groups[[out_name]])  # sort ensures chunk order
  cat("  ", out_name, ":", length(files), "chunks...")

  chunks <- lapply(files, readRDS)
  merged <- do.call(rbind, chunks)

  out_path <- file.path(output_dir, out_name)
  saveRDS(merged, out_path)
  cat(" merged", nrow(merged), "rows\n")
}

cat("=== DONE ===\n")
