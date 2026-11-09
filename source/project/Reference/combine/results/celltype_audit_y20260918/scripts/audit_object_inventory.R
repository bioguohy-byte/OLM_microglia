#!/usr/bin/env Rscript

# Read-only inventory for the original migrated project.  This script never
# changes source Seurat objects; it only records object structure/metadata.
suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: audit_object_inventory.R <manifest.tsv> <out_dir>")
}
manifest <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_rows <- list()
column_rows <- list()
value_rows <- list()

candidate_column <- function(x) {
  grepl("cell.?type|cluster|snn|ident|seurat|sample|condition|disease|group|batch", x,
        ignore.case = TRUE)
}

for (i in seq_len(nrow(manifest))) {
  study <- manifest$study[[i]]
  path <- manifest$object_path[[i]]
  message("Inspecting ", study, ": ", path)
  obj <- qread(path)
  if (!inherits(obj, "Seurat")) {
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      study = study, object_path = path, class = paste(class(obj), collapse = ";"),
      cells = NA_integer_, features = NA_integer_, default_assay = NA_character_,
      reductions = NA_character_, metadata_columns = NA_character_, stringsAsFactors = FALSE
    )
    rm(obj); gc()
    next
  }
  meta <- obj@meta.data
  assay_name <- DefaultAssay(obj)
  feature_n <- nrow(obj[[assay_name]])
  cols <- colnames(meta)
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    study = study, object_path = path, class = paste(class(obj), collapse = ";"),
    cells = ncol(obj), features = feature_n, default_assay = assay_name,
    reductions = paste(Reductions(obj), collapse = ";"),
    metadata_columns = paste(cols, collapse = ";"), stringsAsFactors = FALSE
  )
  for (col in cols) {
    column_rows[[length(column_rows) + 1L]] <- data.frame(
      study = study, column = col, class = paste(class(meta[[col]]), collapse = ";"),
      nonmissing = sum(!is.na(meta[[col]])), unique_values = length(unique(meta[[col]])),
      stringsAsFactors = FALSE
    )
    if (candidate_column(col) && length(unique(meta[[col]])) <= 60L) {
      tab <- sort(table(meta[[col]], useNA = "ifany"), decreasing = TRUE)
      value_rows[[length(value_rows) + 1L]] <- data.frame(
        study = study, column = col, value = names(tab), cells = as.integer(tab),
        stringsAsFactors = FALSE, check.names = FALSE
      )
    }
  }
  rm(meta, obj); gc()
}

write.csv(do.call(rbind, summary_rows), file.path(out_dir, "object_inventory.csv"), row.names = FALSE)
write.csv(do.call(rbind, column_rows), file.path(out_dir, "metadata_column_inventory.csv"), row.names = FALSE)
if (length(value_rows)) {
  write.csv(do.call(rbind, value_rows), file.path(out_dir, "candidate_label_value_counts.csv"), row.names = FALSE)
}
