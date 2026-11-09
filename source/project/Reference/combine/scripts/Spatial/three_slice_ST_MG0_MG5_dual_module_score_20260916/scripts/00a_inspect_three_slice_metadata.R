#!/usr/bin/env Rscript

# Inspect only the metadata schemas of the three fixed spatial microglia inputs.
# This diagnostic does not modify, subset or save any Seurat object.

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/three_slice_ST_MG0_MG5_dual_module_score_20260916")
result_dir <- file.path(project_dir, "results", "00a_metadata_inspection")
log_dir <- file.path(project_dir, "logs")
tmp_dir <- file.path(project_dir, "tmp")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)
.libPaths(c(file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library"), file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_libs"), .libPaths()))

suppressPackageStartupMessages(library(Seurat))

slice_paths <- c(
  s2r1 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r1_microglia_ST_MG0_MG5_from_v2.rds"),
  s2r2 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r2_microglia_ST_MG0_MG5_from_v2.rds"),
  s2r3 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r3_microglia_ST_MG0_MG5_from_v2.rds")
)

metadata_summary <- data.frame()
metadata_values <- data.frame()
for (slice_id in names(slice_paths)) {
  slice_object <- readRDS(slice_paths[[slice_id]])
  metadata <- slice_object[[]]
  for (metadata_column in colnames(metadata)) {
    column_values <- metadata[[metadata_column]]
    non_missing_values <- column_values[!is.na(column_values)]
    distinct_values <- unique(as.character(non_missing_values))
    metadata_summary <- rbind(
      metadata_summary,
      data.frame(
        slice = slice_id,
        metadata_column = metadata_column,
        class = paste(class(column_values), collapse = ";"),
        missing_cells = sum(is.na(column_values)),
        unique_non_missing_values = length(distinct_values),
        stringsAsFactors = FALSE
      )
    )
    if (length(distinct_values) <= 30L) {
      value_counts <- as.data.frame(table(as.character(column_values), useNA = "ifany"), stringsAsFactors = FALSE)
      colnames(value_counts) <- c("value", "cells")
      value_counts$slice <- slice_id
      value_counts$metadata_column <- metadata_column
      metadata_values <- rbind(metadata_values, value_counts[, c("slice", "metadata_column", "value", "cells")])
    }
  }
}
write.table(metadata_summary, file.path(result_dir, "metadata_schema_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(metadata_values, file.path(result_dir, "metadata_low_cardinality_value_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "00a_metadata_inspection_sessionInfo.txt"))
message("Metadata inspection completed.")
