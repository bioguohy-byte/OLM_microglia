#!/usr/bin/env Rscript
# Extract the exact ST_MG0-ST_MG5 cells from each full-slice MERFISH v2 object.
# Source objects are read only; no existing object is overwritten.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

spatial_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
source_dir <- file.path(spatial_root, "data/processed_data")
output_dir <- file.path(source_dir, "per_slice_microglia_ST_MG0_MG5")
audit_dir <- file.path(spatial_root, "results/per_slice_microglia_ST_MG0_MG5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

slices <- c("s2r1", "s2r2", "s2r3")
mg_labels <- paste0("ST_MG", 0:5)
audit_rows <- list()

for (slice_id in slices) {
  source_file <- file.path(source_dir, paste0("merfish_", slice_id, "_v2.rds"))
  output_file <- file.path(
    output_dir,
    paste0("merfish_", slice_id, "_microglia_ST_MG0_MG5_from_v2.rds")
  )

  if (!file.exists(source_file)) stop("Source object was not found: ", source_file)
  full_obj <- readRDS(source_file)
  if (!"CellType" %in% colnames(full_obj@meta.data)) {
    stop("CellType is absent from ", source_file)
  }

  keep_cells <- colnames(full_obj)[as.character(full_obj$CellType) %in% mg_labels]
  if (!length(keep_cells)) stop("No ST_MG0-ST_MG5 cells found in ", source_file)

  # Seurat subsetting retains the sparse RNA assay, metadata, and the spatial
  # reduction for only the selected cells. No new normalization or clustering.
  microglia_obj <- subset(full_obj, cells = keep_cells)
  validObject(microglia_obj)
  saveRDS(microglia_obj, output_file)

  label_counts <- table(factor(as.character(microglia_obj$CellType), levels = mg_labels))
  audit_rows[[slice_id]] <- tibble(
    slice_id = slice_id,
    source_object = source_file,
    output_object = output_file,
    source_cells = ncol(full_obj),
    selected_microglia_cells = ncol(microglia_obj),
    selected_label = names(label_counts),
    selected_label_cells = as.integer(label_counts),
    retained_assays = paste(names(microglia_obj@assays), collapse = ";"),
    retained_reductions = paste(names(microglia_obj@reductions), collapse = ";"),
    selection_rule = "CellType %in% ST_MG0,ST_MG1,ST_MG2,ST_MG3,ST_MG4,ST_MG5"
  )

  rm(full_obj, microglia_obj)
  gc(verbose = FALSE)
}

audit <- bind_rows(audit_rows)
total_selected_cells <- audit %>%
  distinct(slice_id, selected_microglia_cells) %>%
  summarise(total = sum(selected_microglia_cells)) %>%
  pull(total)
write_csv(audit, file.path(audit_dir, "per_slice_microglia_ST_MG0_MG5_audit.csv"))
writeLines(
  c(
    "Per-slice MERFISH microglia objects",
    "Selection rule: CellType %in% ST_MG0-ST_MG5 in the corresponding full-slice v2 object.",
    "No filtering, normalization, reclustering, or coordinate change was performed.",
    paste0("Total selected cells: ", total_selected_cells)
  ),
  file.path(audit_dir, "README.txt")
)

message("Created three per-slice microglia objects in: ", output_dir)
message("Total selected cells: ", total_selected_cells)
