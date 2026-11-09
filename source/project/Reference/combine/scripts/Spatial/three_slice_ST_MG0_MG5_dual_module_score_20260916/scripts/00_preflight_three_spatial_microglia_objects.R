#!/usr/bin/env Rscript

# Preflight for three already selected spatial microglia objects.
# No QC, filtering, clustering, coordinate recalculation or cell-type relabelling
# is done in this script. It only verifies the fixed inputs required for the two
# module-score branches to follow.

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/three_slice_ST_MG0_MG5_dual_module_score_20260916")
result_dir <- file.path(project_dir, "results", "00_preflight")
log_dir <- file.path(project_dir, "logs")
tmp_dir <- file.path(project_dir, "tmp")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir, XDG_CACHE_HOME = file.path(project_dir, "cache"))
dir.create(Sys.getenv("XDG_CACHE_HOME"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library"), file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
})

slice_paths <- c(
  s2r1 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r1_microglia_ST_MG0_MG5_from_v2.rds"),
  s2r2 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r2_microglia_ST_MG0_MG5_from_v2.rds"),
  s2r3 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r3_microglia_ST_MG0_MG5_from_v2.rds")
)
expected_cells <- c(s2r1 = 2994L, s2r2 = 3060L, s2r3 = 3578L)
mouse_marker_qs <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/saved_data/mouse_mg_markers.qs")
human_top50_tsv <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/OLM_MG0_top50_cross_species_20260915/results/01_human_MG0_top50_score_cross_species_OLM/tables/Human_MG0_FindMarkers_top50_reference_signature.tsv")

if (!all(file.exists(slice_paths))) stop("One or more slice RDS files are missing.")
if (!file.exists(mouse_marker_qs)) stop("Mouse cluster2 marker QS is missing: ", mouse_marker_qs)
if (!file.exists(human_top50_tsv)) stop("Human MG0 top-50 marker TSV is missing: ", human_top50_tsv)

mouse_markers <- qread(mouse_marker_qs)
if (!all(c("cluster", "gene") %in% colnames(mouse_markers))) {
  stop("Mouse marker QS must contain cluster and gene columns.")
}
mouse_cluster2 <- mouse_markers[as.character(mouse_markers$cluster) == "cluster2", , drop = FALSE]
if (nrow(mouse_cluster2) == 0L) stop("Mouse marker QS does not contain cluster2 markers.")

human_top50 <- read.delim(human_top50_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (!"gene" %in% colnames(human_top50) || nrow(human_top50) != 50L || length(unique(human_top50$gene)) != 50L) {
  stop("Human MG0 reference table must contain exactly 50 unique genes in a gene column.")
}

object_summary <- data.frame()
celltype_summary <- data.frame()

for (slice_id in names(slice_paths)) {
  slice_object <- readRDS(slice_paths[[slice_id]])
  if (!inherits(slice_object, "Seurat")) stop(slice_id, " is not a Seurat object.")
  if (ncol(slice_object) != expected_cells[[slice_id]]) {
    stop(slice_id, " has ", ncol(slice_object), " cells; expected ", expected_cells[[slice_id]], ".")
  }
  if (!"RNA" %in% Assays(slice_object)) stop(slice_id, " lacks an RNA assay.")
  metadata_names <- colnames(slice_object[[]])
  available_celltype_column <- if ("plot_celltype" %in% metadata_names) {
    "plot_celltype"
  } else if ("CellType" %in% metadata_names) {
    "CellType"
  } else if ("cell.type_manual" %in% metadata_names) {
    "cell.type_manual"
  } else {
    NA_character_
  }
  if (is.na(available_celltype_column)) {
    stop(slice_id, " lacks plot_celltype, CellType and cell.type_manual metadata.")
  }
  object_summary <- rbind(
    object_summary,
    data.frame(
      slice = slice_id,
      input_rds = slice_paths[[slice_id]],
      cells = ncol(slice_object),
      features_RNA = nrow(slice_object[["RNA"]]),
      assays = paste(Assays(slice_object), collapse = ";"),
      reductions = paste(Reductions(slice_object), collapse = ";"),
      default_assay = DefaultAssay(slice_object),
      celltype_column_for_analysis = available_celltype_column,
      preexisting_myelin1 = "myelin1" %in% metadata_names,
      stringsAsFactors = FALSE
    )
  )
  celltype_values <- as.character(slice_object[[available_celltype_column]][, 1])
  celltype_counts <- as.data.frame(table(celltype_values, useNA = "ifany"), stringsAsFactors = FALSE)
  colnames(celltype_counts) <- c("celltype", "cells")
  celltype_counts$slice <- slice_id
  celltype_counts$source_metadata_column <- available_celltype_column
  celltype_summary <- rbind(celltype_summary, celltype_counts[, c("slice", "source_metadata_column", "celltype", "cells")])
}

write.table(object_summary, file.path(result_dir, "three_slice_object_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(celltype_summary, file.path(result_dir, "three_slice_celltype_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mouse_cluster2, file.path(result_dir, "mouse_cluster2_reference_markers_full.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(human_top50, file.path(result_dir, "Human_MG0_top50_reference_markers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parameter_record <- data.frame(
  item = c("original_spatial_reference", "mouse_reference", "human_reference", "planned_module_score_method", "planned_control_gene_number", "planned_violin_method", "cell_retention"),
  value = c(
    file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/scripts/spatial_module.r"),
    "mouse_mg_markers cluster2 genes intersected with the RNA assay",
    "Human MG0 versus Human MG1-MG5 FindMarkers top-50 genes; table copied from the completed cross-species analysis",
    "Seurat::AddModuleScore",
    "10",
    "SCP::FeatureStatPlot with ST_MG0 first and each other ST_MG cluster compared to it",
    "No cells removed; each source slice is retained intact"
  ),
  stringsAsFactors = FALSE
)
write.table(parameter_record, file.path(result_dir, "preflight_parameter_record.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "00_preflight_sessionInfo.txt"))
writeLines(as.character(Sys.time()), file.path(log_dir, "00_preflight_completed_at.txt"))
message("Preflight completed for ", length(slice_paths), " spatial microglia slice objects.")
