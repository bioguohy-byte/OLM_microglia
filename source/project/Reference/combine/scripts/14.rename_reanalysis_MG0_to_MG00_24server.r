#!/usr/bin/env Rscript

# Rename ONLY the reanalysis manual cluster-2 label from MG0 to MG00.
# Original MG0--MG5 source labels remain unchanged because they document the
# earlier analysis.  New result objects and figures are written to MG00-named
# paths, so prior MG0-labelled outputs remain an auditable historical record.

rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
invisible(gc())

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(qs)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
old_stage2_root <- file.path(project_root, "combine/results/full_count_reanalysis/02_microglia_integration")
old_stage3_root <- file.path(project_root, "combine/results/full_count_reanalysis/03_current_MG0_doublet_diagnosis")
result_root <- file.path(project_root, "combine/results/full_count_reanalysis/02_microglia_integration_MG00")
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
object_dir <- file.path(result_root, "object")
runtime_dir <- file.path(result_root, "runtime")
for (directory in c(figure_dir, table_dir, object_dir, runtime_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

reference_file <- file.path(old_stage2_root, "object/integrated_Microglia_processed.qs")
diagnostic_file <- file.path(old_stage3_root, "object/current_MG0_QC_before_doublet_removal_three_method_partial.qs")
if (any(!file.exists(c(reference_file, diagnostic_file)))) {
  stop("Required reanalysis object is missing.")
}

save_plot_both <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
    dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
    device = grDevices::pdf, bg = "white", limitsize = FALSE)
}
remove_old_reanalysis_label_fields <- function(object) {
  for (field in intersect(c("MG0_status", "MG0_transfer_score"), colnames(object[[]]))) {
    object[[field]] <- NULL
  }
  object
}

message("Loading completed reanalysis integrated Microglia object ...")
integrated_MG00 <- qs::qread(reference_file)
if (!"integrated_cluster" %in% colnames(integrated_MG00[[]])) {
  stop("The integrated Microglia object lacks integrated_cluster.")
}
integrated_MG00 <- remove_old_reanalysis_label_fields(integrated_MG00)
integrated_MG00$MG00_status <- factor(
  ifelse(as.character(integrated_MG00$integrated_cluster) == "2", "MG00", "Other integrated clusters"),
  levels = c("Other integrated clusters", "MG00")
)
integrated_MG00$manual_MG_state <- ifelse(
  as.character(integrated_MG00$integrated_cluster) == "2", "MG00", NA_character_
)
integrated_MG00$manual_MG_state_note <- ifelse(
  as.character(integrated_MG00$integrated_cluster) == "2",
  "User-selected integrated cluster 2; renamed from reanalysis MG0 to MG00 (2026-09-07)",
  NA_character_
)
qs::qsave(
  integrated_MG00,
  file.path(object_dir, "integrated_Microglia_MG00_processed.qs"),
  preset = "high"
)

DefaultAssay(integrated_MG00) <- "RNA"
requested_markers <- c("TRIM2", "RNF220", "ST18", "OLIG2", "MBP", "MOG")
available_markers <- intersect(requested_markers, rownames(integrated_MG00))
if (!length(available_markers)) stop("None of the six requested MG00 markers is present.")
marker_availability <- data.frame(
  requested_marker = requested_markers,
  available_in_RNA = requested_markers %in% rownames(integrated_MG00),
  stringsAsFactors = FALSE
)
write.csv(marker_availability, file.path(table_dir, "MG00_requested_marker_availability.csv"), row.names = FALSE)
cluster_counts <- as.data.frame(table(integrated_cluster = integrated_MG00$integrated_cluster), stringsAsFactors = FALSE)
names(cluster_counts)[2] <- "n_cells"
cluster_counts$manual_label <- ifelse(cluster_counts$integrated_cluster == "2", "MG00", NA_character_)
write.csv(cluster_counts, file.path(table_dir, "MG00_integrated_cluster_counts.csv"), row.names = FALSE)

p_cluster <- DimPlot(integrated_MG00, reduction = "umap", group.by = "integrated_cluster", label = TRUE, repel = TRUE) +
  ggtitle("MG00 selection: integrated Microglia Harmony UMAP")
p_MG00 <- DimPlot(
  integrated_MG00, reduction = "umap", group.by = "MG00_status",
  cols = c("Other integrated clusters" = "#D9D9D9", MG00 = "#C51B7D")
) + ggtitle("MG00: user-selected integrated cluster 2")
p_dot <- DotPlot(integrated_MG00, features = available_markers, group.by = "integrated_cluster") +
  RotatedAxis() + ggtitle("MG00 selection: investigator-specified marker dot plot")
p_feature <- FeaturePlot(integrated_MG00, features = available_markers, ncol = 3,
  min.cutoff = "q05", max.cutoff = "q95") +
  plot_annotation(title = "MG00 selection: investigator-specified marker expression")
selection_plots <- list(
  MG00_selection_cluster_UMAP = list(plot = p_cluster, width = 11, height = 8),
  MG00_definition_UMAP = list(plot = p_MG00, width = 10, height = 7),
  MG00_selection_marker_dotplot = list(plot = p_dot, width = 11, height = 7),
  MG00_selection_marker_UMAP = list(plot = p_feature, width = 13, height = 10)
)
for (plot_name in names(selection_plots)) {
  item <- selection_plots[[plot_name]]
  save_plot_both(item$plot, file.path(figure_dir, plot_name), item$width, item$height)
  print(item$plot)
}

integrated_audit <- integrated_MG00[[]]
integrated_audit$Cell <- rownames(integrated_audit)
integrated_audit <- integrated_audit[, c("Cell", "Study", "integrated_cluster", "MG00_status", "manual_MG_state", "manual_MG_state_note"), drop = FALSE]
write.csv(integrated_audit, file.path(table_dir, "integrated_Microglia_MG00_cell_metadata.csv"), row.names = FALSE)
rm(integrated_MG00)
invisible(gc())

message("Loading pre-removal diagnostic object and renaming its transferred label ...")
diagnostic_MG00 <- qs::qread(diagnostic_file)
if (!all(c("transferred_integrated_cluster", "MG0_status", "MG0_transfer_score") %in% colnames(diagnostic_MG00[[]]))) {
  stop("The diagnostic object does not contain the prior MG0 transfer fields.")
}
old_status <- as.character(diagnostic_MG00$MG0_status)
old_score <- as.numeric(diagnostic_MG00$MG0_transfer_score)
diagnostic_MG00 <- remove_old_reanalysis_label_fields(diagnostic_MG00)
diagnostic_MG00$MG00_transfer_score <- old_score
diagnostic_MG00$MG00_status <- factor(
  ifelse(old_status == "MG0 (manual cluster 2 transfer)",
    "MG00 (manual cluster 2 transfer)", "Other predicted cluster"),
  levels = c("Other predicted cluster", "MG00 (manual cluster 2 transfer)")
)
qs::qsave(
  diagnostic_MG00,
  file.path(object_dir, "current_MG00_QC_before_doublet_removal_three_method_partial.qs"),
  preset = "high"
)
diagnostic_audit <- diagnostic_MG00[[]]
diagnostic_audit$Cell <- rownames(diagnostic_audit)
diagnostic_audit <- diagnostic_audit[, c(
  "Cell", "Study", "original_barcode", "transferred_integrated_cluster",
  "MG00_transfer_score", "MG00_status"
), drop = FALSE]
write.csv(diagnostic_audit, file.path(table_dir, "diagnostic_MG00_label_cell_metadata.csv"), row.names = FALSE)

writeLines(c(
  "Scope: rename the reanalysis manual integrated cluster-2 label from MG0 to MG00.",
  "Original CellType1 MG0--MG5 labels are historical source-pipeline labels and are not edited.",
  "New object: integrated_Microglia_MG00_processed.qs contains MG00_status and manual_MG_state.",
  "New object: current_MG00_QC_before_doublet_removal_three_method_partial.qs contains MG00_status and MG00_transfer_score.",
  "Prior MG0_status and MG0_transfer_score fields were removed from both newly written MG00 objects.",
  "No counts, QC flags, UMAP coordinates, clusters, or doublet calls changed."
), file.path(runtime_dir, "MG00_rename_manifest.txt"))
writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "MG00_rename_R_sessionInfo.txt"))
message("MG00 rename complete. New objects and four selection figures saved under: ", result_root)
