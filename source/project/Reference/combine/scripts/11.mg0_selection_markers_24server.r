#!/usr/bin/env Rscript

# Display-only evidence for manual MG0 selection.  This script never changes
# the integrated Seurat object or its cluster labels.  The marker panel is
# specified by the investigator: TRIM2, RNF220, ST18, OLIG2, MBP and MOG.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(qs)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
stage2_root <- file.path(project_root, "combine/results/full_count_reanalysis/02_microglia_integration")
object_file <- file.path(stage2_root, "object/integrated_Microglia_processed.qs")
figure_dir <- file.path(stage2_root, "figure")
table_dir <- file.path(stage2_root, "table")
runtime_dir <- file.path(stage2_root, "runtime")
for (d in c(figure_dir, table_dir, runtime_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(object_file)) stop("Missing completed stage-2 object: ", object_file)

save_plot_both <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
    dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
    device = grDevices::pdf, bg = "white", limitsize = FALSE)
}

mg <- qs::qread(object_file)
DefaultAssay(mg) <- "RNA"
if (!"integrated_cluster" %in% colnames(mg[[]])) stop("integrated_cluster is absent from the completed stage-2 object.")
requested_markers <- c("TRIM2", "RNF220", "ST18", "OLIG2", "MBP", "MOG")
available_markers <- intersect(requested_markers, rownames(mg))
marker_availability <- data.frame(
  requested_marker = requested_markers,
  available_in_RNA = requested_markers %in% rownames(mg),
  stringsAsFactors = FALSE
)
utils::write.csv(marker_availability,
  file.path(table_dir, "MG0_selection_requested_marker_availability.csv"), row.names = FALSE)
if (!length(available_markers)) stop("None of the investigator-specified MG0 markers are available in the RNA feature axis.")
if (length(setdiff(requested_markers, available_markers))) {
  warning("Unavailable requested MG0 marker(s): ",
    paste(setdiff(requested_markers, available_markers), collapse = ", "))
}

cluster_counts <- as.data.frame(table(integrated_cluster = mg$integrated_cluster), stringsAsFactors = FALSE)
names(cluster_counts)[2] <- "n_cells"
utils::write.csv(cluster_counts, file.path(table_dir, "MG0_selection_integrated_cluster_counts.csv"), row.names = FALSE)

p_cluster <- DimPlot(mg, group.by = "integrated_cluster", label = TRUE, repel = TRUE) +
  ggtitle("MG0 selection: integrated Microglia Harmony UMAP")
p_dot <- DotPlot(mg, features = available_markers, group.by = "integrated_cluster") +
  RotatedAxis() + ggtitle("MG0 selection: investigator-specified marker dot plot")
p_feature <- FeaturePlot(mg, features = available_markers, ncol = 3,
  min.cutoff = "q05", max.cutoff = "q95") +
  plot_annotation(title = "MG0 selection: investigator-specified marker expression")

selection_plots <- list(
  MG0_selection_cluster_UMAP = list(plot = p_cluster, width = 11, height = 8),
  MG0_selection_marker_dotplot = list(plot = p_dot, width = 11, height = 7),
  MG0_selection_marker_UMAP = list(plot = p_feature, width = 13, height = 10)
)
for (plot_name in names(selection_plots)) {
  item <- selection_plots[[plot_name]]
  save_plot_both(item$plot, file.path(figure_dir, plot_name), item$width, item$height)
  print(item$plot)
}

writeLines(c(
  "Manual MG0 selection evidence only; the integrated object and cluster assignments are not changed.",
  paste0("Requested marker panel: ", paste(requested_markers, collapse = ", ")),
  paste0("Markers present in the integrated RNA feature axis: ", paste(available_markers, collapse = ", ")),
  "Dot size represents percent of cells expressing the gene; colour represents scaled average expression.",
  "Choose MG0 manually after reviewing the UMAP, dotplot and Feature UMAP."
), file.path(runtime_dir, "MG0_selection_marker_evidence_manifest.txt"))
message("MG0 marker evidence complete. Review the three MG0_selection_* figures and specify the chosen integrated cluster.")
