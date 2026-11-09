#!/usr/bin/env Rscript

# Spatial microglia: original reference branch.
# This reproduces the module-score part of spatial_module.r for each supplied
# ST_MG0--ST_MG5 slice, using mouse_mg_markers cluster2 and AddModuleScore
# with ctrl = 10. No input object is overwritten.

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/three_slice_ST_MG0_MG5_dual_module_score_20260916")
result_dir <- file.path(project_dir, "results", "01_mouse_cluster2_module_score")
figure_dir <- file.path(result_dir, "figures")
table_dir <- file.path(result_dir, "tables")
object_dir <- file.path(result_dir, "objects")
log_dir <- file.path(project_dir, "logs")
tmp_dir <- file.path(project_dir, "tmp")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir, XDG_CACHE_HOME = file.path(project_dir, "cache"))
dir.create(Sys.getenv("XDG_CACHE_HOME"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library"), file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(SCP)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260916)

slice_paths <- c(
  s2r1 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r1_microglia_ST_MG0_MG5_from_v2.rds"),
  s2r2 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r2_microglia_ST_MG0_MG5_from_v2.rds"),
  s2r3 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r3_microglia_ST_MG0_MG5_from_v2.rds")
)
mouse_marker_qs <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/saved_data/mouse_mg_markers.qs")
target_order <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")
st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35", "ST_MG1" = "#4DBBD5", "ST_MG2" = "#00A087",
  "ST_MG3" = "#3C5488", "ST_MG4" = "#F39B7F", "ST_MG5" = "#7E6148"
)

mouse_mg_markers <- qread(mouse_marker_qs)
mouse_cluster2_markers <- mouse_mg_markers[as.character(mouse_mg_markers$cluster) == "cluster2", , drop = FALSE]
mouse_cluster2_genes <- unique(as.character(mouse_cluster2_markers$gene))
if (length(mouse_cluster2_genes) == 0L) stop("No cluster2 genes were found in mouse_mg_markers.")

all_slice_intersection <- data.frame()
all_slice_cluster_summary <- data.frame()
all_slice_comparisons <- data.frame()

for (slice_id in names(slice_paths)) {
  slice_object <- readRDS(slice_paths[[slice_id]])
  DefaultAssay(slice_object) <- "RNA"
  if (!"CellType" %in% colnames(slice_object[[]])) stop(slice_id, " lacks CellType.")
  slice_object$plot_celltype <- factor(as.character(slice_object$CellType), levels = target_order)
  if (any(is.na(slice_object$plot_celltype))) stop(slice_id, " contains a CellType outside ST_MG0--ST_MG5.")
  if (!all(target_order %in% levels(droplevels(slice_object$plot_celltype)))) stop(slice_id, " lacks one or more ST_MG groups.")

  slice_gene_intersection <- intersect(mouse_cluster2_genes, rownames(slice_object[["RNA"]]))
  if (length(slice_gene_intersection) == 0L) stop(slice_id, " has no overlap with mouse cluster2 markers.")
  slice_intersection_table <- data.frame(
    slice = slice_id,
    mouse_cluster2_marker = mouse_cluster2_genes,
    present_in_slice_RNA = mouse_cluster2_genes %in% rownames(slice_object[["RNA"]]),
    stringsAsFactors = FALSE
  )
  all_slice_intersection <- rbind(all_slice_intersection, slice_intersection_table)

  slice_object <- AddModuleScore(
    object = slice_object,
    assay = "RNA",
    features = list(slice_gene_intersection),
    ctrl = 10,
    name = "myelin",
    seed = 20260916
  )
  score_column <- "myelin1"
  if (!score_column %in% colnames(slice_object[[]])) stop(slice_id, " did not produce myelin1.")

  comparison_groups <- levels(slice_object$plot_celltype)
  my_comparisons <- lapply(comparison_groups[-1], function(x) c(x, comparison_groups[1]))
  slice_violin <- FeatureStatPlot(
    srt = slice_object,
    assay = "RNA",
    slot = "data",
    group.by = "plot_celltype",
    stat.by = score_column,
    add_box = FALSE,
    comparisons = my_comparisons,
    bg.by = "plot_celltype",
    palcolor = st_mg_colors_npg,
    bg_palcolor = st_mg_colors_npg,
    legend.position = "right",
    xlab = "",
    theme_args = list(base_size = 10)
  )
  slice_violin <- slice_violin & theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.margin = ggplot2::margin(t = 16, r = 4, b = 4, l = 4, unit = "pt")
  )
  slice_violin <- slice_violin + patchwork::plot_annotation(title = paste0(slice_id, ": mouse cluster2 module score"))
  ggsave(file.path(figure_dir, paste0(slice_id, "_mouse_cluster2_myelin1_violin.png")), slice_violin, width = 8.5, height = 6.2, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, paste0(slice_id, "_mouse_cluster2_myelin1_violin.pdf")), slice_violin, width = 8.5, height = 6.2, device = grDevices::cairo_pdf, bg = "white")

  score_values <- as.numeric(slice_object[[score_column]][, 1])
  celltype_values <- as.character(slice_object$plot_celltype)
  slice_cluster_summary <- do.call(rbind, lapply(comparison_groups, function(group_name) {
    group_scores <- score_values[celltype_values == group_name]
    data.frame(
      slice = slice_id,
      cluster = group_name,
      cells = length(group_scores),
      mean_myelin1 = mean(group_scores),
      median_myelin1 = median(group_scores),
      sd_myelin1 = sd(group_scores),
      stringsAsFactors = FALSE
    )
  }))
  all_slice_cluster_summary <- rbind(all_slice_cluster_summary, slice_cluster_summary)

  slice_comparison_statistics <- do.call(rbind, lapply(comparison_groups[-1], function(comparison_cluster) {
    mg0_scores <- score_values[celltype_values == "ST_MG0"]
    other_scores <- score_values[celltype_values == comparison_cluster]
    test_result <- wilcox.test(other_scores, mg0_scores, alternative = "two.sided", exact = FALSE)
    data.frame(
      slice = slice_id,
      reference_cluster = "ST_MG0",
      comparison_cluster = comparison_cluster,
      ST_MG0_cells = length(mg0_scores),
      comparison_cells = length(other_scores),
      ST_MG0_mean_myelin1 = mean(mg0_scores),
      comparison_mean_myelin1 = mean(other_scores),
      wilcoxon_W = unname(test_result$statistic),
      p_value_unadjusted = test_result$p.value,
      stringsAsFactors = FALSE
    )
  }))
  slice_comparison_statistics$p_value_BH_within_slice <- p.adjust(slice_comparison_statistics$p_value_unadjusted, method = "BH")
  slice_comparison_statistics$plot_significance <- ifelse(slice_comparison_statistics$p_value_unadjusted <= 1e-4, "****", ifelse(slice_comparison_statistics$p_value_unadjusted <= 1e-3, "***", ifelse(slice_comparison_statistics$p_value_unadjusted <= 1e-2, "**", ifelse(slice_comparison_statistics$p_value_unadjusted <= 0.05, "*", "ns"))))
  all_slice_comparisons <- rbind(all_slice_comparisons, slice_comparison_statistics)

  slice_cell_audit <- data.frame(
    cell_barcode = colnames(slice_object),
    slice = slice_id,
    CellType = as.character(slice_object$CellType),
    plot_celltype = as.character(slice_object$plot_celltype),
    myelin1 = score_values,
    stringsAsFactors = FALSE
  )
  write.table(slice_cell_audit, file.path(table_dir, paste0(slice_id, "_mouse_cluster2_myelin1_per_cell.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(slice_object, file.path(object_dir, paste0(slice_id, "_ST_MG0_MG5_mouse_cluster2_myelin1_scored.rds")), compress = "gzip")
}

write.table(all_slice_intersection, file.path(table_dir, "mouse_cluster2_marker_intersection_by_slice.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_slice_cluster_summary, file.path(table_dir, "mouse_cluster2_myelin1_summary_by_slice_and_cluster.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_slice_comparisons, file.path(table_dir, "mouse_cluster2_myelin1_Wilcoxon_ST_MG0_vs_other_clusters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parameter_record <- data.frame(
  item = c("reference_script", "marker_reference", "reference_cluster", "score_method", "score_assay", "control_genes", "score_column", "grouping", "reference_group", "violin_method", "test"),
  value = c(
    file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/scripts/spatial_module.r"),
    mouse_marker_qs, "cluster2", "Seurat::AddModuleScore", "RNA", "10", "myelin1", "CellType (ST_MG0--ST_MG5)", "ST_MG0", "SCP::FeatureStatPlot", "two-sided Wilcoxon; raw p displayed as stars, BH p saved per slice"
  ),
  stringsAsFactors = FALSE
)
write.table(parameter_record, file.path(result_dir, "01_mouse_cluster2_parameter_record.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "01_mouse_cluster2_module_score_sessionInfo.txt"))
writeLines(as.character(Sys.time()), file.path(log_dir, "01_mouse_cluster2_module_score_completed_at.txt"))
message("Mouse cluster2 module-score branch completed.")
