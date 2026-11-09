#!/usr/bin/env Rscript

# Spatial microglia: Human MG0 top-50 replacement-reference branch.
# The Human MG0 top-50 table was generated from Human MG0 versus Human MG1--MG5
# FindMarkers. The human symbols are matched to each MERFISH RNA assay without
# changing case-sensitive source feature names, then scored with the same
# AddModuleScore(ctrl = 10) and FeatureStatPlot workflow as spatial_module.r.

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/three_slice_ST_MG0_MG5_dual_module_score_20260916")
result_dir <- file.path(project_dir, "results", "02_Human_MG0_top50_module_score")
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
human_top50_tsv <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/OLM_MG0_top50_cross_species_20260915/results/01_human_MG0_top50_score_cross_species_OLM/tables/Human_MG0_FindMarkers_top50_reference_signature.tsv")
target_order <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")
st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35", "ST_MG1" = "#4DBBD5", "ST_MG2" = "#00A087",
  "ST_MG3" = "#3C5488", "ST_MG4" = "#F39B7F", "ST_MG5" = "#7E6148"
)

human_top50 <- read.delim(human_top50_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (!"gene" %in% colnames(human_top50) || nrow(human_top50) != 50L) {
  stop("The Human MG0 reference file must have 50 rows and a gene column.")
}
human_top50_genes <- as.character(human_top50$gene)
if (length(unique(human_top50_genes)) != 50L) stop("Human MG0 top-50 genes are not unique.")

all_slice_intersection <- data.frame()
all_slice_cluster_summary <- data.frame()
all_slice_comparisons <- data.frame()

for (slice_id in names(slice_paths)) {
  slice_object <- readRDS(slice_paths[[slice_id]])
  DefaultAssay(slice_object) <- "RNA"
  if (!"CellType" %in% colnames(slice_object[[]])) stop(slice_id, " lacks CellType.")
  slice_object$plot_celltype <- factor(as.character(slice_object$CellType), levels = target_order)
  if (any(is.na(slice_object$plot_celltype))) stop(slice_id, " contains a CellType outside ST_MG0--ST_MG5.")

  slice_RNA_genes <- rownames(slice_object[["RNA"]])
  slice_RNA_upper <- toupper(slice_RNA_genes)
  human_top50_upper <- toupper(human_top50_genes)
  matched_feature_index <- match(human_top50_upper, slice_RNA_upper)
  slice_human_top50_genes <- unique(slice_RNA_genes[matched_feature_index[!is.na(matched_feature_index)]])
  if (length(slice_human_top50_genes) == 0L) stop(slice_id, " has no case-normalized overlap with the Human MG0 top-50 signature.")
  slice_intersection_table <- data.frame(
    slice = slice_id,
    Human_MG0_top50_gene = human_top50_genes,
    case_normalized_symbol = human_top50_upper,
    matched_slice_RNA_gene = ifelse(is.na(matched_feature_index), NA_character_, slice_RNA_genes[matched_feature_index]),
    present_in_slice_RNA = !is.na(matched_feature_index),
    stringsAsFactors = FALSE
  )
  all_slice_intersection <- rbind(all_slice_intersection, slice_intersection_table)

  slice_object <- AddModuleScore(
    object = slice_object,
    assay = "RNA",
    features = list(slice_human_top50_genes),
    ctrl = 10,
    name = "Human_MG0_top50",
    seed = 20260916
  )
  score_column_from_AddModuleScore <- "Human_MG0_top501"
  if (!score_column_from_AddModuleScore %in% colnames(slice_object[[]])) {
    stop(slice_id, " did not produce ", score_column_from_AddModuleScore, ".")
  }
  score_column <- "Human_MG0_top50_score"
  slice_object[[score_column]] <- slice_object[[score_column_from_AddModuleScore]][, 1]

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
  slice_violin <- slice_violin + patchwork::plot_annotation(title = paste0(slice_id, ": Human MG0 top-50 module score"))
  ggsave(file.path(figure_dir, paste0(slice_id, "_Human_MG0_top50_violin.png")), slice_violin, width = 8.5, height = 6.2, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, paste0(slice_id, "_Human_MG0_top50_violin.pdf")), slice_violin, width = 8.5, height = 6.2, device = grDevices::cairo_pdf, bg = "white")

  score_values <- as.numeric(slice_object[[score_column]][, 1])
  celltype_values <- as.character(slice_object$plot_celltype)
  slice_cluster_summary <- do.call(rbind, lapply(comparison_groups, function(group_name) {
    group_scores <- score_values[celltype_values == group_name]
    data.frame(
      slice = slice_id,
      cluster = group_name,
      cells = length(group_scores),
      mean_Human_MG0_top50_score = mean(group_scores),
      median_Human_MG0_top50_score = median(group_scores),
      sd_Human_MG0_top50_score = sd(group_scores),
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
      ST_MG0_mean_Human_MG0_top50_score = mean(mg0_scores),
      comparison_mean_Human_MG0_top50_score = mean(other_scores),
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
    Human_MG0_top50_score = score_values,
    stringsAsFactors = FALSE
  )
  write.table(slice_cell_audit, file.path(table_dir, paste0(slice_id, "_Human_MG0_top50_per_cell.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(slice_object, file.path(object_dir, paste0(slice_id, "_ST_MG0_MG5_Human_MG0_top50_scored.rds")), compress = "gzip")
}

write.table(all_slice_intersection, file.path(table_dir, "Human_MG0_top50_marker_case_normalized_intersection_by_slice.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_slice_cluster_summary, file.path(table_dir, "Human_MG0_top50_summary_by_slice_and_cluster.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_slice_comparisons, file.path(table_dir, "Human_MG0_top50_Wilcoxon_ST_MG0_vs_other_clusters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parameter_record <- data.frame(
  item = c("reference_script", "replacement_marker_reference", "Human_marker_contrast", "feature_matching", "score_method", "score_assay", "control_genes", "score_column", "grouping", "reference_group", "violin_method", "test"),
  value = c(
    file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/scripts/spatial_module.r"),
    human_top50_tsv, "Human MG0 versus Human MG1-MG5 in RNA data", "case-normalized gene symbols matched to the native spatial RNA feature names", "Seurat::AddModuleScore", "RNA", "10", "Human_MG0_top50_score", "CellType (ST_MG0--ST_MG5)", "ST_MG0", "SCP::FeatureStatPlot", "two-sided Wilcoxon; raw p displayed as stars, BH p saved per slice"
  ),
  stringsAsFactors = FALSE
)
write.table(parameter_record, file.path(result_dir, "02_Human_MG0_top50_parameter_record.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "02_Human_MG0_top50_module_score_sessionInfo.txt"))
writeLines(as.character(Sys.time()), file.path(log_dir, "02_Human_MG0_top50_module_score_completed_at.txt"))
message("Human MG0 top-50 replacement-reference module-score branch completed.")
