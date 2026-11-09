#!/usr/bin/env Rscript

# Three mouse MERFISH microglia slices: RNA FindAllMarkers by existing
# CellType labels ST_MG0--ST_MG5. Each slice is analysed independently.  Only
# genes present in the native 483-gene MERFISH panel can appear in the result.
# Input RDS objects are read only and never overwritten.

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/three_slice_ST_MG0_MG5_dual_module_score_20260916")
result_dir <- file.path(project_dir, "results", "04_three_slice_RNA_FindAllMarkers")
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
  library(ggplot2)
})

set.seed(20260916)

slice_paths <- c(
  s2r1 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r1_microglia_ST_MG0_MG5_from_v2.rds"),
  s2r2 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r2_microglia_ST_MG0_MG5_from_v2.rds"),
  s2r3 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r3_microglia_ST_MG0_MG5_from_v2.rds")
)
target_order <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")
st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35", "ST_MG1" = "#4DBBD5", "ST_MG2" = "#00A087",
  "ST_MG3" = "#3C5488", "ST_MG4" = "#F39B7F", "ST_MG5" = "#7E6148"
)

all_slice_markers <- data.frame()
all_slice_top_markers <- data.frame()
all_slice_cell_counts <- data.frame()
normalization_audit <- data.frame()

for (slice_id in names(slice_paths)) {
  slice_object <- readRDS(slice_paths[[slice_id]])
  DefaultAssay(slice_object) <- "RNA"
  if (!"CellType" %in% colnames(slice_object[[]])) stop(slice_id, " lacks the CellType metadata column.")
  if (nrow(slice_object[["RNA"]]) != 483L) stop(slice_id, " does not contain the expected 483-gene MERFISH RNA panel.")

  slice_object$plot_celltype <- factor(as.character(slice_object$CellType), levels = target_order)
  if (any(is.na(slice_object$plot_celltype))) stop(slice_id, " contains a CellType outside ST_MG0--ST_MG5.")
  if (!all(target_order %in% unique(as.character(slice_object$plot_celltype)))) stop(slice_id, " lacks one or more required ST_MG groups.")
  Idents(slice_object) <- "plot_celltype"

  rna_data <- tryCatch(GetAssayData(slice_object, assay = "RNA", slot = "data"), error = function(error_message) NULL)
  if (is.null(rna_data) || ncol(rna_data) != ncol(slice_object) || nrow(rna_data) == 0L) {
    slice_object <- NormalizeData(slice_object, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    normalization_action <- "RNA data layer absent; LogNormalize(scale.factor=10000) calculated in analysis copy"
  } else {
    normalization_action <- "existing RNA data layer used; no renormalization"
  }
  normalization_audit <- rbind(normalization_audit, data.frame(slice = slice_id, normalization_action = normalization_action, stringsAsFactors = FALSE))

  slice_cell_counts <- as.data.frame(table(slice_object$plot_celltype), stringsAsFactors = FALSE)
  colnames(slice_cell_counts) <- c("cluster", "cells")
  slice_cell_counts$slice <- slice_id
  slice_cell_counts <- slice_cell_counts[, c("slice", "cluster", "cells")]
  all_slice_cell_counts <- rbind(all_slice_cell_counts, slice_cell_counts)
  write.table(slice_cell_counts, file.path(table_dir, paste0(slice_id, "_CellType_cell_counts.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  slice_markers <- FindAllMarkers(
    object = slice_object,
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    min.pct = 0.10,
    logfc.threshold = 0.25,
    only.pos = FALSE,
    verbose = FALSE
  )
  slice_markers$gene <- rownames(slice_markers)
  rownames(slice_markers) <- NULL
  slice_markers$slice <- slice_id
  marker_columns <- c("slice", "cluster", "gene", setdiff(colnames(slice_markers), c("slice", "cluster", "gene")))
  slice_markers <- slice_markers[, marker_columns]
  write.table(slice_markers, file.path(table_dir, paste0(slice_id, "_FindAllMarkers_RNA.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.csv(slice_markers, file.path(table_dir, paste0(slice_id, "_FindAllMarkers_RNA.csv")), row.names = FALSE, quote = TRUE)
  all_slice_markers <- rbind(all_slice_markers, slice_markers)

  logfc_column <- intersect(c("avg_log2FC", "avg_logFC"), colnames(slice_markers))[1]
  if (is.na(logfc_column)) stop(slice_id, " FindAllMarkers output lacks an average log-fold-change column.")
  slice_positive_markers <- slice_markers[slice_markers[[logfc_column]] > 0, , drop = FALSE]
  slice_top_markers <- do.call(rbind, lapply(target_order, function(cluster_name) {
    cluster_markers <- slice_positive_markers[as.character(slice_positive_markers$cluster) == cluster_name, , drop = FALSE]
    cluster_markers <- cluster_markers[order(cluster_markers$p_val_adj, -cluster_markers[[logfc_column]], cluster_markers$gene), , drop = FALSE]
    head(cluster_markers, n = 20L)
  }))
  slice_top_markers$rank_within_cluster <- ave(
    seq_len(nrow(slice_top_markers)),
    slice_top_markers$cluster,
    FUN = seq_along
  )
  write.table(slice_top_markers, file.path(table_dir, paste0(slice_id, "_FindAllMarkers_RNA_top20_positive_by_cluster.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.csv(slice_top_markers, file.path(table_dir, paste0(slice_id, "_FindAllMarkers_RNA_top20_positive_by_cluster.csv")), row.names = FALSE, quote = TRUE)
  all_slice_top_markers <- rbind(all_slice_top_markers, slice_top_markers)

  dotplot_markers <- do.call(rbind, lapply(target_order, function(cluster_name) {
    cluster_markers <- slice_top_markers[as.character(slice_top_markers$cluster) == cluster_name, , drop = FALSE]
    head(cluster_markers, n = 5L)
  }))
  dotplot_genes <- unique(as.character(dotplot_markers$gene))
  dotplot_genes <- dotplot_genes[dotplot_genes %in% rownames(slice_object[["RNA"]])]
  if (length(dotplot_genes) > 0L) {
    marker_dotplot <- DotPlot(slice_object, features = dotplot_genes, assay = "RNA", group.by = "plot_celltype", cols = c("#B8B8B8", "#B2182B")) +
      RotatedAxis() +
      ggtitle(paste0(slice_id, ": top positive RNA FindAllMarkers genes")) +
      theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9)
      )
    ggsave(file.path(figure_dir, paste0(slice_id, "_FindAllMarkers_top_positive_DotPlot.png")), marker_dotplot, width = 10, height = 7.5, dpi = 300, bg = "white")
    ggsave(file.path(figure_dir, paste0(slice_id, "_FindAllMarkers_top_positive_DotPlot.pdf")), marker_dotplot, width = 10, height = 7.5, device = grDevices::cairo_pdf, bg = "white")
  }

  saveRDS(slice_object, file.path(object_dir, paste0(slice_id, "_ST_MG0_MG5_RNA_FindAllMarkers_ready.rds")), compress = "gzip")
}

write.table(all_slice_markers, file.path(table_dir, "three_slices_RNA_FindAllMarkers_combined.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(all_slice_markers, file.path(table_dir, "three_slices_RNA_FindAllMarkers_combined.csv"), row.names = FALSE, quote = TRUE)
write.table(all_slice_top_markers, file.path(table_dir, "three_slices_RNA_FindAllMarkers_top20_positive_by_cluster_combined.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(all_slice_top_markers, file.path(table_dir, "three_slices_RNA_FindAllMarkers_top20_positive_by_cluster_combined.csv"), row.names = FALSE, quote = TRUE)
write.table(all_slice_cell_counts, file.path(table_dir, "three_slices_CellType_cell_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(normalization_audit, file.path(table_dir, "three_slices_FindAllMarkers_normalization_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parameter_record <- data.frame(
  item = c("input_objects", "species", "assay", "identity_column", "groups", "marker_test", "marker_layer", "min_pct", "logfc_threshold", "only_positive", "spatial_panel", "input_mutation"),
  value = c(
    paste(slice_paths, collapse = ";"), "Mouse", "RNA", "CellType", "ST_MG0--ST_MG5",
    "Seurat::FindAllMarkers Wilcoxon", "RNA data (log-normalized)", "0.10", "0.25", "FALSE; full up/down table exported; top20 table is positive markers only", "483 genes", "None; original RDS read only; separate analysis-copy RDS written"
  ),
  stringsAsFactors = FALSE
)
write.table(parameter_record, file.path(result_dir, "04_three_slice_RNA_FindAllMarkers_parameter_record.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "04_three_slice_RNA_FindAllMarkers_sessionInfo.txt"))
writeLines(as.character(Sys.time()), file.path(log_dir, "04_three_slice_RNA_FindAllMarkers_completed_at.txt"))
message("Three-slice RNA FindAllMarkers analysis completed.")
