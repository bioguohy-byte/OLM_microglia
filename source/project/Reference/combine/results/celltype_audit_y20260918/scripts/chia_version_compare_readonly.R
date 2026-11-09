#!/usr/bin/env Rscript

# Read-only comparison of the two historical Chia full-cell objects.

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: chia_version_compare_readonly.R <Chia.qs> <Chia_v2.qs> <output_dir>")
old_path <- args[[1]]
v2_path <- args[[2]]
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

describe_object <- function(obj, version) {
  md <- obj@meta.data
  candidate_columns <- intersect(c("CellType", "cell.type", "SCT_snn_res.0.8", "RNA_snn_res.0.8", "seurat_clusters"), colnames(md))
  stats <- data.frame(
    object_version = version,
    n_cells = ncol(obj),
    assays = paste(Assays(obj), collapse = ";"),
    reductions = paste(Reductions(obj), collapse = ";"),
    default_assay = DefaultAssay(obj),
    graph_names = paste(Graphs(obj), collapse = ";"),
    metadata_columns = paste(colnames(md), collapse = ";"),
    stringsAsFactors = FALSE
  )
  labels <- do.call(rbind, lapply(candidate_columns, function(k) {
    vals <- as.character(md[[k]])
    data.frame(
      object_version = version,
      metadata_column = k,
      n_levels = length(unique(vals)),
      level_values = paste(sort(unique(vals)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  list(stats = stats, labels = labels, cells = colnames(obj))
}

v2 <- qread(v2_path)
v2_info <- describe_object(v2, "Chia_v2")
v2_celltype <- if ("CellType" %in% colnames(v2@meta.data)) as.character(v2$CellType) else NULL
v2_umap <- if ("umap" %in% Reductions(v2)) Embeddings(v2, reduction = "umap") else NULL
rm(v2); gc()

old <- qread(old_path)
old_info <- describe_object(old, "Chia")
old_celltype <- if ("cell.type" %in% colnames(old@meta.data)) as.character(old$cell.type) else NULL
old_umap <- if ("umap" %in% Reductions(old)) Embeddings(old, reduction = "umap") else NULL
same_celltype <- !is.null(v2_celltype) && !is.null(old_celltype) && identical(v2_celltype, old_celltype)
umap_same_dimensions <- !is.null(v2_umap) && !is.null(old_umap) && identical(dim(v2_umap), dim(old_umap))
umap_max_abs_difference <- if (umap_same_dimensions) max(abs(v2_umap - old_umap)) else NA_real_
umap_exactly_identical <- !is.null(v2_umap) && !is.null(old_umap) && identical(v2_umap, old_umap)
rm(old); gc()

object_summary <- rbind(old_info$stats, v2_info$stats)
label_summary <- rbind(old_info$labels, v2_info$labels)
cell_set_summary <- data.frame(
  comparison = c("same_cell_count", "identical_cell_order", "same_cell_set", "only_in_Chia", "only_in_Chia_v2",
                 "same_author_celltype_values", "UMAP_same_dimensions", "UMAP_exactly_identical", "UMAP_max_abs_coordinate_difference"),
  value = c(
    length(old_info$cells) == length(v2_info$cells),
    identical(old_info$cells, v2_info$cells),
    setequal(old_info$cells, v2_info$cells),
    length(setdiff(old_info$cells, v2_info$cells)),
    length(setdiff(v2_info$cells, old_info$cells)),
    same_celltype, umap_same_dimensions, umap_exactly_identical, umap_max_abs_difference
  ),
  stringsAsFactors = FALSE
)
write.csv(object_summary, file.path(out_dir, "chia_object_version_summary.csv"), row.names = FALSE)
write.csv(label_summary, file.path(out_dir, "chia_object_version_metadata_summary.csv"), row.names = FALSE)
write.csv(cell_set_summary, file.path(out_dir, "chia_object_version_cell_set_comparison.csv"), row.names = FALSE)
writeLines(c(
  "# Chia historical object-version comparison (read-only)", "",
  "Both source .qs files were opened read-only. No object was saved or modified.",
  "This comparison reports only object structure, candidate label/cluster columns, and cell barcode agreement."
), file.path(out_dir, "VERSION_COMPARISON_README.md"))
