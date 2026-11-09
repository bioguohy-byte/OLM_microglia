#!/usr/bin/env Rscript

# Read-only structural/provenance diagnosis for the historical Chia object.
# This script does not save a Seurat object or alter source metadata.

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: chia_readonly_integrity_diagnosis.R <Chia.qs> <output_dir>")
}
source_path <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- qread(source_path)
md <- obj@meta.data
cell_ids <- colnames(obj)

as_text <- function(x) {
  if (is.factor(x)) as.character(x) else as.character(x)
}

make_crosswalk <- function(x, y, x_name, y_name) {
  keep <- !is.na(x) & !is.na(y)
  out <- as.data.frame(table(as_text(x[keep]), as_text(y[keep])), stringsAsFactors = FALSE)
  colnames(out) <- c(x_name, y_name, "n_cells")
  out <- out[out$n_cells > 0L, , drop = FALSE]
  out[order(out[[x_name]], -out$n_cells), , drop = FALSE]
}

cluster_purity <- function(cluster, label) {
  keep <- !is.na(cluster) & !is.na(label)
  tab <- table(as_text(cluster[keep]), as_text(label[keep]))
  total <- rowSums(tab)
  top_n <- apply(tab, 1L, max)
  top_label <- apply(tab, 1L, function(z) colnames(tab)[which.max(z)])
  p <- sweep(tab, 1L, total, "/")
  entropy <- apply(p, 1L, function(z) {
    z <- z[z > 0]
    -sum(z * log2(z))
  })
  data.frame(
    SCT_snn_res_0_8 = rownames(tab),
    n_cells = as.integer(total),
    top_author_cell_type = top_label,
    top_author_fraction = as.numeric(top_n / total),
    author_label_entropy_bits = as.numeric(entropy),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

required <- c("SCT_snn_res.0.8", "RNA_snn_res.0.8", "seurat_clusters", "cell.type")
missing_cols <- setdiff(required, colnames(md))
if (length(missing_cols)) stop("Missing expected metadata column(s): ", paste(missing_cols, collapse = ", "))

umap <- Embeddings(obj, reduction = "umap")
rna_data <- GetAssayData(obj, assay = "RNA", layer = "data")
sct_data <- GetAssayData(obj, assay = "SCT", layer = "data")

marker_panel <- c(
  "AQP4","ITGB4","ETNPPL","SLC1A2","RFX4","GLI3","GRAMD2B","GFAP",
  "VCAN","C1QL1","SUSD5","COL9A1",
  "MOG","MOBP","SLC5A11","PICK1","HAPLN2","CNDP1","MYRF","CPM","TF","TMC6","CDK18",
  "PTPRC","CSF1R","MRC1","IL18","C3","SLC2A5","BIN2","DOCK2","ENTPD1",
  "CFAP54","DNAH11","DCN","LAMA2","RERGL","ABCG2","EMCN","LEF1","FLT1","SLC6A12","CFH","SLCO1B3",
  "RBFOX3","SNAP25","SYT1","GRIN1","CD247"
)
rna_feature_upper <- toupper(rownames(rna_data))
marker_present <- marker_panel %in% rna_feature_upper

summary_rows <- data.frame(
  check = c(
    "source_path", "n_cells", "n_features_RNA", "n_features_SCT", "assays", "default_assay",
    "reductions", "graphs", "metadata_rows_match_cells", "metadata_rownames_match_cells",
    "UMAP_rows_match_cells", "UMAP_rownames_match_cells", "UMAP_finite_coordinates",
    "duplicate_cell_barcodes", "duplicate_metadata_barcodes", "SCT_cluster_levels",
    "RNA_cluster_levels", "seurat_clusters_levels", "active_ident_matches_seurat_clusters",
    "Yao_marker_panel_present_in_RNA"
  ),
  value = c(
    source_path, ncol(obj), nrow(rna_data), nrow(sct_data), paste(Assays(obj), collapse = ";"), DefaultAssay(obj),
    paste(Reductions(obj), collapse = ";"), paste(Graphs(obj), collapse = ";"), nrow(md) == length(cell_ids),
    identical(rownames(md), cell_ids), nrow(umap) == length(cell_ids), identical(rownames(umap), cell_ids),
    all(is.finite(umap)), anyDuplicated(cell_ids), anyDuplicated(rownames(md)),
    length(unique(as_text(md[["SCT_snn_res.0.8"]]))), length(unique(as_text(md[["RNA_snn_res.0.8"]]))),
    length(unique(as_text(md[["seurat_clusters"]]))),
    identical(as_text(Idents(obj)), as_text(md[["seurat_clusters"]])),
    paste0(sum(marker_present), "/", length(marker_panel))
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_rows, file.path(out_dir, "chia_object_integrity_summary.csv"), row.names = FALSE)

cluster_counts <- rbind(
  transform(as.data.frame(table(as_text(md[["SCT_snn_res.0.8"]])), stringsAsFactors = FALSE), cluster_column = "SCT_snn_res.0.8"),
  transform(as.data.frame(table(as_text(md[["RNA_snn_res.0.8"]])), stringsAsFactors = FALSE), cluster_column = "RNA_snn_res.0.8"),
  transform(as.data.frame(table(as_text(md[["seurat_clusters"]])), stringsAsFactors = FALSE), cluster_column = "seurat_clusters")
)
colnames(cluster_counts) <- c("cluster", "n_cells", "cluster_column")
write.csv(cluster_counts, file.path(out_dir, "chia_cluster_count_by_metadata_column.csv"), row.names = FALSE)

write.csv(
  make_crosswalk(md[["SCT_snn_res.0.8"]], md[["seurat_clusters"]], "SCT_snn_res_0_8", "seurat_clusters"),
  file.path(out_dir, "chia_SCT25_vs_seurat_clusters_crosswalk.csv"), row.names = FALSE
)
write.csv(
  make_crosswalk(md[["SCT_snn_res.0.8"]], md[["RNA_snn_res.0.8"]], "SCT_snn_res_0_8", "RNA_snn_res_0_8"),
  file.path(out_dir, "chia_SCT25_vs_RNA_res08_crosswalk.csv"), row.names = FALSE
)
write.csv(
  make_crosswalk(md[["SCT_snn_res.0.8"]], md[["cell.type"]], "SCT_snn_res_0_8", "author_cell_type"),
  file.path(out_dir, "chia_SCT25_vs_author_cell_type_crosswalk.csv"), row.names = FALSE
)
write.csv(
  make_crosswalk(md[["seurat_clusters"]], md[["cell.type"]], "seurat_clusters", "author_cell_type"),
  file.path(out_dir, "chia_seurat_clusters_vs_author_cell_type_crosswalk.csv"), row.names = FALSE
)
write.csv(
  cluster_purity(md[["SCT_snn_res.0.8"]], md[["cell.type"]]),
  file.path(out_dir, "chia_SCT25_author_cell_type_purity.csv"), row.names = FALSE
)
write.csv(
  data.frame(marker = marker_panel, present_in_RNA = marker_present, stringsAsFactors = FALSE),
  file.path(out_dir, "chia_yao_marker_panel_presence.csv"), row.names = FALSE
)

writeLines(c(
  "# Chia historical object: read-only integrity diagnosis", "",
  "Source object was opened with qs::qread and not saved or modified.",
  "The report checks barcode and UMAP alignment, assay availability, clustering-column relationships, and agreement with existing author cell.type metadata.",
  "SCT_snn_res.0.8 is the 25-cluster field used in the initial marker audit; seurat_clusters is preserved separately so any mismatch is explicit.",
  "No CellType mapping, QC filtering, UMAP recomputation, or source-object mutation was performed."
), file.path(out_dir, "README.md"))

message("Read-only diagnosis complete: ", normalizePath(out_dir))
