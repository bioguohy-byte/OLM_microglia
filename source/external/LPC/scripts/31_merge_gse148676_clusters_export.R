#!/usr/bin/env Rscript

# GSE148676 microglia: merge selected Seurat clusters without changing the
# original seurat_clusters column, then export group-level expression and UMAP.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(Seurat)
})

project <- Sys.getenv("MG_LPC_ROOT")
input_object <- file.path(
  project,
  "objects",
  "simple_microglia",
  "mouse_GSE148676_simple_microglia_seurat_v5.rds"
)
output_object <- file.path(
  project,
  "objects",
  "simple_microglia",
  "mouse_GSE148676_merged_groups_seurat_v5.rds"
)
report_dir <- file.path(project, "reports", "simple_microglia")
figure_dir <- file.path(project, "figures", "simple_microglia")

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

merged_3_9_10 <- c("3", "9", "10")
merged_0_1_2_4_5_8 <- c("0", "1", "2", "4", "5", "8")
all_expected_clusters <- as.character(0:17)
remaining_clusters <- setdiff(
  all_expected_clusters,
  c(merged_3_9_10, merged_0_1_2_4_5_8)
)
new_cluster_levels <- as.character(0:10)
old_to_new_cluster <- c(
  # Final IDs are ordered by merged-group cell count, from largest to smallest.
  setNames(rep("1", length(merged_3_9_10)), merged_3_9_10),
  setNames(rep("0", length(merged_0_1_2_4_5_8)), merged_0_1_2_4_5_8),
  setNames(as.character(2:10), remaining_clusters)
)

message("Reading ", input_object)
microglia <- readRDS(input_object)
DefaultAssay(microglia) <- "RNA"

original_cluster <- as.character(microglia$seurat_clusters)
unexpected_clusters <- setdiff(unique(original_cluster), all_expected_clusters)
if (length(unexpected_clusters) > 0) {
  stop(
    "Unexpected original cluster labels: ",
    paste(unexpected_clusters, collapse = ", ")
  )
}

merged_cluster_source <- original_cluster
merged_cluster_source[original_cluster %in% merged_3_9_10] <- "3_9_10"
merged_cluster_source[original_cluster %in% merged_0_1_2_4_5_8] <- "0_1_2_4_5_8"
microglia$merged_cluster <- factor(
  unname(old_to_new_cluster[original_cluster]),
  levels = new_cluster_levels
)
microglia$merged_cluster_source <- merged_cluster_source
Idents(microglia) <- "merged_cluster"

definition_table <- data.table(
  original_cluster = all_expected_clusters,
  merged_cluster = unname(old_to_new_cluster[all_expected_clusters]),
  merged_cluster_source = ifelse(
    all_expected_clusters %in% merged_3_9_10,
    "3_9_10",
    ifelse(
      all_expected_clusters %in% merged_0_1_2_4_5_8,
      "0_1_2_4_5_8",
      all_expected_clusters
    )
  )
)
definition_table[, merged_cluster_order := as.integer(merged_cluster)]
definition_table[, original_cluster_order := as.integer(original_cluster)]
setorder(definition_table, merged_cluster_order, original_cluster_order)
definition_table[, c("merged_cluster_order", "original_cluster_order") := NULL]
fwrite(
  definition_table,
  file.path(report_dir, "mouse_GSE148676_merged_groups_definition.tsv"),
  sep = "\t"
)

group_counts <- data.table(
  merged_cluster = as.character(microglia$merged_cluster)
)[, .(n_cells = .N), by = merged_cluster]
group_counts[, merged_cluster_order := as.integer(merged_cluster)]
setorder(group_counts, merged_cluster_order)
group_counts[, merged_cluster_order := NULL]
fwrite(
  group_counts,
  file.path(report_dir, "mouse_GSE148676_merged_groups_cell_counts.tsv"),
  sep = "\t"
)
message("Merged group cell counts:")
print(group_counts)

# All source UMAP coordinates are retained; only the display identity changes.
umap_plot <- DimPlot(
  microglia,
  reduction = "umap",
  group.by = "merged_cluster",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.18,
  cols = c(
    "0" = "#B2182B",
    "1" = "#2166AC",
    "2" = "#4DAF4A",
    "3" = "#984EA3",
    "4" = "#FF7F00",
    "5" = "#A65628",
    "6" = "#F781BF",
    "7" = "#999999",
    "8" = "#66C2A5",
    "9" = "#E6AB02",
    "10" = "#7570B3"
  )
) +
  labs(
    title = "GSE148676 microglia: merged and renumbered clusters",
    color = "Merged cluster"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.key.height = grid::unit(0.45, "cm"),
    legend.text = element_text(size = 8)
  )

png_path <- file.path(
  figure_dir,
  "mouse_GSE148676_merged_groups_UMAP.png"
)
pdf_path <- file.path(
  figure_dir,
  "mouse_GSE148676_merged_groups_UMAP.pdf"
)
ggsave(png_path, umap_plot, width = 12, height = 8, dpi = 300, bg = "white")
ggsave(
  pdf_path,
  umap_plot,
  width = 12,
  height = 8,
  device = grDevices::cairo_pdf,
  bg = "white"
)

# Calculate mean log-normalized expression for every merged group without
# converting the full sparse gene-by-cell layer into a dense matrix.
normalized_expression <- LayerData(microglia[["RNA"]], layer = "data")
cell_group <- factor(
  as.character(microglia$merged_cluster),
  levels = new_cluster_levels
)
membership <- sparseMatrix(
  i = seq_len(ncol(microglia)),
  j = as.integer(cell_group),
  x = 1,
  dims = c(ncol(microglia), length(new_cluster_levels)),
  dimnames = list(colnames(microglia), new_cluster_levels)
)
mean_expression <- normalized_expression %*% membership
mean_expression <- sweep(
  mean_expression,
  MARGIN = 2,
  STATS = as.numeric(colSums(membership)),
  FUN = "/"
)

mean_expression_table <- data.table(gene = rownames(mean_expression))
mean_expression_table <- cbind(
  mean_expression_table,
  as.data.table(as.matrix(mean_expression))
)
fwrite(
  mean_expression_table,
  file.path(
    report_dir,
    "mouse_GSE148676_merged_groups_average_log_normalized_expression.csv"
  )
)

# Cell-level mapping exposes original clusters, merged clusters, metadata, and
# existing UMAP coordinates. The full cell-level matrix remains in output_object.
cell_metadata <- as.data.table(microglia[[]], keep.rownames = "cell")
cell_metadata[, original_cluster := original_cluster[match(cell, colnames(microglia))]]
cell_metadata[, merged_cluster := as.character(
  microglia$merged_cluster[match(cell, colnames(microglia))]
)]
cell_metadata[, merged_cluster_source := microglia$merged_cluster_source[
  match(cell, colnames(microglia))
]]
umap_coordinates <- as.data.table(Embeddings(microglia, reduction = "umap"), keep.rownames = "cell")
setnames(umap_coordinates, 2:3, c("UMAP_1", "UMAP_2"))
cell_metadata <- merge(cell_metadata, umap_coordinates, by = "cell", all.x = TRUE, sort = FALSE)
setcolorder(
  cell_metadata,
  c(
    "cell",
    "original_cluster",
    "merged_cluster",
    "merged_cluster_source",
    "UMAP_1",
    "UMAP_2"
  )
)
fwrite(
  cell_metadata,
  file.path(report_dir, "mouse_GSE148676_merged_groups_cell_metadata.csv")
)

saveRDS(microglia, output_object, compress = FALSE)

message("Wrote ", png_path)
message("Wrote ", pdf_path)
message("Wrote merged-group average expression matrix and cell metadata to ", report_dir)
message("Wrote full Seurat object with merged_cluster to ", output_object)
