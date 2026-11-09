#!/usr/bin/env Rscript
# Per-slice microglia UMAPs using the recovered original integrated-object logic.
# The original combined object used an integrated assay. Each single-slice object
# has only RNA, so standard RNA normalization is performed before the identical
# PCA(50) -> neighbors(PC1-10) -> clusters(0.1:2) -> UMAP(PC1-10) sequence.
# Input objects are read only; processed copies are saved separately.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

spatial_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
input_dir <- file.path(spatial_root, "data/processed_data/per_slice_microglia_ST_MG0_MG5")
output_dir <- file.path(spatial_root, "data/processed_data/per_slice_microglia_ST_MG0_MG5_RNA_PCA_UMAP")
figure_dir <- file.path(spatial_root, "results/per_slice_microglia_ST_MG0_MG5_RNA_PCA_UMAP")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Recovered from merfish_immune.combined_with4.rds command log.
seed_use <- 42L
n_pcs <- 50L
umap_dims <- 1:10
cluster_resolutions <- seq(0.1, 2, 0.1)
target_order <- paste0("ST_MG", 0:5)
st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35",
  "ST_MG1" = "#4DBBD5",
  "ST_MG2" = "#00A087",
  "ST_MG3" = "#3C5488",
  "ST_MG4" = "#F39B7F",
  "ST_MG5" = "#7E6148"
)

slices <- c("s2r1", "s2r2", "s2r3")
audit_rows <- list()

for (slice_id in slices) {
  input_file <- file.path(input_dir, paste0("merfish_", slice_id, "_microglia_ST_MG0_MG5_from_v2.rds"))
  output_file <- file.path(output_dir, paste0("merfish_", slice_id, "_microglia_RNA_PCA50_UMAP_PC1_10.rds"))
  pdf_file <- file.path(figure_dir, paste0(slice_id, "_microglia_ST_MG0_MG5_UMAP.pdf"))
  png_file <- file.path(figure_dir, paste0(slice_id, "_microglia_ST_MG0_MG5_UMAP.png"))

  if (!file.exists(input_file)) stop("Input object not found: ", input_file)
  obj <- readRDS(input_file)
  if (!"CellType" %in% colnames(obj@meta.data)) stop("CellType is absent: ", input_file)
  if (!all(as.character(obj$CellType) %in% target_order)) {
    stop("Input contains a non-ST_MG0-ST_MG5 label: ", input_file)
  }

  # There is no integrated assay in a single brain slice; RNA is the valid
  # within-slice analysis assay. The original per-slice object is not modified.
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = min(2000L, nrow(obj)), verbose = FALSE)
  obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)

  # Exact downstream dimensionality and random seeds recovered from the
  # original integrated object's Seurat command log.
  pca_features <- VariableFeatures(obj)
  n_pcs_actual <- min(n_pcs, length(pca_features), ncol(obj) - 1L)
  obj <- RunPCA(obj, features = pca_features, npcs = n_pcs_actual,
                verbose = FALSE, seed.use = seed_use)
  dims_actual <- seq_len(min(max(umap_dims), ncol(Embeddings(obj, "pca"))))
  obj <- FindNeighbors(obj, reduction = "pca", dims = dims_actual, verbose = FALSE)
  obj <- FindClusters(obj, resolution = cluster_resolutions,
                      random.seed = 0L, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = dims_actual,
                 seed.use = seed_use, verbose = FALSE)

  # The original UMAP drawing call is retained; plot_celltype is a copied,
  # ordered alias of the original CellType labels for compatibility with it.
  obj$plot_celltype <- factor(as.character(obj$CellType), levels = target_order)
  p_umap <- DimPlot(obj, group.by = "plot_celltype", cols = st_mg_colors_npg)
  ggsave(pdf_file, p_umap, width = 7, height = 6, bg = "white")
  ggsave(png_file, p_umap, width = 7, height = 6, dpi = 300, bg = "white")
  validObject(obj)
  saveRDS(obj, output_file)

  label_counts <- table(factor(as.character(obj$CellType), levels = target_order))
  audit_rows[[slice_id]] <- tibble(
    slice_id = slice_id,
    input_object = input_file,
    output_object = output_file,
    figure_pdf = pdf_file,
    cells = ncol(obj),
    genes = nrow(obj),
    variable_features = length(pca_features),
    pca_dimensions = ncol(Embeddings(obj, "pca")),
    umap_dimensions = ncol(Embeddings(obj, "umap")),
    selected_label = names(label_counts),
    selected_label_cells = as.integer(label_counts),
    normalization = "LogNormalize; scale.factor=10000",
    pca = paste0("RNA; all variable features; npcs=", n_pcs_actual, "; seed=", seed_use),
    neighbors = paste0("pca; dims=", paste(dims_actual, collapse = ",")),
    clusters = "resolution=0.1,0.2,...,2.0; random.seed=0",
    umap = paste0("pca; dims=", paste(dims_actual, collapse = ","), "; seed=", seed_use,
                  "; Seurat defaults: uwot, n.neighbors=30, min.dist=0.3, metric=cosine")
  )

  rm(obj)
  gc(verbose = FALSE)
}

audit <- bind_rows(audit_rows)
write_csv(audit, file.path(figure_dir, "per_slice_microglia_RNA_PCA_UMAP_audit.csv"))
writeLines(
  c(
    "Three independent brain-slice microglia UMAP analyses.",
    "Input: existing per-slice ST_MG0-ST_MG5 objects, one object per MERFISH brain slice.",
    "No filtering or doublet calling was performed: this is an embedding-only analysis of the already extracted cells.",
    "Single-slice objects lack the combined object's integrated assay; RNA preprocessing is therefore used before the recovered PCA/neighbor/cluster/UMAP parameters.",
    paste0("Seurat version: ", as.character(packageVersion("Seurat")))
  ),
  file.path(figure_dir, "README.txt")
)

message("Completed per-slice microglia PCA/UMAP analysis: ", figure_dir)
