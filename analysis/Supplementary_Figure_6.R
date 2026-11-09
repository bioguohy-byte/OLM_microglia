#!/usr/bin/env Rscript
# Figure-ordered code from the reviewed manuscript analysis package.
# Each panel contains its original quoted R code.
# Source this chapter, then use mg_list_panels() or mg_show_panel("PANEL_ID").
# Source paths in panel records identify the original analysis scripts.
.mg_panels <- new.env(parent=emptyenv())
.mg_order <- character()
mg_register_panel <- function(id, code) {
  if (!id %in% .mg_order) .mg_order <<- c(.mg_order,id)
  assign(id, code, envir=.mg_panels)
}
mg_list_panels <- function() .mg_order
mg_panel_code <- function(id) {
  if (!exists(id,envir=.mg_panels,inherits=FALSE)) stop("Unknown panel: ",id)
  get(id,envir=.mg_panels,inherits=FALSE)
}
mg_show_panel <- function(id) {
  blocks <- mg_panel_code(id)
  if (is.null(blocks)) {cat(id, ": No R block is included for this panel.\n");return(invisible(NULL))}
  for (block in blocks) {
    cat("# Source: ",block$source," [actual ",block$ranges,"; requested ",block$requested_ranges,"]\n# Role: ",block$role,"\n",sep="")
    cat(paste(deparse(block$code,width.cutoff=120L),collapse="\n"),"\n\n")
  }
  invisible(blocks)
}


# Supplementary Figure 6

# S6A  Whole s2r2 slice: OLG, OPC, ST-OLM and other cells
# Kind: computational
mg_register_panel("S6A", list(
  list(source="source/project/Reference/combine/scripts/Spatial/scripts/full_slice_four_colour_center_coordinates_1to1.R", ranges="10-60", requested_ranges="10-57", role="Shared three-slice script; this panel uses s2r2", code=quote({
project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
data_dir <- file.path(project_dir, "data/processed_data")
out_dir <- file.path(project_dir, "results/figure/code_verification_full_slice_four_colour_center_coordinates_1to1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_colours <- c(
  "Other cells" = "grey80",
  "OLG" = "blue",
  "OPC" = "yellow",
  "ST_MG0" = "red"
)
plot_levels <- names(plot_colours)

for (sample_id in c("s2r1", "s2r2", "s2r3")) {
  obj <- readRDS(file.path(data_dir, paste0("merfish_", sample_id, "_v2.rds")))
  stopifnot(all(c("CellType", "center_x", "center_y") %in% colnames(obj@meta.data)))

  # The original full-slice display labels horizontal axis as y and vertical axis
  # as x. No rotation, mirroring, or use of the Seurat spatial reduction occurs.
  df <- obj@meta.data %>%
    transmute(
      center_x = as.numeric(center_x),
      center_y = as.numeric(center_y),
      four_colour_celltype = if_else(
        CellType %in% c("OLG", "OPC", "ST_MG0"), as.character(CellType), "Other cells"
      )
    ) %>%
    mutate(four_colour_celltype = factor(four_colour_celltype, levels = plot_levels))

  p <- ggplot(df, aes(x = center_y, y = center_x, colour = four_colour_celltype)) +
    geom_point(size = 0.25) +
    scale_colour_manual(values = plot_colours, drop = FALSE) +
    coord_fixed(ratio = 1, expand = FALSE) +
    labs(
      title = paste0("Whole-slice raw-coordinate plot: ", sample_id),
      subtitle = "Original center_y (horizontal) and center_x (vertical); 1:1 coordinate scale",
      x = "y (center_y)", y = "x (center_x)", colour = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )

  stem <- file.path(out_dir, paste0(sample_id, "_whole_slice_OLG_OPC_ST_MG0_other_center_coordinates_1to1"))
  ggsave(paste0(stem, ".pdf"), p, width = 8, height = 8, bg = "white")
  ggsave(paste0(stem, ".png"), p, width = 8, height = 8, dpi = 300, bg = "white")
  rm(obj, df, p)
  gc()
}
  }))
))

# S6B  s2r1 microglial UMAP
# Kind: computational
mg_register_panel("S6B", list(
  list(source="source/project/Reference/combine/scripts/Spatial/scripts/run_per_slice_microglia_umap_from_integrated_logic.R", ranges="15-105", requested_ranges="15-44;46-80", role="Shared three-slice script; this panel uses s2r1", code=quote({
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
  }))
))

# S6C  s2r1 myelin1 module scores
# Kind: computational
mg_register_panel("S6C", list(
  list(source="source/project/Reference/combine/scripts/Spatial/three_slice_ST_MG0_MG5_dual_module_score_20260916/scripts/01_mouse_cluster2_module_score_per_slice.R", ranges="1-176", requested_ranges="1-176", role="Shared three-slice script; this panel uses s2r1", code=quote({
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
  }))
))

# S6D  Whole s2r1 slice: OLG, OPC, ST-OLM and other cells
# Kind: computational
mg_register_panel("S6D", list(
  list(source="source/project/Reference/combine/scripts/Spatial/scripts/full_slice_four_colour_center_coordinates_1to1.R", ranges="10-60", requested_ranges="10-57", role="Shared three-slice script; this panel uses s2r1", code=quote({
project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
data_dir <- file.path(project_dir, "data/processed_data")
out_dir <- file.path(project_dir, "results/figure/code_verification_full_slice_four_colour_center_coordinates_1to1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_colours <- c(
  "Other cells" = "grey80",
  "OLG" = "blue",
  "OPC" = "yellow",
  "ST_MG0" = "red"
)
plot_levels <- names(plot_colours)

for (sample_id in c("s2r1", "s2r2", "s2r3")) {
  obj <- readRDS(file.path(data_dir, paste0("merfish_", sample_id, "_v2.rds")))
  stopifnot(all(c("CellType", "center_x", "center_y") %in% colnames(obj@meta.data)))

  # The original full-slice display labels horizontal axis as y and vertical axis
  # as x. No rotation, mirroring, or use of the Seurat spatial reduction occurs.
  df <- obj@meta.data %>%
    transmute(
      center_x = as.numeric(center_x),
      center_y = as.numeric(center_y),
      four_colour_celltype = if_else(
        CellType %in% c("OLG", "OPC", "ST_MG0"), as.character(CellType), "Other cells"
      )
    ) %>%
    mutate(four_colour_celltype = factor(four_colour_celltype, levels = plot_levels))

  p <- ggplot(df, aes(x = center_y, y = center_x, colour = four_colour_celltype)) +
    geom_point(size = 0.25) +
    scale_colour_manual(values = plot_colours, drop = FALSE) +
    coord_fixed(ratio = 1, expand = FALSE) +
    labs(
      title = paste0("Whole-slice raw-coordinate plot: ", sample_id),
      subtitle = "Original center_y (horizontal) and center_x (vertical); 1:1 coordinate scale",
      x = "y (center_y)", y = "x (center_x)", colour = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )

  stem <- file.path(out_dir, paste0(sample_id, "_whole_slice_OLG_OPC_ST_MG0_other_center_coordinates_1to1"))
  ggsave(paste0(stem, ".pdf"), p, width = 8, height = 8, bg = "white")
  ggsave(paste0(stem, ".png"), p, width = 8, height = 8, dpi = 300, bg = "white")
  rm(obj, df, p)
  gc()
}
  }))
))

# S6E  s2r3 microglial UMAP
# Kind: computational
mg_register_panel("S6E", list(
  list(source="source/project/Reference/combine/scripts/Spatial/scripts/run_per_slice_microglia_umap_from_integrated_logic.R", ranges="15-105", requested_ranges="15-44;46-80", role="Shared three-slice script; this panel uses s2r3", code=quote({
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
  }))
))

# S6F  s2r3 myelin1 module scores
# Kind: computational
mg_register_panel("S6F", list(
  list(source="source/project/Reference/combine/scripts/Spatial/three_slice_ST_MG0_MG5_dual_module_score_20260916/scripts/01_mouse_cluster2_module_score_per_slice.R", ranges="1-176", requested_ranges="1-176", role="Shared three-slice script; this panel uses s2r3", code=quote({
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
  }))
))

# S6G  Whole s2r3 slice: OLG, OPC, ST-OLM and other cells
# Kind: computational
mg_register_panel("S6G", list(
  list(source="source/project/Reference/combine/scripts/Spatial/scripts/full_slice_four_colour_center_coordinates_1to1.R", ranges="10-60", requested_ranges="10-57", role="Shared three-slice script; this panel uses s2r3", code=quote({
project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
data_dir <- file.path(project_dir, "data/processed_data")
out_dir <- file.path(project_dir, "results/figure/code_verification_full_slice_four_colour_center_coordinates_1to1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_colours <- c(
  "Other cells" = "grey80",
  "OLG" = "blue",
  "OPC" = "yellow",
  "ST_MG0" = "red"
)
plot_levels <- names(plot_colours)

for (sample_id in c("s2r1", "s2r2", "s2r3")) {
  obj <- readRDS(file.path(data_dir, paste0("merfish_", sample_id, "_v2.rds")))
  stopifnot(all(c("CellType", "center_x", "center_y") %in% colnames(obj@meta.data)))

  # The original full-slice display labels horizontal axis as y and vertical axis
  # as x. No rotation, mirroring, or use of the Seurat spatial reduction occurs.
  df <- obj@meta.data %>%
    transmute(
      center_x = as.numeric(center_x),
      center_y = as.numeric(center_y),
      four_colour_celltype = if_else(
        CellType %in% c("OLG", "OPC", "ST_MG0"), as.character(CellType), "Other cells"
      )
    ) %>%
    mutate(four_colour_celltype = factor(four_colour_celltype, levels = plot_levels))

  p <- ggplot(df, aes(x = center_y, y = center_x, colour = four_colour_celltype)) +
    geom_point(size = 0.25) +
    scale_colour_manual(values = plot_colours, drop = FALSE) +
    coord_fixed(ratio = 1, expand = FALSE) +
    labs(
      title = paste0("Whole-slice raw-coordinate plot: ", sample_id),
      subtitle = "Original center_y (horizontal) and center_x (vertical); 1:1 coordinate scale",
      x = "y (center_y)", y = "x (center_x)", colour = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )

  stem <- file.path(out_dir, paste0(sample_id, "_whole_slice_OLG_OPC_ST_MG0_other_center_coordinates_1to1"))
  ggsave(paste0(stem, ".pdf"), p, width = 8, height = 8, bg = "white")
  ggsave(paste0(stem, ".png"), p, width = 8, height = 8, dpi = 300, bg = "white")
  rm(obj, df, p)
  gc()
}
  }))
))

# S6H  OLM proportion in three observations
# Kind: experimental_quantification
mg_register_panel("S6H", list(
  list(source="source/recovered/experimental/S4D_OLM_proportion.R", ranges="3-5;7-8;11-15;17-58;60-60;62-76;78-79", requested_ranges="1-79", role="Original three-value OLM plot", code=quote({
suppressPackageStartupMessages({
  library(ggplot2)
})

output_dir <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "S4D")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

olm_data <- data.frame(
  group = "OLM",
  value = c(0.08244, 0.10910, 0.08240),
  point_x = c(0.94, 1.00, 1.06)
)

plot_olm <- ggplot(olm_data, aes(x = 1, y = value)) +
  geom_boxplot(
    aes(group = 1),
    width = 0.38,
    outlier.shape = NA,
    fill = "#1786C0",
    alpha = 0.35,
    color = "#1786C0",
    linewidth = 0.9,
    staplewidth = 0.55
  ) +
  geom_point(
    aes(x = point_x),
    shape = 21,
    size = 2.6,
    stroke = 0.9,
    fill = "#1786C0",
    color = "black"
  ) +
  scale_x_continuous(
    breaks = 1,
    labels = "OLM",
    limits = c(0.5, 1.5),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    name = expression(Trim2^"+" * Iba1^"+" * OLM/Iba1^"+" ~ MG),
    breaks = seq(0, 0.12, 0.02),
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL) +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(linewidth = 0.9),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    plot.margin = margin(8, 8, 7, 7)
  ) +
  coord_cartesian(ylim = c(0, 0.12))

print(plot_olm)

ggsave(
  file.path(output_dir, "OLM_三值箱线图.png"),
  plot_olm,
  width = 3.0,
  height = 4.3,
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(output_dir, "OLM_三值箱线图.pdf"),
  plot_olm,
  width = 3.0,
  height = 4.3,
  device = cairo_pdf
)

cat("输入值（%）：", paste(olm_data$value * 100, collapse = ", "), "\n")
cat("输出目录：", output_dir, "\n", sep = "")
  }))
))
