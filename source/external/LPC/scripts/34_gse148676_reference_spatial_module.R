# MODULE : Cell Type Annotation and Module Scoring -------------------
#
# Direct adaptation of spatial_module.r for merged GSE148676 cells.
# The original package calls and plotting structure are retained. Only the
# input/output paths and the GSE148676 cluster mapping differ. The saveRDS()
# directly after AddModuleScore makes the original reload-for-plot step work
# when this script is sourced from a fresh R session.

## 1 Parameters --------------------------------------------------------------
# Load dependencies
library(Seurat)
library(qs)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SCP)

# File paths
project <- Sys.getenv("MG_LPC_ROOT")
file_in_merfish <- file.path(
  project,
  "objects/simple_microglia/mouse_GSE148676_merged_groups_seurat_v5.rds"
)
file_in_markers <- file.path(
  project,
  "raw_data/reference/mouse_mg_markers.qs"
)
file_out_merfish <- file.path(
  project,
  "objects/simple_microglia/mouse_GSE148676_merged_groups_reference_module_score.rds"
)
file_out_umap_pdf <- file.path(
  project,
  "figures/simple_microglia/mouse_GSE148676_reference_module_score_UMAP.pdf"
)
file_out_umap_png <- sub("\\.pdf$", ".png", file_out_umap_pdf)
file_out_violin_pdf <- file.path(
  project,
  "figures/simple_microglia/mouse_GSE148676_reference_module_score_violin.pdf"
)
file_out_violin_png <- sub("\\.pdf$", ".png", file_out_violin_pdf)
dir.create(dirname(file_out_merfish), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(file_out_umap_pdf), recursive = TRUE, showWarnings = FALSE)

# Cell-type mapping dictionary (old merged-cluster label -> plot label)
col_old_celltype <- "merged_cluster"
col_new_celltype <- "plot_celltype"
celltype_dict <- c(
  "0" = "cluster0",
  "1" = "cluster1",
  "2" = "cluster2",
  "3" = "cluster3",
  "4" = "cluster4",
  "5" = "cluster5",
  "6" = "cluster6",
  "7" = "cluster7",
  "8" = "cluster8",
  "9" = "cluster9",
  "10" = "cluster10"
)

# Marker extraction and module-score parameters
marker_cluster_col <- "cluster"
target_cluster <- "cluster2"
marker_gene_col <- "gene"

score_module_name <- "myelin"
score_feature_name <- paste0(score_module_name, "1")
score_ctrl_size <- 10

# Plotting parameter retained from the reference script
plot_pt_size <- 0


## 2 Process Data ------------------------------------------------------------
# 2.1 Read data
merfish <- readRDS(file_in_merfish)
mouse_mg_markers <- qread(file_in_markers)

# 2.2 Replace cell-type names
merfish@meta.data[[col_new_celltype]] <- as.character(
  merfish@meta.data[[col_old_celltype]]
)
for (old_name in names(celltype_dict)) {
  new_name <- celltype_dict[[old_name]]
  merfish@meta.data[[col_new_celltype]][
    merfish@meta.data[[col_new_celltype]] %in% old_name
  ] <- new_name
}

# 2.3 Extract target-cluster genes and intersect with GSE148676 genes
filtered_data <- mouse_mg_markers %>%
  filter(.data[[marker_cluster_col]] == target_cluster)
int_gene <- intersect(filtered_data[[marker_gene_col]], rownames(merfish))

# 2.4 Calculate module score
merfish <- AddModuleScore(
  object = merfish,
  features = list(int_gene),
  ctrl = score_ctrl_size,
  name = score_module_name
)

# The reference SCP function accesses the Seurat v4 RNA@data slot. This
# conversion preserves the existing normalized RNA layer for plotting only.
merfish[["RNA"]] <- as(merfish[["RNA"]], "Assay")

# The reference script re-opens a previously saved scored object for plotting.
# Save here so the same operation is valid in a new session.
saveRDS(merfish, file = file_out_merfish)


## 3 Generate Plots ----------------------------------------------------------
merfish <- readRDS(file_out_merfish)

st_mg_colors_npg <- c(
  "cluster1" = "#E64B35",
  "cluster0" = "#4DBBD5",
  "cluster2" = "#00A087",
  "cluster3" = "#3C5488",
  "cluster4" = "#F39B7F",
  "cluster5" = "#7E6148",
  "cluster6" = "#8E6BBE",
  "cluster7" = "#7F7F7F",
  "cluster8" = "#66C2A5",
  "cluster9" = "#E6AB02",
  "cluster10" = "#A65628"
)

(p1 <- DimPlot(merfish, group.by = col_new_celltype, cols = st_mg_colors_npg))

# The first group is the merged candidate group (old clusters 3, 9, and 10),
# preserving the reference comparison pattern of first group versus all others.
target_order <- c(
  "cluster1", "cluster0", "cluster2", "cluster3", "cluster4", "cluster5",
  "cluster6", "cluster7", "cluster8", "cluster9", "cluster10"
)
merfish$plot_celltype <- factor(merfish$plot_celltype, levels = target_order)

celltype <- levels(merfish$plot_celltype)
my_comparisons <- lapply(celltype[-1], function(x) c(x, celltype[1]))

p2 <- FeatureStatPlot(
  srt = merfish, group.by = "plot_celltype",
  stat.by = c(score_feature_name), add_box = FALSE,
  comparisons = my_comparisons,
  bg.by = "plot_celltype",
  palcolor = st_mg_colors_npg,
  bg_palcolor = st_mg_colors_npg
)


# Keep the reference layout object available in the session, but export its two
# component panels separately for figure use.
combined_plot <- p1 | p2


## 4 Save Results ------------------------------------------------------------
# Save updated Seurat object (with the new cell-type and module score)
saveRDS(merfish, file = file_out_merfish)

# Save standalone UMAP and violin panels
ggsave(filename = file_out_umap_pdf, plot = p1, width = 5, height = 5)
ggsave(filename = file_out_umap_png, plot = p1, width = 5, height = 5, dpi = 300)
ggsave(filename = file_out_violin_pdf, plot = p2, width = 7, height = 5)
ggsave(filename = file_out_violin_png, plot = p2, width = 7, height = 5, dpi = 300)
