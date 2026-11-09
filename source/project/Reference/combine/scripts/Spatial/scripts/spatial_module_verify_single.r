# Verification copy of spatial_module.r: single-plot reproduction only.
# The original merged object retained in this migration already contains
# plot_celltype and myelin1. No label, score, or source object is changed.

scp_lib <- file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library")
.libPaths(c(scp_lib, .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(SCP)
})

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
file_in_merfish <- file.path(project_dir, "data/processed_data/merfish_immune.combined_with4.rds")
dir_out_fig <- file.path(project_dir, "results/figure/code_verification_module")
dir.create(dir_out_fig, recursive = TRUE, showWarnings = FALSE)

# Original script settings retained verbatim.
col_new_celltype <- "plot_celltype"
score_feature_name <- "myelin1"
st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35",
  "ST_MG1" = "#4DBBD5",
  "ST_MG2" = "#00A087",
  "ST_MG3" = "#3C5488",
  "ST_MG4" = "#F39B7F",
  "ST_MG5" = "#7E6148"
)
target_order <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")

merfish <- readRDS(file_in_merfish)
stopifnot(all(c(col_new_celltype, score_feature_name) %in% colnames(merfish@meta.data)))
merfish$plot_celltype <- factor(merfish$plot_celltype, levels = target_order)

# Blue-box panel D: original DimPlot call, saved as a standalone figure.
p_umap <- DimPlot(merfish, group.by = col_new_celltype, cols = st_mg_colors_npg) +
  ggtitle("Spatial MG subclusters")

# Blue-box panel E: original FeatureStatPlot call, saved as a standalone figure.
celltype <- levels(merfish$plot_celltype)
my_comparisons <- lapply(celltype[-1], function(x) c(x, celltype[1]))
p_myelin <- SCP::FeatureStatPlot(
  srt = merfish,
  group.by = "plot_celltype",
  stat.by = c(score_feature_name),
  add_box = FALSE,
  comparisons = my_comparisons,
  bg.by = "plot_celltype",
  palcolor = st_mg_colors_npg,
  bg_palcolor = st_mg_colors_npg
)

ggsave(file.path(dir_out_fig, "spatial_MG_subclusters_UMAP_verified.pdf"), p_umap,
       width = 7, height = 6, bg = "white")
ggsave(file.path(dir_out_fig, "spatial_MG_subclusters_UMAP_verified.png"), p_umap,
       width = 7, height = 6, dpi = 300, bg = "white")
ggsave(file.path(dir_out_fig, "myelin_module_score_FeatureStatPlot_verified.pdf"), p_myelin,
       width = 7, height = 6, bg = "white")
ggsave(file.path(dir_out_fig, "myelin_module_score_FeatureStatPlot_verified.png"), p_myelin,
       width = 7, height = 6, dpi = 300, bg = "white")

message("Single-plot verification outputs: ", dir_out_fig)
