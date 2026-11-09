#!/usr/bin/env Rscript
# Figure-specific grouping of the original human microglia plotting code.
suppressPackageStartupMessages({library(Seurat);library(qs);library(dplyr);library(ggplot2);library(forcats);library(patchwork);library(ggplotify);library(tidydr);library(RColorBrewer);library(pheatmap)})
input_obj_v2 <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs")
output_root <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "S3B", "plots")
figure_dir <- file.path(output_root, "pdf")
png_dir <- file.path(output_root, "png")
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(png_dir, recursive=TRUE, showWarnings=FALSE)
group_col <- "Region"
type_col <- "CellType1"
tmp_obj <- qread(input_obj_v2)

p_group_dim <- DimPlot(
  tmp_obj,
  reduction = "umap",
  group.by = group_col,
  raster = FALSE,
  label = FALSE,
  pt.size = 0.1
) +
  ggtitle(group_col) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2.5))) +
  tidydr::theme_dr(
    xlength = 0.2,
    ylength = 0.2,
    arrow = grid::arrow(length = grid::unit(0.15, "inches"), type = "closed")
  ) &
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none",
    axis.title.x = element_text(hjust = 0),
    axis.title.y = element_text(hjust = 0)
  )

p_group_split <- DimPlot(tmp_obj, group.by = group_col, split.by = group_col) +
  NoAxes() +
  theme(plot.title = element_blank())

combined_group_dim <- (as.ggplot(p_group_dim) | as.ggplot(p_group_split)) +
  plot_layout(widths = c(1, 6))

ggsave(
  file.path(figure_dir, paste0("UMAP_Split_", group_col, ".pdf")),
  combined_group_dim,
  width = 18,
  height = 3
)
ggsave(file.path(png_dir, paste0("UMAP_Split_", group_col, ".png")), combined_group_dim, width=18, height=3, dpi=300)

