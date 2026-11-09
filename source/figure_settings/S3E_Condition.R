#!/usr/bin/env Rscript
# Figure-specific grouping of the original human microglia plotting code.
suppressPackageStartupMessages({library(Seurat);library(qs);library(dplyr);library(ggplot2);library(forcats);library(patchwork);library(ggplotify);library(tidydr);library(RColorBrewer);library(pheatmap)})
input_obj_v2 <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs")
output_root <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "S3E", "plots")
figure_dir <- file.path(output_root, "pdf")
png_dir <- file.path(output_root, "png")
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(png_dir, recursive=TRUE, showWarnings=FALSE)
# Fig.09: ROE Heatmap ----------------------------------------------------------
var1 <- "CellType1"
var2 <- "Condition"
roe_output_path <- file.path(figure_dir, paste0("ROE_", var1, "_", var2, ".pdf"))

tmp_obj <- qread(input_obj_v2)

tab <- table(tmp_obj[[var1, drop = TRUE]], tmp_obj[[var2, drop = TRUE]])

if ("Unknown" %in% colnames(tab)) {
  tab <- tab[, colnames(tab) != "Unknown", drop = FALSE]
}

row_sum <- rowSums(tab)
col_sum <- colSums(tab)
total <- sum(tab)
expected <- outer(row_sum, col_sum) / total
roe <- as.matrix(tab / expected)

roe_plot_pal <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(100))

pdf(roe_output_path, width = 6, height = 6)
pheatmap(
  roe,
  display_numbers = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = roe_plot_pal,
  number_format = "%.1f",
  number_color = "white",
  cellwidth = 35,
  cellheight = 25,
  main = paste0("ROE: ", var1, " vs ", var2),
  angle_col = "45"
)
dev.off()

message("Fig.09 done.")
png(file.path(png_dir, paste0("ROE_", var1, "_", var2, ".png")),width=6,height=6,units="in",res=300)
pheatmap(
  roe,
  display_numbers = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = roe_plot_pal,
  number_format = "%.1f",
  number_color = "white",
  cellwidth = 35,
  cellheight = 25,
  main = paste0("ROE: ", var1, " vs ", var2),
  angle_col = "45"
)
dev.off()

