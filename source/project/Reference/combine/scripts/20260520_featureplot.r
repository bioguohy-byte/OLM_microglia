setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/"))
rm(list=ls())

# 0. Global Setup --------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(qs)
library(openxlsx)
library(future)
library(pheatmap)
library(ggplotify)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(patchwork)
library(stringr)

input_file      <- 'data/processed_data/mg_cca_meta_v2.qs'
output_fig_dir  <- 'results/figure/'
save_path   <- paste0(output_fig_dir, '20260520_featurePlot.pdf')

dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

gene_list <- c('TRIM2','RNF220','MBP','MOG')

plot_pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(10)
tmp_obj <- qread(input_file)

plot <- FeaturePlot(tmp_obj,features = gene_list,cols = plot_pal, raster = FALSE)
ggsave(save_path, plot, width = 7.5, height = 6, limitsize = FALSE)
