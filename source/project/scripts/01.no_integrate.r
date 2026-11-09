# no integrate------------------------------
setwd(Sys.getenv("MG_PROJECT_ROOT"))
rm(list=ls())

# 0. library --------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Seurat)
library(Matrix)
library(cowplot)
library(reshape2)
library(SingleCellExperiment)
library(umap)
library(ggthemes)
library(dplyr)
library(patchwork)
library(ggsci)
library(scales)
library(tidyverse)
library(magick)
library(magrittr)
library(devtools)
library(ComplexHeatmap)
library(ggplotify)
library(circlize)
library(pheatmap)

# 1. load data----------------
data.filt <- readRDS('data/processed_data/02.24samples_after_filter_data.filt.rds')
seu.obj <- data.filt

# 2. seurat 流程---------------
seu.obj <- NormalizeData(seu.obj, normalization.method = "LogNormalize", scale.factor = 10000)
seu.obj <- FindVariableFeatures(seu.obj, selection.method = "vst", nfeatures = 2000)
seu.obj <- ScaleData(seu.obj, features = rownames(seu.obj),verbose = TRUE) 
seu.obj <- RunPCA(seu.obj, features = VariableFeatures(object = seu.obj), npcs = 30, verbose = T)

saveRDS(seu.obj,file = 'data/processed_data/03.no_integrate_before_cluster.rds')

seu.obj <- FindNeighbors(seu.obj, dims = 1:30)
seu.obj <- FindClusters(seu.obj, resolution = seq(0.1,2,0.1))
seu.obj <- RunUMAP(seu.obj, dims = 1:30)

saveRDS(seu.obj,file = 'data/processed_data/04.no_integrate_after_cluster.rds')

# batch-------------------
(p_batch <- DimPlot(seu.obj,group.by = 'orig.ident',label = T,raster = F))
ggsave('results/figure/no_integrate_batch_umap.png',plot = p_batch,width = 12,height = 8)
ggsave('results/figure/no_integrate_batch_umap.pdf',plot = p_batch,width = 12,height = 8)

#  3. Cluster dimplot ---------------------------------------------------------
# — 3.1 by_resolution -----------------------------------------------------------
seu.obj <- readRDS('data/processed_data/04.no_integrate_after_cluster.rds')

for (i in seq(0.1,2,0.1)) {
  group<-paste("RNA_snn_res.",i,sep = "")
  title<-paste("resolution ",i,sep = "")
  p<-DimPlot(seu.obj, reduction = "umap", label = TRUE,group.by = group) + labs(title = title)
  assign(paste("p", i*10, sep=""),p)
}
p<-(p1|p2|p3|p4|p5)/(p6|p7|p8|p9|p10)
ggsave("results/figure/no_integrate_cluster_umap_0.1-1.png",p,width = 25,height = 8)
p<-(p11|p12|p13|p14|p15)/(p16|p17|p18|p19|p20)
ggsave("results/figure/no_integrate_cluster_umap_1.1-2.png",p,width = 25,height = 8)

