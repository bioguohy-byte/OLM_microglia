# cca方法整合 ------------------------
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

library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)

library(future) # 新增
plan("multicore", workers = 4) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限


# data.filt <- readRDS("data/processed_data/08.no_integrate_after_qc_round2.rds")
# 
# data.filt.list <- SplitObject(data.filt, split.by = "sample")
# data.filt.list <- lapply(X = data.filt.list, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst",
#                             nfeatures = 2000)
# })
# 
# features <- SelectIntegrationFeatures(object.list = data.filt.list)
# 
# # perform integration
# data.filt.anchors <- FindIntegrationAnchors(object.list = data.filt.list, anchor.features = features)
# data.filt.combined <- IntegrateData(anchorset = data.filt.anchors)
# saveRDS(data.filt.combined,file = 'data/processed_data/09.integrate_cca_before_cluster.rds')
# 
# data.filt.combined <- readRDS('data/processed_data/09.integrate_cca_before_cluster.rds')
# 
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# 
# DefaultAssay(data.filt.combined) <- "RNA"
# data.filt.combined <- CellCycleScoring(data.filt.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# 
# DefaultAssay(data.filt.combined) <- "integrated"
# data.filt.combined <- ScaleData(data.filt.combined,
#                        vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
#                        verbose = TRUE)
# 
# saveRDS(data.filt.combined,file = 'data/processed_data/10.integrate_cca_after_cluster.rds')

data.filt.combined <- readRDS('data/processed_data/10.integrate_cca_after_cluster.rds')

data.filt.combined <- RunPCA(data.filt.combined, npcs = 70, verbose = FALSE)
data.filt.combined <- RunUMAP(data.filt.combined, reduction = "pca", dims = 1:50)
data.filt.combined <- FindNeighbors(data.filt.combined, reduction = "pca", dims = 1:50)
data.filt.combined <- FindClusters(data.filt.combined, resolution = seq(0.1,2,0.1))
saveRDS(data.filt.combined,file = 'data/processed_data/10.integrate_cca_after_cluster.rds')

DefaultAssay(data.filt.combined) <- "RNA"
p1 <- DimPlot(data.filt.combined,label = T,raster=F)
ggsave(Sys.getenv("MG_EXTERNAL_068"),p1,width = 6,height = 8)