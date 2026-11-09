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
library(harmony)
library(qs)

# 多线程运行------------------------
# # 终端运行时选择 multicore
# library(future) # 新增
# plan("multicore", workers = 4) # 根据你的服务器负载调整 workers 数量
# options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限

# rstudio server运行时选择 multisession
library(future)
plan(multisession, workers = 6)
nbrOfWorkers()
options(future.globals.maxSize = 1999999 * 1924^2)


# subset mg -------------------
all_data <- readRDS('data/processed_data/12.integrate_cca_CellType_meta_v1.rds')
all_data <- tmp_obj
tmp_obj <- subset(all_data,subset = CellType == "MG")

# harmony--------------------------------
tmp_obj <- NormalizeData(tmp_obj)
tmp_obj <- FindVariableFeatures(tmp_obj, selection.method = "vst", nfeatures = 2000)
tmp_obj <- ScaleData(tmp_obj,
                       vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"),
                       verbose = TRUE)
tmp_obj <- RunPCA(tmp_obj, npcs = 70, verbose = FALSE)

qsave(tmp_obj,'data/processed_data/12.5.mg_before_cluster.qs',nthreads = 8)
# tmp_obj <- qread('data/processed_data/12.5.mg_before_cluster.qs')

tmp_obj <- RunHarmony(tmp_obj, 
                        group.by.vars = "Sample", 
                        plot_convergence = TRUE, 
                        dims.use = 1:50,  # 可以先用较少的PC尝试
                        verbose = TRUE)
cat("开始UMAP降维和聚类...\n")
tmp_obj <- RunUMAP(tmp_obj, reduction = "harmony", dims = 1:30)
tmp_obj <- FindNeighbors(tmp_obj, reduction = "harmony", dims = 1:30)
tmp_obj <- FindClusters(tmp_obj, resolution = 0.1)
DimPlot(tmp_obj)

qsave(tmp_obj,'data/processed_data/13.mg_combined.qs',nthreads = 8)

#mg_combined meta -----------------------------------
mg_combined <- qread('data/processed_data/13.mg_combined.qs')

mg_combined$sub_clusters <- as.character(mg_combined$seurat_clusters)
mg_combined$sub_clusters[mg_combined$sub_clusters%in% c('1')]<-"M0"
mg_combined$sub_clusters[mg_combined$sub_clusters%in% c('4')]<-"M1"
mg_combined$sub_clusters[mg_combined$sub_clusters%in% c('2')]<-"M2"
mg_combined$sub_clusters[mg_combined$sub_clusters%in% c('5')]<-"M3"
mg_combined$sub_clusters[mg_combined$sub_clusters%in% c('0','3')]<-"M4"

sub_clusters_levels <- c("M0", "M1", "M2", "M3", "M4")
mg_combined$sub_clusters <- factor(mg_combined$sub_clusters, levels = sub_clusters_levels)
qsave(mg_combined,'data/processed_data/14.mg_combined_meta.qs',nthreads = 8)
