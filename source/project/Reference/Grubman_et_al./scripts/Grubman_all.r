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
# 终端运行时选择 multicore
library(future) # 新增
plan("multicore", workers = 30) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限

# 创建seurat对象-----------------
counts <- ReadMtx(
  mtx = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Grubman_et_al./data/raw_data/matrix.mtx.gz"),       # 对应 matrix 文件
  features = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Grubman_et_al./data/raw_data/features.tsv.gz"), # 对应基因/特征文件 (features)
  cells = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Grubman_et_al./data/raw_data/barcodes.tsv.gz"),     # 对应细胞条形码文件 (barcodes)
  feature.column = 1
)


# 创建 Seurat 对象
seu_obj <- CreateSeuratObject(
  counts = counts,
  project = "Grubman",       # 你可以给项目起个名字
  min.cells = 3,          # 可选：过滤掉少于3个细胞表达的基因
  min.features = 200      # 可选：过滤掉检测到少于200个基因的细胞
)

metadata <- read.table(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Grubman_et_al./data/raw_data/metadata.txt"), header = TRUE, row.names = 1, sep = "\t")

# 将元数据添加到 Seurat 对象中
seu_obj <- AddMetaData(object = seu_obj, metadata = metadata)
seu_obj$percent.mt <- seu_obj$pctMT


# cca 整合--------------------------------

data.filt <- seu_obj
data.filt.list <- SplitObject(data.filt, split.by = "batch")
data.filt.list <- lapply(X = data.filt.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = data.filt.list)

# perform integration
data.filt.anchors <- FindIntegrationAnchors(object.list = data.filt.list, anchor.features = features)
data.filt.combined <- IntegrateData(anchorset = data.filt.anchors)

DefaultAssay(data.filt.combined) <- "integrated"
data.filt.combined <- ScaleData(data.filt.combined,verbose = TRUE)
data.filt.combined <- RunPCA(data.filt.combined, npcs = 70, verbose = FALSE)
data.filt.combined <- RunUMAP(data.filt.combined, reduction = "pca", dims = 1:30)
data.filt.combined <- FindNeighbors(data.filt.combined, reduction = "pca", dims = 1:30)
data.filt.combined <- FindClusters(data.filt.combined, resolution = seq(0.1,1,0.1))


# redefine celltype----------
data.filt.combined$CellType<-as.character(data.filt.combined$integrated_snn_res.0.1)

data.filt.combined@meta.data %<>% 
  mutate(CellType=case_when(integrated_snn_res.0.1%in%c("0")~"Oligo",
                            integrated_snn_res.0.1%in%c("1")~"Astro",
                            integrated_snn_res.0.1%in%c("2")~"OPC",
                            integrated_snn_res.0.1%in%c("3")~"Unknown",
                            integrated_snn_res.0.1%in%c("4")~"Neuron",
                            integrated_snn_res.0.1%in%c("5")~"Microglia",
                            integrated_snn_res.0.1%in%c("6")~"Endo"))



qsave(data.filt.combined, file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Grubman_et_al./data/processed_data/Grubman.qs"))
