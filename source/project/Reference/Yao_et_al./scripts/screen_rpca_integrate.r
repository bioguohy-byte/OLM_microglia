setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Yao_et_al./"))
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
library(DoubletFinder)

# 多线程运行
# # 终端运行时选择 multicore
library(future) # 新增
plan("multicore", workers = 20) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限

# 9. 单细胞整合------------------
# — 9.2 RPCA 整合-------------------
data.filt <- qread("data/processed_data/Data_after_filtered.v1.qs")

# 准备 Meta data 和创建 Seurat 对象
tmp_meta <- data.filt@meta.data %>% 
  .[,setdiff(colnames(.), grep("RNA_snn_res.", colnames(.), value = T))] %>% 
  dplyr::select(., -seurat_clusters)

tmp_obj <- data.filt %>% 
  {.@assays[["RNA"]]@counts} %>% 
  CreateSeuratObject(., min.features = 200, min.cell = 3, meta.data = tmp_meta)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
Batch_name <- "orig.ident"

# 拆分对象
tmp_obj.list <- SplitObject(tmp_obj, split.by = Batch_name)

# 标准化和寻找高变基因
tmp_obj.list <- lapply(X = tmp_obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# 选择整合特征
features <- SelectIntegrationFeatures(object.list = tmp_obj.list, nfeatures = 2000)

# RPCA 核心步骤：整合前需对每个对象单独 Scale 和 PCA
tmp_obj.list <- lapply(X = tmp_obj.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# 使用 RPCA 寻找锚点
seurat.anchors <- FindIntegrationAnchors(object.list = tmp_obj.list, 
                                         anchor.features = features, 
                                         reduction = "rpca")

tmp_Sobj <- IntegrateData(anchorset = seurat.anchors)

# 切换回 RNA assay 进行细胞周期评分
DefaultAssay(tmp_Sobj) <- "RNA"
tmp_Sobj <- ScaleData(tmp_Sobj, verbose = T)
tmp_Sobj <- CellCycleScoring(tmp_Sobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# 切换到 integrated assay 进行后续降维
DefaultAssay(tmp_Sobj) <- "integrated"
tmp_Sobj <- ScaleData(tmp_Sobj, verbose = T)
tmp_Sobj_pca <- RunPCA(tmp_Sobj, npcs = 100, verbose = T)
qsave(tmp_Sobj_pca, "data/processed_data/tmp_Sobj.rpca.after_PCA.v1.qs")

# UMAP 和 聚类
tmp_Sobj <- RunUMAP(tmp_Sobj_pca, seed.use = 42, reduction = "pca", dims = 1:70) %>% 
  FindNeighbors(., reduction = "pca",dims = 1:70) %>% 
  FindClusters(., resolution = seq(0.1, 2, 0.1))
qsave(tmp_Sobj, "data/processed_data/tmp_Sobj.rpca.after_FindClusters.v1.qs")

# 绘图部分 1
tmp_Plist <- list()
for (i in grep("integrated_snn_res", colnames(tmp_Sobj@meta.data), value = T)[1:10]) {
  title <- gsub(".*res.", "", i) %>% paste("Resolution ", ., sep = "")
  p <- DimPlot(tmp_Sobj, reduction = "umap", label = TRUE, group.by = i, raster = T) + 
    labs(title = title) + NoLegend()
  tmp_Plist[[as.character(i)]] <- p
}
(p1 <- Reduce("+", tmp_Plist) + plot_layout(ncol = 5))
ggsave("results/figure/tmp_Sobj.rpca.Resolution.Umap.1.v1.pdf", p1, width = 18, height = 9)

# 绘图部分 2
tmp_Plist <- list()
for (i in grep("integrated_snn_res", colnames(tmp_Sobj@meta.data), value = T)[11:20]) {
  title <- gsub(".*res.", "", i) %>% paste("Resolution ", ., sep = "")
  p <- DimPlot(tmp_Sobj, reduction = "umap", label = TRUE, group.by = i, raster = T) + 
    labs(title = title) + NoLegend()
  tmp_Plist[[as.character(i)]] <- p
}
(p2 <- Reduce("+", tmp_Plist) + plot_layout(ncol = 5))
ggsave("results/figure/tmp_Sobj.rpca.Resolution.Umap.2.v1.pdf", p2, width = 18, height = 9)

