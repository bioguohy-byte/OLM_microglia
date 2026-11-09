setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/"))
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
plan("multicore", workers = 10) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限

# 9. 单细胞整合------------------
tryCatch({
  Su_mg <- qread('Su_et_al./data/processed_data/Su_mg.qs')
  Chia_mg <- qread('Chia_et_al./data/processed_data/Chia_mg.qs')
  Fullard_mg <- qread('Fullard_et_al./data/processed_data/Fullard_mg.qs')
  Jakel_mg <- qread('Jakel_et_al./data/processed_data/Jakel_mg.qs')
  Grubman_mg <- qread('Grubman_et_al./data/processed_data/Grubman_mg.qs')
  Yao_mg <- qread('Yao_et_al./data/processed_data/Yao_mg.qs')
  
}, error = function(e) { stop("Error loading files: ", e$message) })


DefaultAssay(Su_mg) <- 'RNA'
DefaultAssay(Chia_mg) <- 'RNA'
DefaultAssay(Fullard_mg) <- 'RNA'
DefaultAssay(Jakel_mg) <- 'RNA'
DefaultAssay(Grubman_mg) <- 'RNA'
DefaultAssay(Yao_mg) <- 'RNA'


data.filt <- merge(x=Su_mg,
                 y=c(Chia_mg,Fullard_mg,Jakel_mg,Grubman_mg,Yao_mg), 
                 add.cell.ids=c("Su_mg","Chia_mg","Fullard_mg","Jakel_mg","Grubman_mg",
                                "Yao_mg"))


tmp_obj <- data.filt %>% 
  NormalizeData() %>% 
  FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(., features = VariableFeatures(.), npcs = 30, verbose = T) %>% 
  FindNeighbors(.,dims = 1:30) %>% 
  FindClusters(.,resolution = seq(0.1,1,0.1)) %>%
  RunUMAP(., seed.use = 42, reduction = "pca", dims = 1:30)


DimPlot(tmp_obj,group.by = 'Study')

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

Batch_name <- "Study"
tmp_obj.list <- SplitObject(tmp_obj, split.by = Batch_name)
tmp_obj.list <- lapply(X = tmp_obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = tmp_obj.list,nfeatures = 2000)
seurat.anchors <- FindIntegrationAnchors(object.list = tmp_obj.list,anchor.features = features)
tmp_Sobj <- IntegrateData(anchorset = seurat.anchors)
DefaultAssay(tmp_Sobj) <- "RNA"
tmp_Sobj <- ScaleData(tmp_Sobj,verbose = T)
tmp_Sobj <- CellCycleScoring(tmp_Sobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(tmp_Sobj) <- "integrated"
tmp_Sobj <- ScaleData(tmp_Sobj,verbose = T)
tmp_Sobj_pca <- RunPCA(tmp_Sobj, npcs = 30, verbose = T)
qsave(tmp_Sobj_pca,"combine/data/processed_data/tmp_Sobj.after_PCA.v1.qs")

tmp_Sobj <- RunUMAP(tmp_Sobj_pca, seed.use = 42, reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(., reduction = "pca") %>% 
  FindClusters(., resolution = seq(0.1,1,0.1))
qsave(tmp_Sobj,"combine/data/processed_data/tmp_Sobj.after_FindClusters.v1.qs")

tmp_Plist<-list()
for (i in grep("integrated_snn_res",colnames(tmp_Sobj@meta.data),value = T)[1:10]) {
  title<-gsub(".*res.","",i) %>% paste("Resolution ",.,sep = "")
  p<-DimPlot(tmp_Sobj, reduction = "umap", label = TRUE,group.by = i,raster = T) + labs(title = title) + NoLegend()
  tmp_Plist[[as.character(i)]]<-p
}
(p1<-Reduce("+",tmp_Plist)+plot_layout(ncol = 5))
ggsave(paste("results/figure/tmp_Sobj.Resolution.Umap.1.v1.pdf",sep = ""),p1,width = 18,height = 9)

tmp_Plist<-list()
for (i in grep("integrated_snn_res",colnames(tmp_Sobj@meta.data),value = T)[11:20]) {
  title<-gsub(".*res.","",i) %>% paste("Resolution ",.,sep = "")
  p<-DimPlot(tmp_Sobj, reduction = "umap", label = TRUE,group.by = i,raster = T) + labs(title = title) + NoLegend()
  tmp_Plist[[as.character(i)]]<-p
}
(p2<-Reduce("+",tmp_Plist)+plot_layout(ncol = 5))
ggsave(paste("results/figure/tmp_Sobj.Resolution.Umap.2.v1.pdf",sep = ""),p2,width = 18,height = 9)
# 10. 细胞类型定义------------------