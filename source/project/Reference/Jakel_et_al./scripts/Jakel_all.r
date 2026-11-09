setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Jakel_et_al."))
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

# 多线程运行------------------------
# 终端运行时选择 multicore
library(future) # 新增
plan("multicore", workers = 50) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限

# seurat_obj <- readRDS(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Jakel_et_al./data/processed_data/Jakel.rds"))
# 
# # nCount nFeature
# counts <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")
# seurat_obj$nCount_RNA <- colSums(counts)
# seurat_obj$nFeature_RNA <- colSums(counts > 0)
# 
# # percent.mt
# seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
# seurat_obj$orig.ident <- seurat_obj$sample
# 
# # 7. DoubletFinder------------------------------------
# tmp_obj <- seurat_obj
# 
# raw_type <- "Celltypes"
# tmp_Idents<-"RNA_snn_res.0.2"
# tmp_meta<-tmp_obj@meta.data
# all_DL<-data.frame()
# summary<-data.frame()
# tmp_Plist<-list()
# 
# for (i in unique(tmp_obj$orig.ident)) {
#   # i<-"NTES"
#   subobj<-subset(tmp_obj, subset = orig.ident %in% i)
#   DoubletRate = ncol(subobj)*8*1e-6-0.002
#   # DoubletRate = 0.04
#   
#   ## Pre-process Seurat object (standard)
#   seu_obj <- NormalizeData(subobj)
#   seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
#   seu_obj <- ScaleData(seu_obj)
#   seu_obj <- RunPCA(seu_obj)
#   seu_obj <- RunUMAP(seu_obj, dims = 1:10)
#   
#   ## pK Identification (no ground-truth)
#   sweep.res.list <- paramSweep_v3(seu_obj, PCs = 1:10, sct = FALSE)
#   sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
#   bcmvn <- find.pK(sweep.stats)
#   pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
#   
#   ## Homotypic Doublet Proportion Estimate
#   annotations <- tmp_meta[tmp_meta$orig.ident%in%i,][[raw_type]]
#   homotypic.prop <- modelHomotypic(annotations)
#   nExp_poi <- round(DoubletRate*nrow(seu_obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
#   nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#   
#   ## Run DoubletFinder with varying classification stringencies
#   seu_obj <- doubletFinder_v3(seu_obj, PCs = 1:10, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
#   
#   p<-DimPlot(seu_obj, reduction = "umap", group.by = grep("^DF",colnames(seu_obj@meta.data),value = T))+labs(title = i)
#   tmp_Plist[[i]]<-p
#   
#   DL<-seu_obj@meta.data[,c(grep("^DF",colnames(seu_obj@meta.data),"orig.ident"),value = T),drop=F] %>%  
#     {colnames(.)<-c("if_DL","orig.ident");.}
#   
#   all_DL<-rbind(all_DL,DL)
#   summary<-rbind(summary,data.frame(sample=i,nCells=nrow(seu_obj@meta.data),
#                                     Raw_Doublet_rate=DoubletRate,Raw_Doublet=nExp_poi,
#                                     Homotypic_rate=homotypic.prop,Hetertypic_Doublet=nExp_poi.adj))
#   
# }
# qsavem(all_DL, summary, tmp_Plist, file = "data/processed_data/DoubletFinder.result.qs")
# 
# Reduce("+",tmp_Plist)+plot_layout(ncol = 4,guides = "collect")
# 
# tmp_obj@meta.data %<>% {.$DoubletFinder<-"Singlet";.} %<>% 
#   {.[rownames(all_DL[all_DL$if_DL%in%"Doublet",]),]$DoubletFinder<-"Doublet";.}
# 
# qsave(tmp_obj,"data/processed_data/Raw_data.DoubletFinder.qs")

# 8. 第一轮质控--------------------

tmp_obj <- qread('data/processed_data/Raw_data.DoubletFinder.qs')
alldata<-tmp_obj

# step1
data.filt <- alldata %>% subset(.,subset=DoubletFinder%in%"Singlet")

# 9. 单细胞整合------------------

tmp_meta<-data.filt@meta.data %>% 
  .[,setdiff(colnames(.),grep("RNA_snn_res.",colnames(.),value = T))]

tmp_obj<-data.filt %>% 
  {.@assays[["RNA"]]@counts} %>% 
  CreateSeuratObject(.,min.features=200,min.cell=3,meta.data = tmp_meta)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

Batch_name <- "donor_ID"
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
qsave(tmp_Sobj_pca,"data/processed_data/tmp_Sobj.after_PCA.v1.qs")

tmp_Sobj <- RunUMAP(tmp_Sobj_pca, seed.use = 42, reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(., reduction = "pca") %>% 
  FindClusters(., resolution = seq(0.1,2,0.1))
qsave(tmp_Sobj,"data/processed_data/tmp_Sobj.after_FindClusters.v1.qs")
qsave(tmp_Sobj,"data/processed_data/Jakel.qs")
