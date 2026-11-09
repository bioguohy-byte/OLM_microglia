setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/"))
rm(list=ls())

# 0. Libraries -----------------------------------------------------------------
library(Seurat)
library(harmony)
library(qs)
library(dplyr)
library(patchwork)
library(ggplot2)

library(future) # 新增
plan("multicore", workers = 20) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限

# 1. Load Data & Prepare Metadata ----------------------------------------------
message("\n[Step 1] Loading Datasets...")

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

list <- list(Su_mg,Chia_mg,Fullard_mg,Jakel_mg,Grubman_mg,Yao_mg)

# 2. SCT + CCA----------------------
list <- lapply(X = list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 2000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",anchor.features = features)
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
qsave(combined.sct,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_before_inte.qs"))

combined.sct <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_before_inte.qs"))

DefaultAssay(combined.sct)<- "integrated"
combined.sct <- RunPCA(combined.sct, verbose = FALSE)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:10)
combined.sct <- FindNeighbors(combined.sct,  reduction = "pca",dims = 1:10)
combined.sct <- FindClusters(combined.sct, resolution = seq(0.1, 1, 0.1))
# combined.sct <- FindClusters(combined.sct, resolution = 0.25)

DefaultAssay(combined.sct)<-"RNA"
combined.sct <- NormalizeData(combined.sct)

qsave(combined.sct,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca.qs"))

