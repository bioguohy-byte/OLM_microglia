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
plan("multicore", workers = 50) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限

# # subset mg------------------
# Su_v2 <- qread('Su_et_al./data/processed_data/Su_CSC_alldata_v2.qs')
# Su_mg <- subset(Su_v2,subset = CellType %in% c('MG'))
# qsave(Su_mg,'Su_et_al./data/processed_data/Su_mg.qs')
# 
# Chia_v2 <- qread('Chia_et_al./data/processed_data/Chia_v2.qs')
# Chia_mg <- subset(Chia_v2,subset = CellType %in% c('Microglia'))
# qsave(Chia_mg,'Chia_et_al./data/processed_data/Chia_mg.qs')
# 
# Fullard_v2 <- qread('Fullard_et_al./data/processed_data/Fullard_v2.qs')
# Fullard_mg <- subset(Fullard_v2,subset = CellType %in% c('Mic'))
# qsave(Fullard_mg,'Fullard_et_al./data/processed_data/Fullard_mg.qs')
# 
# Jakel_v2 <- qread('Jakel_et_al./data/processed_data/Jakel_v2.qs')
# Jakel_mg <- subset(Jakel_v2,subset = CellType %in% c('Immune_cells','ImOlGs','Microglia_Macrophages'))
# qsave(Jakel_mg,'Jakel_et_al./data/processed_data/Jakel_mg.qs')
# 
# Grubman_v2 <- qread('Grubman_et_al./data/processed_data/Grubman_v2.qs')
# Grubman_mg <- subset(Grubman_v2,subset = CellType %in% c('Microglia'))
# qsave(Grubman_mg,'Grubman_et_al./data/processed_data/Grubman_mg.qs')
# 
# Yao_v2 <- qread('Yao_et_al./data/processed_data/Yao_v2.qs')
# Yao_mg <- subset(Yao_v2,subset = CellType %in% c('MG'))
# qsave(Yao_mg,'Yao_et_al./data/processed_data/Yao_mg.qs')


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
# DefaultAssay(Jakel_mg) <- 'RNA'
DefaultAssay(Grubman_mg) <- 'RNA'
DefaultAssay(Yao_mg) <- 'RNA'
# 确保所有对象放入 list
list <- list(Su_mg, Chia_mg, Fullard_mg, Jakel_mg, Grubman_mg, Yao_mg)

# ------------------------------------------------------------------------------
# 2. SCTransform + RPCA Integration
# ------------------------------------------------------------------------------
message("\n[Step 2] Running SCTransform and RPCA Integration...")

# 1. 对每个对象运行 SCTransform
# 注意：如果数据量非常大，建议在该步骤增加 vst.flavor="v2" 或 method="glmGamPoi" (如果安装了的话) 以提速
list <- lapply(X = list, FUN = SCTransform, method = "glmGamPoi", verbose = FALSE) 

# 2. 选择整合特征 (Integration Features)
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000) # 建议 RPCA 稍微增加特征数到 3000

# 3. 准备 SCT 整合 (确保所有对象都有这些特征的计算残差)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)

# 4. [关键步骤] 对每个对象运行 PCA
# RPCA 需要预先计算 PCA，且必须基于选定的 integration features
list <- lapply(X = list, FUN = RunPCA, features = features, verbose = FALSE)

# 5. 寻找锚点 (使用 RPCA)
# reduction = "rpca" 是核心参数
# k.anchor 默认是 5，如果数据集批次效应差异巨大，可以适当调整，但通常默认即可
anchors <- FindIntegrationAnchors(
  object.list = list, 
  normalization.method = "SCT", 
  anchor.features = features, 
  dims = 1:30, 
  reduction = "rpca", 
  k.anchor = 20 # RPCA 通常可以处理更多锚点，有助于提高稳健性，Seurat v4/v5 推荐在 RPCA 时适当通过 k.anchor 调整
)

# 6. 整合数据
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)

# 保存中间结果
qsave(combined.sct, file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/after_inte_rpca.qs"))

# ------------------------------------------------------------------------------
# 3. Downstream Analysis
# ------------------------------------------------------------------------------
message("\n[Step 3] Downstream Processing (PCA, UMAP, Clustering)...")

DefaultAssay(combined.sct) <- "integrated"

# 整合后的数据再次运行 PCA (这是在整合空间上运行)
combined.sct <- RunPCA(combined.sct, verbose = FALSE)

# Run UMAP
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30)

# Find Neighbors & Clusters
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:30)
combined.sct <- FindClusters(combined.sct, resolution = seq(0.1, 1, 0.1))

# 切换回 RNA assay 进行标准化 (用于可视化/差异分析)
DefaultAssay(combined.sct) <- "RNA"
combined.sct <- NormalizeData(combined.sct)

# Save Final Object
qsave(combined.sct, file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_rpca.qs"))

message("Done!")