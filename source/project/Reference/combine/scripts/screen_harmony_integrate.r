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
plan("multicore", workers = 50) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限

# 9. 单细胞整合------------------
# — 9.3 Harmony 整合--------------
# ==============================================================================
# 0. 数据读取与预处理 (Data Loading & Preprocessing)
# ==============================================================================

# 1. 读取数据
Su_v2      <- qread('Su_et_al./data/processed_data/Su_CSC_alldata_v2.qs')
Yao_v2     <- qread('Yao_et_al./data/processed_data/Yao_v2.qs')

# 2. 构建列表并处理元数据 (Metadata Handling)
# 将对象放入列表方便循环处理
obj_list <- list(
  Su = Su_v2, 
  Yao = Yao_v2
)

# 循环处理：添加 Study 标识，并构建全局唯一的 Batch ID
# 假设原始 Batch 信息存储在 meta.data 的 "Batch" 列中
for (study_name in names(obj_list)) {
  # 1. 标记来源于哪个 Study
  obj_list[[study_name]]$Study <- study_name
  
  # 2. 构建唯一的整合 Batch ID (Study名_原始Batch名)
  # 这样可以防止不同文章都有 "Batch1" 导致的混淆
  obj_list[[study_name]]$Unique_Batch <- paste(study_name, obj_list[[study_name]]$Batch, sep = "_")
}

# 3. 合并为一个大对象
# add.cell.ids 非常重要，防止不同数据集的细胞 Barcode 重复
tmp_obj <- merge(
  x = obj_list[[1]], 
  y = obj_list[2:length(obj_list)], 
  add.cell.ids = names(obj_list), 
  project = "Human_Brain_Integration"
)

# 4. 清理内存 (可选，防止内存溢出)
rm(Su_v2, Yao_v2, obj_list)
gc()

# ==============================================================================
# 下面开始套用您的流程 (Pipeline)
# ==============================================================================

# 设置细胞周期基因
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# *** 关键设置 *** # 这里指定我们刚刚生成的唯一 ID 作为整合的批次变量
Batch_name <- "Unique_Batch" 

# 1. 筛选高质量整合特征
# 注意：SplitObject 这里会根据 Unique_Batch 把大对象拆开，分别找高变基因
tmp_obj.list <- SplitObject(tmp_obj, split.by = Batch_name)

tmp_obj.list <- lapply(X = tmp_obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = tmp_obj.list, nfeatures = 2000)

# 2. 数据处理与 PCA
# 将整合特征赋给主对象
VariableFeatures(tmp_obj) <- features
tmp_obj <- NormalizeData(tmp_obj)
tmp_obj <- ScaleData(tmp_obj, verbose = T)

# 细胞周期评分
tmp_obj <- CellCycleScoring(tmp_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# 再次 Scale (基于评分后的数据) 并运行 PCA
# 如果想去除细胞周期影响，可以在这里加 vars.to.regress = c("S.Score", "G2M.Score")
tmp_obj <- ScaleData(tmp_obj, verbose = T)
tmp_Sobj_pca <- RunPCA(tmp_obj, npcs = 50, verbose = T)

# 3. 运行 Harmony
# group.by.vars 使用我们定义的 Unique_Batch (即 Study + Sample)
tmp_Sobj_pca <- RunHarmony(tmp_Sobj_pca, group.by.vars = Batch_name)

# 保存 PCA 结果
qsave(tmp_Sobj_pca, "combine/data/processed_data/tmp_Sobj.harmony.after_PCA.v1.qs")

# 4. UMAP 和 聚类 (使用 harmony reduction)
tmp_Sobj <- RunUMAP(tmp_Sobj_pca, reduction = "harmony", dims = 1:50, seed.use = 42) %>% 
  FindNeighbors(., reduction = "harmony", dims = 1:50) %>% 
  FindClusters(., resolution = seq(0.1, 2, 0.1))

# 保存聚类结果
qsave(tmp_Sobj, "combine/data/processed_data/tmp_Sobj.harmony.after_FindClusters.v1.qs")

# 5. 绘图
cluster_cols <- grep("_snn_res", colnames(tmp_Sobj@meta.data), value = T)

# 绘图部分 1
if(length(cluster_cols) > 0) {
  tmp_Plist <- list()
  # 限制循环次数，防止只有不到10个resolution时报错
  loop_end <- min(10, length(cluster_cols))
  
  for (i in cluster_cols[1:loop_end]) {
    title <- gsub(".*res.", "", i) %>% paste("Resolution ", ., sep = "")
    p <- DimPlot(tmp_Sobj, reduction = "umap", label = TRUE, group.by = i, raster = T) + 
      labs(title = title) + NoLegend()
    tmp_Plist[[as.character(i)]] <- p
  }
  (p1 <- Reduce("+", tmp_Plist) + plot_layout(ncol = 5))
  ggsave("combine/results/figure/tmp_Sobj.harmony.Resolution.Umap.1.v1.pdf", p1, width = 18, height = 9)
}

# 绘图部分 2
if(length(cluster_cols) > 10) {
  tmp_Plist <- list()
  loop_end <- min(20, length(cluster_cols))
  
  for (i in cluster_cols[11:loop_end]) {
    title <- gsub(".*res.", "", i) %>% paste("Resolution ", ., sep = "")
    p <- DimPlot(tmp_Sobj, reduction = "umap", label = TRUE, group.by = i, raster = T) + 
      labs(title = title) + NoLegend()
    tmp_Plist[[as.character(i)]] <- p
  }
  (p2 <- Reduce("+", tmp_Plist) + plot_layout(ncol = 5))
  ggsave("combine/results/figure/tmp_Sobj.harmony.Resolution.Umap.2.v1.pdf", p2, width = 18, height = 9)
}