# harmony整合
# harmony方法整合 ------------------------
setwd(Sys.getenv("MG_PROJECT_ROOT"))
rm(list=ls())

# 0. library --------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Matrix)
library(cowplot)
library(reshape2)
library(SingleCellExperiment)
library(umap)
library(ggthemes)
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
# 新增 harmony 包
library(harmony)

library(future) # 新增
plan("multicore", workers = 5) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限

data.filt <- readRDS(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/merged_data.rds"))

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
data.filt <- CellCycleScoring(data.filt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

data.filt <- NormalizeData(data.filt)
data.filt <- FindVariableFeatures(data.filt, selection.method = "vst", nfeatures = 2000)
data.filt <- ScaleData(data.filt, vars.to.regress = c("nCount_RNA","percent.mt","S.Score","G2M.Score"),verbose = T)
data.filt <- RunPCA(data.filt, npcs = 70, verbose = FALSE)
data.filt <- RunHarmony(data.filt, 
                        group.by.vars = "Batch", 
                        plot_convergence = TRUE, 
                        dims.use = 1:50,  # 可以先用较少的PC尝试
                        verbose = TRUE)

saveRDS(data.filt,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/harmony_before_cluster.rds"))

cat("开始UMAP降维和聚类...\n")
data.filt <- RunUMAP(data.filt, reduction = "harmony", dims = 1:50)
data.filt <- FindNeighbors(data.filt, reduction = "harmony", dims = 1:50)
data.filt <- FindClusters(data.filt, resolution = seq(0.1,1,0.1))

cat("保存最终聚类结果...\n")
saveRDS(data.filt, file = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/harmony_after_cluster.rds"))

# 可视化------------------------------------------------------------------------

# 加载数据用于可视化
# load('data/processed_data/10.integrate_harmony_after_cluster.RData')

# 绘制批次分布图，检查批次效应是否被有效去除
cat("生成批次效应评估图...\n")
batch_plot <- DimPlot(data.filt, reduction = "umap", group.by = "Batch", 
                      label = FALSE, pt.size = 0.5) +
  ggtitle("Batch Distribution after Harmony Integration") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

celltype_plot <- DimPlot(data.filt, reduction = "umap", group.by = "seurat_clusters", 
                         label = TRUE, pt.size = 0.5) +
  ggtitle("Cell Clusters after Harmony Integration") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# 并排显示批次和聚类图
combined_plot <- batch_plot | celltype_plot
ggsave("results/figure/harmony_integration_evaluation.png", combined_plot, 
       width = 12, height = 5, dpi = 300)

# 根据resolution画dimplot（取消注释的代码）
cat("生成不同分辨率的聚类图...\n")
for (i in seq(0.1,1,0.1)) {
  group <- paste("RNA_snn_res.", i, sep = "") # 注意：如果没有integrated assay，这里通常是 RNA_snn_res
  title <- paste("resolution ", i, sep = "")
  p <- DimPlot(data.filt, reduction = "umap", label = TRUE, group.by = group) + 
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  assign(paste("p", i*10, sep=""), p)
}

# 分批保存图片
p1 <- (p1|p2|p3|p4|p5)/(p6|p7|p8|p9|p10)
ggsave("results/figure/50dims_harmony_cluster1_umap_0.1-1.png", p1, width = 25, height = 8, dpi = 300)


cat("Harmony整合流程完成！\n")
cat("整合后的数据维度:", dim(data.filt), "\n")
cat("识别的细胞类型数量:", length(unique(data.filt$seurat_clusters)), "\n")