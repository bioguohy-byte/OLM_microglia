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
library(openxlsx)

# 多线程
library(future)
plan("multisession", workers = 6)
options(future.globals.maxSize = 100 * 1024^3)

# 加载数据--------
tmp_obj <- qread('data/processed_data/14.mg_combined_meta.qs')

# FindAllMarkers-----------------
Idents(tmp_obj) <- 'sub_clusters'
sub_clusters_markers <- FindAllMarkers(tmp_obj,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
qsave(sub_clusters_markers,'data/saved_data/sub_clusters_markers.qs')
write.xlsx(sub_clusters_markers,file = 'results/table/sub_clusters_markers.xlsx',overwrite = T)

# dotplot---------------
# 取各亚群top10的基因作图
sub_clusters_markers <- qread('data/saved_data/sub_clusters_markers.qs')

top10_sub_clusters_markers <- sub_clusters_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

# gene list
dotplot_gene_list <- top10_sub_clusters_markers$gene
dotplot_gene_list <- unique(dotplot_gene_list)


#' 自定义美化点图函数
#' @param obj Seurat对象
#' @param features 基因列表
#' @param group_var 分组的metadata列名（如 "sub_clusters"）
#' @param plot_title 图表标题
#' @param assay 使用的Assay，默认 'RNA'
#' @return 返回一个ggplot对象
plot_mg_dotplot <- function(obj, features, group_var, plot_title = "Gene Expression", assay = 'RNA') {
  
  # 1. 确保环境设置正确
  DefaultAssay(obj) <- assay
  
  # 2. 调用Seurat内置DotPlot获取基础数据
  # 使用 group.by 自动处理 idents，避免之前遇到的 factor level 报错
  p_base <- DotPlot(object = obj, 
                    features = unique(features), 
                    group.by = group_var)
  
  # 3. 提取数据并按你的风格进行 ggplot 绘制
  p_custom <- p_base$data %>%
    # 统一重命名列名，方便后续 aes 调用
    dplyr::rename(
      `Average Expression` = avg.exp.scaled, # 或者用 avg.exp，取决于你是否需要标准化
      `Percent Expressed` = pct.exp,
      Gene = features.plot,
      Cluster = id
    ) %>%
    ggplot(aes(x = Gene, y = Cluster, size = `Percent Expressed`, colour = `Average Expression`)) +
    geom_point() +
    # 使用你指定的配色方案
    scale_color_gradient2(low = "#000000", mid = "#007965", high = "#ffcc29", midpoint = 0) +
    theme_classic() +
    ggtitle(plot_title) +
    theme(
      axis.text = element_text(size = 11, color = "black"),
      plot.margin = margin(t = 4, r = 0, b = 1, l = 1, unit = "cm"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) + 
    RotatedAxis() # 调用Seurat内置的坐标轴旋转
  
  return(p_custom)
}


# 1. 准备基因列表
dotplot_gene_list <- top10_sub_clusters_markers$gene %>% unique()

# 2. 一行代码生成图片
p_final <- plot_mg_dotplot(
  obj = tmp_obj, 
  features = dotplot_gene_list, 
  group_var = "sub_clusters", 
  plot_title = "MG sub_clusters Gene Expression"
)

# 3. 展示并保存
p_final
ggsave("results/figure/MG_dotplot_v2.pdf", plot = p_final, width = 20, height = 5)
