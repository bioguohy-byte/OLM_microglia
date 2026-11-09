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
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)

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

# 开始作图
# 常规作图
DefaultAssay(tmp_obj) <- 'RNA'
p<-DotPlot(object = tmp_obj,  #指定使用哪个seurat object
           features = dotplot_gene_list,  #指定使用哪些基因
           idents = tmp_obj$sub_clusters) + RotatedAxis() #横坐标倾斜

# 利用数据自己绘制
(p2<-p$data %>% 
    {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.} %>%
    ggplot(.,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
    geom_point()+
    scale_color_gradient2(low="#000000",mid="#007965",high ="#ffcc29",midpoint = 0)+
    theme_classic()+
    ggtitle("MG sub_clusters Gene Expression") +  # 添加图形标题
    theme(axis.text=element_text(size=11, color="black"),
          plot.margin = margin(t = 4,  # Top margin
                               r = 0,  # Right margin
                               b = 1,  # Bottom margin
                               l = 1,  # Left margin
                               unit = "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16)  # 设置标题居中并增加字体大小
    ) + RotatedAxis())
p2

# 保存图片
ggsave("results/figure/MG_sub_clusters_top10_dotplot.png", plot = p2,width = 18, height = 5)
ggsave("results/figure/MG_sub_clusters_top10_dotplot.pdf", plot = p2,width = 18, height = 5)

# heatmap--------------------
# 取top50的基因
top50_sub_clusters_markers <- sub_clusters_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 50) %>%
  ungroup()


# Plot the top 50 markers using mean expression within each group.
Exp_av <- AverageExpression(tmp_obj,
                            group.by = "sub_clusters",
                            assays = "RNA")

Exp_av <- Exp_av[[1]]
Exp_av <- as.data.frame(Exp_av)
Exp_av$gene <- rownames(Exp_av)
result <- Exp_av %>%
  semi_join(top50_sub_clusters_markers, by = "gene")

#按照top50$gene的顺序来排列result
result <- result[match(top50_sub_clusters_markers$gene, result$gene),]

#现在的result中有一列是gene列，所以不能直接放在pheatmap中进行绘图，#因为里面不全是数字了，所以要删除那一列
result <- subset(result, select = -gene)

pheatmap(result)

p <- pheatmap(result, cluster_cols = FALSE, cluster_rows = FALSE,scale = "row",show_rownames = F,legend = FALSE)
p
# cluster_cols = FALSE: 这个参数决定是否对列进行聚类。如果设置为FALSE，那么列就不会被聚类。
# cluster_rows = FALSE: 这个参数决定是否对行进行聚类。如果设置为FALSE，那么行就不会被聚类。
# scale = "row": 这个参数决定是否对数据进行缩放，以及如何进行缩放。如果设置为"row"，那么每一行的数据都会被缩放，使得每一行的平均值为0，标准差为1。这样可以帮助你更好地比较不同基因的表达模式。


# 保存图片
ggsave("results/figure/MG_sub_clusters_top50_heatmap.png", plot = p,width = 5, height = 18)
ggsave("results/figure/MG_sub_clusters_top50_heatmap.pdf", plot = p,width = 5, height = 18)

# GO----------------------
sub_clusters_markers <- qread('data/saved_data/sub_clusters_markers.qs')

filt_0.5_sub_clusters_markers <- sub_clusters_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  ungroup()

#写循环，分别对每个cluster中的基因做GO分析，放到一个list当中
object <- filt_0.5_sub_clusters_markers

res_go <- list()
cc <- levels(object$cluster)
for (i in cc) {
  g1 <- object$gene[which(object$cluster == i)]
  print(g1)
  gg <- bitr(g1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") #记得对应人还是鼠
  ego <- enrichGO(
    gene = gg$ENTREZID,
    OrgDb = org.Hs.eg.db, #记得对应人还是鼠
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  res_go[[i]] <- ego
}

#保存结果
qsave(res_go, file = "data/saved_data/filt_0.5_sub_clusters_markers_res_go.qs")
res_go <- qread("data/saved_data/filt_0.5_sub_clusters_markers_res_go.qs")

(p1 <- merge_result(res_go) %>%
    dotplot(., showCategory = 5,includeAll = FALSE))

ggsave('results/figure/res_go_5_clusters_dotplot.png',p1,width = 7,height = 12,bg='#ffffff')
ggsave('results/figure/res_go_5_clusters_dotplot.pdf',p1,width = 7,height = 12,bg='#ffffff')

#把list中的结果合并在一起，保存为表格
merge_res_go <- merge_result(res_go)
qsave(merge_res_go, "data/saved_data/filt_0.5_sub_clusters_markers_merge_res_go.qs")
write.xlsx(merge_res_go,file = 'results/table/filt_0.5_sub_clusters_markers_merge_res_go.xlsx',overwrite = T)

# 筛选出的各个cluster的GO分别作图
# 读取数据
merge_res_go <- qread("data/saved_data/filt_0.5_sub_clusters_markers_merge_res_go.qs")
filter_GO <- as.data.frame(merge_res_go)

# 2. 定义 5 个 Cluster 的专属配色
# M0-M4 对应不同的颜色，方便在论文中区分
cluster_colors <- c(
  "M0" = "#E64B35", # 红色
  "M1" = "#4DBBD5", # 蓝色
  "M2" = "#00A087", # 绿色
  "M3" = "#3C5488", # 深蓝
  "M4" = "#F39B7F"  # 橙色
)

# 3. 循环生成并保存每个 Cluster 的 Top 10 柱状图
for (i in names(cluster_colors)) {
  
  # 筛选当前 Cluster 的数据并取 Top 10 (按 p.adjust 排序)
  plot_data <- filter_GO %>%
    filter(Cluster == i) %>%
    arrange(p.adjust) %>%  # p.adjust 越小越显著
    head(10) %>%
    mutate(logP = -log10(p.adjust)) # 计算用于绘图的纵坐标
  
  # 如果该簇没有结果则跳过
  if(nrow(plot_data) == 0) next
  
  # 设置因子水平，确保绘图顺序（最显著的在最上方）
  plot_data$Description <- factor(plot_data$Description, levels = rev(plot_data$Description))
  
  # 创建融合风格的图表
  p <- ggplot(plot_data, aes(x = logP, y = Description)) +
    # 动态映射该 Cluster 的颜色
    geom_bar(stat = "identity", width = 0.8, alpha = 0.7, fill = cluster_colors[i]) +
    # 融合 p1 的文字风格：将 GO 条目名称直接写在条形图起始位置
    geom_text(aes(x = 0.05, label = Description), size = 4.5, hjust = 0, color = "black") +
    theme_minimal() +
    labs(
      title = paste0("Top 10 GO Terms (BP): Cluster ", i),
      x = "-log10(p.adjust)",
      y = ""
    ) +
    theme(
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"), 
      axis.text.y = element_blank(), # 隐藏 Y 轴刻度文字（因为文字已在图中）
      axis.ticks.x = element_line(color = "black"),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 12)
    ) +
    # 留出 expansion 空间，防止文字标签被切掉
    scale_x_continuous(expand = c(0, 0, 0.2, 0)) 
  
  # 4. 保存为矢量 PDF
  # 确保目录存在
  dir.create("results/figure/cluster_go/", recursive = TRUE, showWarnings = FALSE)
  
  file_name <- paste0('results/figure/cluster_go/Cluster_', i, '_top10_GO.pdf')
  ggsave(file_name, plot = p, width = 8, height = 5)
  
  # 同时保存一份预览用的 PNG
  ggsave(gsub(".pdf", ".png", file_name), plot = p, width = 8, height = 5, dpi = 300)
}

cat("所有亚群的 GO Top 10 柱状图已绘制完成，并已按颜色分类保存。\n")

# pySCENIC----------------
# 新版的流程和脚本全都存储在${MG_PROJECT_ROOT}/scripts/SCENIC 这个文件夹下面

# 小胶质亚群4个基因集打分 AddModuleScore --------------
tmp_obj <- qread('data/processed_data/14.mg_combined_meta.qs')

# 假设你的对象名为 tmp_obj
# 定义基因列表（名称与你的 sub_clusters 一一对应）
mg_marker_list <- list(
  M0 = c('IL1RAPL1','RNF220', 'ST18', 'PCDH9', 'SLC24A2'),
  M1 = c('DPP10', 'CTNNA2', 'RORA', 'SORBS1', 'GPC5'),
  M2 = c('RBFOX1', 'CSMD1', 'FAM155A', 'DLGAP1', 'OPCML'),
  M3 = c('F13A1', 'CD163', 'VAV3', 'DAPK1', 'MCTP1'),
  M4 = c('CD83', 'FOXP2', 'FAM149A', 'LPCAT2', 'FRMD4A')
)

# 2. 自动化 AddModuleScore
# 循环打分，避免重复写五遍
for (id in names(mg_marker_list)) {
  tmp_obj <- AddModuleScore(
    object = tmp_obj,
    features = list(mg_marker_list[[id]]),
    name = paste0(id, "_Score")
  )
}

# Seurat 默认会在名字后加 "1"，记录下这些特征名
score_features <- paste0(names(mg_marker_list), "_Score1")

# 3. 自动化批量生成 FeaturePlot
# 设定配色
pal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(10)

# 使用 lapply 批量生成图表
plots <- lapply(names(mg_marker_list), function(id) {
  feature_name <- paste0(id, "_Score1")
  
  FeaturePlot(tmp_obj, features = feature_name, order = TRUE, cols = pal, raster = FALSE) +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "none", # 如果总图不需要每张都有 legend，这里设为 none
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    labs(title = paste0(id, " Score")) +
    NoAxes()
})
names(plots) <- names(mg_marker_list)

# 4. 拼图与保存

# 4.1 保存总图 (5个亚群横向排列)
# 提取一个通用的 Legend
temp_p <- FeaturePlot(tmp_obj, features = score_features[1], cols = pal) + 
  theme(legend.text = element_text(color = "transparent"))
common_legend <- cowplot::get_legend(temp_p + theme(legend.position = "right"))

combined_p <- wrap_plots(plots, nrow = 1)
final_plot <- cowplot::plot_grid(combined_p, common_legend, rel_widths = c(1, 0.05))

ggsave('results/figure/combined_MG_ModuleScore_FeaturePlot.png', 
       final_plot, height = 3, width = 16, bg = '#ffffff')

# 4.2 循环自动保存每一张单图
for (id in names(plots)) {
  # 为单图加上 Legend 后保存
  p_single <- plots[[id]] + 
    theme(legend.position = "right", 
          legend.text = element_text(size = 8))
  
  ggsave(paste0('results/figure/ModuleScore_', id, '.png'), 
         p_single, height = 5, width = 6, bg = '#ffffff')
}




# Ro/e ---------------------

tmp_obj <- qread('data/processed_data/14.mg_combined_meta.qs')

# 1. 设定：在这里指定你想要对比的两列 meta 信息
var1 <- "sub_clusters" # 纵轴：细胞亚群
var2 <- "Stage"        # 横轴：时期/分组

# 2. 计算 ROE (Relative Observed to Expected)
# 利用矩阵运算，无论 var2 有多少列（4列、10列还是2列），逻辑都通用
tab <- table(tmp_obj[[var1, drop=TRUE]],
             tmp_obj[[var2, drop=TRUE]])
row_sum <- rowSums(tab)
col_sum <- colSums(tab)
total <- sum(tab)

# 核心：计算 Expected 矩阵 (行和 * 列和 / 总数)
expected <- outer(row_sum, col_sum) / total
roe <- as.matrix(tab / expected)

# 3. 排序 (可选：确保热图显示的顺序是你想要的)
# 默认按照 factor 的顺序，也可以手动指定顺序
# rows_to_plot <- c("M0", "M1", "M2", "M3", "M4", "M5")
# cols_to_plot <- c("Child", "Adult", "Aging", "Stroke")
rows_to_plot <- rownames(roe)
cols_to_plot <- colnames(roe)

# 4. 绘图 (pheatmap 风格保持一致)
pdf(paste0("results/figure/ROE_", var1, "_", var2, ".pdf"), width = 6, height = 6)

pheatmap::pheatmap(
  roe[rows_to_plot, cols_to_plot], 
  display_numbers = TRUE,
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  color = rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(100)),
  number_format = "%.1f",
  number_color = "white",
  cellwidth = 35, 
  cellheight = 25,
  scale = 'row',
  main = paste0("ROE: ", var1, " vs ", var2),
  angle_col = "45"
)
dev.off()

