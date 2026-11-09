#工作目录------------
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine"))

# Fig.1 拟时序分析 Monocle3-------------
#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(clustree)
library(harmony)
library(SingleR)
library(dplyr)
library(tidyverse)
library(CellChat)
library(SCopeLoomR)
library(SCENIC)
library(AUCell)
library(foreach)
library(KernSmooth)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(pheatmap)
library(ggheatmap)
library(reshape2)
library(SummarizedExperiment)
library(ggpubr)
library(ggrepel)
library(SeuratWrappers)
# library(SeuratData)
library(copykat)
library(limma)
library(GSEABase)
library(GSVA)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(stringi)	
library(GOplot)
R.utils::setOption("clusterProfiler.download.method",'auto')
library(monocle3)
library(ClusterGVis)
library(openxlsx)
library(qs)


#读取数据----------------------
sco <- qread('data/processed_data/mg_cca_meta_v2.qs')

#获取表达矩阵----------------
data <- GetAssayData(sco, assay = 'RNA', layer  = 'counts')

#细胞注释信息------------
cell_metadata <- sco@meta.data

#基因名----------------
gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))

#创建一个新的 cell_data_set 对象-----------------
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

#预处理cds-------------------
cds <- preprocess_cds(cds, method = "PCA")
#降维------------------
cds <- reduce_dimension(cds, reduction_method = "UMAP",preprocess_method = 'PCA')
#画图
# p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="sub_clusters2", show_trajectory_graph = FALSE) + ggtitle('cds.umap')

#将seurat对象的UMAP导入
int.embed <- Embeddings(sco, reduction = "umap")
#排序
int.embed <- int.embed[rownames(cds@int_colData$reducedDims$UMAP),]
#导入-------------------
cds@int_colData$reducedDims$UMAP <- int.embed

#聚类分区，不同分区的细胞会进行单独的轨迹分析
cds <- cluster_cells(cds,k = 100,cluster_method = 'louvain')

#构建细胞轨迹
cds <- learn_graph(cds, learn_graph_control = list(euclidean_distance_ratio = 0.8))

#发育轨迹（拟时序）排列细胞。
cds <- order_cells(cds)

qsave(cds,'scripts/Pseudotime/data/processed_data/save_cds.qs')
cds <- qread('scripts/Pseudotime/data/processed_data/save_cds.qs')

#画拟时序umap图
p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
               label_leaves = FALSE,  label_branch_points = FALSE) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  NoAxes() +
  ggtitle("Trajectory_Pseudotime") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) # 居中标题
p
ggsave("scripts/Pseudotime/results/figure/Trajectory_Pseudotime.pdf", plot = p, width = 5, height = 4)
ggsave("scripts/Pseudotime/results/figure/Trajectory_Pseudotime.png", plot = p, width = 5, height = 4)

#寻找拟时轨迹差异基因
#根据低维嵌入和主图测试基因的差异表达
#Monocle3 引入了一种寻找此类基因的新方法，它借鉴了空间相关性分析中的一项强大技术--莫兰 I 检验。
#莫兰 I 是一种多向、多维空间自相关的测量方法。
#该统计量可以告诉您，轨迹上邻近位置的细胞对被测基因的表达水平是否相似（或不同）。
#虽然皮尔逊相关性和莫兰 I 的范围都在-1 到 1 之间，但对莫兰 I 的解释略有不同：
#+1表示邻近细胞的表达完全相似；0表示没有相关性；-1表示邻近细胞没有相关性。
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=50)
# Console中运行trace('calculateLW', edit = T, where = asNamespace("monocle3"))
# 在出现的文档中Matrix::rBind 改为rbind，或者将rbind 改为Matrix::rBind

#导出
write.csv(Track_genes, "scripts/Pseudotime/results/table/Track_genes_all.csv", row.names = F)
#保存
qsave(Track_genes,'scripts/Pseudotime/data/processed_data/Track_genes.qs')
Track_genes <- qread('scripts/Pseudotime/data/processed_data/Track_genes.qs')

#拟时基因热图
genes <- row.names(subset(Track_genes, q_value == 0 & morans_I > 0.15))

mat <- pre_pseudotime_matrix(cds_obj = cds,
                             gene_list = genes)

#保存
qsave(mat,'scripts/Pseudotime/data/processed_data/mat.qs')
mat <- qread('scripts/Pseudotime/data/processed_data/mat.qs')

head(mat[1:5,1:5])


# # -----------
cds <- qread('scripts/Pseudotime/data/processed_data/save_cds.qs')
mat <- qread('scripts/Pseudotime/data/processed_data/mat.qs')
# 
# #选取合适的clusters
# 
# library(ClusterGVis)
# 
# #kmeans
# ck <- clusterData(exp = mat,
#                   cluster.method = "kmeans",
#                   cluster.num = 5)
# 
# qsave(ck,'scripts/Pseudotime/data/processed_data/ck.qs')
ck <- qread('scripts/Pseudotime/data/processed_data/ck.qs')



my_base_colors <- rev(c("#810000", "#CE1212", "#F05454", "#ffd06f",
                        "#ffe6b7", "#aadce0", "#72bcd5", "#528fad",
                        "#376795", "#1e466e"))


heatmap_cells <- colnames(mat)
my_cell_types <- cds@colData[heatmap_cells, "CellType1"]
my_cell_types <- as.character(my_cell_types)
print(length(my_cell_types) == ncol(mat))
unique_types <- unique(my_cell_types)

my_anno_colors <- c(
  "MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A",
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628"
)

names(my_anno_colors) <- c("MG2", "MG5", "MG0", "MG3", "MG4", "MG1")
my_expanded_cell_colors <- my_anno_colors[my_cell_types]
print(length(my_expanded_cell_colors) == 21244)

# add line annotation
pdf('scripts/Pseudotime/results/figure/monocle3_2.pdf',height = 10,width = 8,onefile = F)
visCluster(object = ck,
           add.sampleanno = F,
           plot.type = "both",
           cluster_rows = FALSE,
           cluster.order = c(1, 2, 5, 4, 3),
           markGenes = sample(rownames(mat),30,replace = F),
           ht.col.list = list(
             col_range = seq(-2, 2, length.out = 10), 
             col_color = my_base_colors
           ),
           line.side = "left",
           )
dev.off()

# 带GO的热图--------------
enrich <- enrichCluster(object = ck,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 100,
                        seed = 24)

# 获取所有不重复的 cluster (group) 名称
groups <- unique(enrich$group)

# 循环遍历每一个 cluster
for (g in groups) {
  # 提取当前 cluster 的数据
  tmp_data <- enrich |> dplyr::filter(group == g)
  
  # 设置要保存的文件名，例如 "GO_Result_C1.csv"
  file_name <- paste0(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Pseudotime/results/table/GO_Result_"), g, ".csv")
  
  # 导出为 CSV 文件
  write.csv(tmp_data, file = file_name, row.names = FALSE)
  
  cat("已保存:", file_name, "\n")
}


# check
head(enrich,3)

palette = c("Grays", "Light Grays", "Blues 2", "Blues 3", "Purples 2")


my_base_colors <- rev(c("#810000", "#CE1212", "#F05454", "#ffd06f",
                        "#ffe6b7", "#aadce0", "#72bcd5", "#528fad",
                        "#376795", "#1e466e"))

# loop
lapply(seq_along(unique(enrich$group)), function(x){
  # go plot
  tmp <- enrich |> dplyr::filter(group == unique(enrich$group)[x]) |>
    dplyr::arrange(desc(pvalue))
  
  tmp$Description <- factor(tmp$Description, levels = tmp$Description)
  
  # plot
  p <-
    ggplot(tmp) +
    geom_col(aes(x = -log10(pvalue), y = Description, fill = -log10(pvalue)),
             width = 0.75) +
    geom_line(aes(x = log10(ratio), y = as.numeric(Description)), color = "grey50") +
    geom_point(aes(x = log10(ratio), y = Description), size = 3, color = "orange") +
    theme_bw() +
    scale_y_discrete(position = "right",
                     labels = function(x) stringr::str_wrap(x, width = 40)) +
    scale_x_continuous(sec.axis = sec_axis(~., name = "log10(ratio)")) +
    colorspace::scale_fill_binned_sequential(palette = palette[x]) +
    ylab("")
  
  return(p)
}) -> gglist

# assign names
names(gglist) <- paste("C", 1:5, sep = "")

pdf('scripts/Pseudotime/results/figure/Trajectory_Heatmap_with_GO.pdf', height = 10, width = 14, onefile = F)

visCluster(object = ck,
           add.sampleanno = F,
           plot.type = "both",
           cluster_rows = FALSE,
           cluster.order = c(1, 2, 5, 4, 3),
           markGenes = sample(rownames(mat), 30, replace = F),
           ht.col.list = list(
             col_range = seq(-2, 2, length.out = 10), 
             col_color = my_base_colors
           ),
           line.side = "left",
           ggplot.panel.arg = c(4.5,0.5,16,"grey90",NA),
           gglist = gglist
)
dev.off()


# 手动画注释的bar----------------------
library(ggplot2)

# 1. 提取所有细胞的拟时间值
pt_values <- pseudotime(cds)

# 2. 过滤掉没有算在轨迹里的细胞（拟时间为 Infinite 的细胞）
# 通常 pre_pseudotime_matrix 内部也会自动丢弃这些细胞
pt_values <- pt_values[is.finite(pt_values)]

# 3. 【核心步骤】按照拟时间从小到大，对细胞条形码（Barcode）进行排序
ordered_barcodes <- names(sort(pt_values))

# 建议运行这行检查一下，应该是 TRUE
print(paste("排序后的细胞数与 mat 列数是否一致:", length(ordered_barcodes) == ncol(mat)))

# 4. 按照这个排序好的细胞顺序，去提取对应的细胞类型
# 这样提取出来的细胞类型，就和 mat 里的列 100% 一一对应了！
sorted_cell_types <- cds@colData[ordered_barcodes, "CellType1"]

# 5. 构建绘图数据框
df_bar <- data.frame(
  pos = 1:length(sorted_cell_types), # x 轴 1 到两万多
  y = 1,
  # 加上 "MG" 前缀，并转化为 factor，匹配你的颜色字典
  CellType =sorted_cell_types) 

# 6. 画图
p_bar <- ggplot(df_bar, aes(x = pos, y = y, fill = CellType)) +
  geom_raster() + # geom_raster 处理两万个点极快
  scale_fill_manual(values = my_anno_colors) + 
  theme_void() + 
  theme(
    legend.position = "none", 
    plot.margin = margin(0, 0, 0, 0) 
  )

# 7. 保存为又长又扁的 PDF
ggsave("scripts/Pseudotime/results/figure/cell_type_bar_only.pdf", 
       plot = p_bar, 
       width = 8,     # 宽度建议和你的热图保持一致
       height = 0.4, 
       limitsize = FALSE)

print("按照拟时间严格排序的颜色条已生成！")

# 指定基因画jitter plot-----------
# genes <- c('TRIM2','RNF220','ST18','MBP','BCAS1','MOBP','PLP1')
genes <- c('ACSL1','CTSL','MOG','CTNNA2','FGF14')

my_anno_colors <- c(
  "MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A",
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628"
)

p <- monocle3::plot_genes_in_pseudotime(
  cds[genes,], 
  color_cells_by = "CellType1", 
  min_expr = 0.5, 
  ncol = 7,
  panel_order = genes  # <--- 核心参数：直接指定分面顺序
) + 
  scale_colour_manual(values = my_anno_colors)

ggsave("scripts/Pseudotime/results/figure/Genes_Jitterplot_2.pdf", plot = p, width = 16, height = 2)


# 提取每一个cluster中的基因--------------------
# =========================================================================
# 1. 加载必要的 R 包
# =========================================================================
library(qs)
library(ClusterGVis)
# library(org.Hs.eg.db) # 不做GO分析，可以不再加载此包

# 读取数据
ck <- qread('scripts/Pseudotime/data/processed_data/ck.qs')

# =========================================================================
# 2. 颜色与注释设置
# =========================================================================
my_base_colors <- rev(c("#810000", "#CE1212", "#F05454", "#ffd06f",
                        "#ffe6b7", "#aadce0", "#72bcd5", "#528fad",
                        "#376795", "#1e466e"))

heatmap_cells <- colnames(mat)
my_cell_types <- cds@colData[heatmap_cells, "CellType1"]
my_cell_types <- as.character(my_cell_types)
print(length(my_cell_types) == ncol(mat))
unique_types <- unique(my_cell_types)

my_anno_colors <- c(
  "MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A",
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628"
)

names(my_anno_colors) <- c("MG2", "MG5", "MG0", "MG3", "MG4", "MG1")
my_expanded_cell_colors <- my_anno_colors[my_cell_types]
print(length(my_expanded_cell_colors) == 21244)

# =========================================================================
# 3. 绘制带有线图聚类的热图并保存
# =========================================================================
pdf('scripts/Pseudotime/results/figure/monocle3_2.pdf', height = 10, width = 8, onefile = F)
visCluster(object = ck,
           add.sampleanno = F,
           plot.type = "both",
           cluster_rows = FALSE,
           cluster.order = c(1, 2, 5, 4, 3),
           markGenes = sample(rownames(mat), 30, replace = F),
           ht.col.list = list(
             col_range = seq(-2, 2, length.out = 10), 
             col_color = my_base_colors
           ),
           line.side = "left")
dev.off()

# =========================================================================
# 4. 提取每个 cluster 的基因，并保存为 CSV 文件
# =========================================================================
# 提取 cluster 和 gene 这两列，并去除重复项
unique_gene_cluster <- unique(ck$long.res[, c("cluster", "gene")])

# 按照 cluster 分组，把基因拆分成 List (列表)
cluster_genes_list <- split(unique_gene_cluster$gene, unique_gene_cluster$cluster)

# 在控制台打印查看每个 cluster 的基因数量
print("各个 Cluster 包含的基因数量：")
print(lapply(cluster_genes_list, length))

# 为了导出为整齐的 CSV 表格，找出基因最多的 cluster 的长度
max_len <- max(sapply(cluster_genes_list, length))

# 用 NA 补齐较短的 cluster 列，将其转换为 DataFrame
genes_df <- data.frame(lapply(cluster_genes_list, function(x) {
  length(x) <- max_len  # 长度不足的会自动填补 NA
  return(x)
}))

# 规范化列名 (例如：Cluster_1, Cluster_2, ...)
colnames(genes_df) <- paste0("Cluster_", names(cluster_genes_list))

# 导出基因列表为 CSV 文件
write.csv(genes_df, "scripts/Pseudotime/results/table/monocle3_cluster_genes.csv", row.names = FALSE)

print("基因提取完成！文件已保存至：scripts/Pseudotime/results/table/monocle3_cluster_genes.csv")
