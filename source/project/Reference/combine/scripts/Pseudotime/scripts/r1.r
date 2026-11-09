# 工作目录
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine"))

# 仅加载该流程实际需要的包
library(Seurat)
library(monocle3)
library(qs)
library(ggplot2)
library(ClusterGVis)

# 1. 读取数据并构建 Monocle3 对象
sco <- qread('data/processed_data/mg_cca_meta_v2.qs')
data <- GetAssayData(sco, assay = 'RNA', layer = 'counts')
cell_metadata <- sco@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

# 2. 降维与 Seurat UMAP 坐标导入
cds <- preprocess_cds(cds, method = "PCA")
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = 'PCA')
int.embed <- Embeddings(sco, reduction = "umap")
int.embed <- int.embed[rownames(cds@int_colData$reducedDims$UMAP), ]
cds@int_colData$reducedDims$UMAP <- int.embed

# 3. 聚类、构建并学习轨迹
cds <- cluster_cells(cds, k = 100, cluster_method = 'louvain')
cds <- learn_graph(cds, learn_graph_control = list(euclidean_distance_ratio = 0.8),minimal_branch_len = 15, ncenter = 250 )
cds <- order_cells(cds) # 此处会弹出交互界面让你选择起点
qsave(cds, 'scripts/Pseudotime/data/processed_data/cds.qs')

# 4. 绘制拟时序轨迹图
p <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                label_leaves = FALSE, label_branch_points = FALSE) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  NoAxes() +
  ggtitle("Trajectory_Pseudotime") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("scripts/Pseudotime/results/figure/Trajectory_Pseudotime.pdf", plot = p, width = 5, height = 4)

# 5. 寻找拟时轨迹差异基因
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=50)
write.csv(Track_genes, "scripts/Pseudotime/results/table/Track_genes_all.csv", row.names = F)
qsave(Track_genes, 'scripts/Pseudotime/data/processed_data/Track_genes.qs')

# 6. 提取显著差异基因并生成热图矩阵
genes <- row.names(subset(Track_genes, q_value == 0 & morans_I > 0.15))
mat <- pre_pseudotime_matrix(cds_obj = cds, gene_list = genes)
qsave(mat, 'scripts/Pseudotime/data/processed_data/mat.qs')

# 7. 基因聚类与热图可视化
ck <- qread('scripts/Pseudotime/data/processed_data/ck.qs') # 读取已跑好的 kmeans 聚类结果
my_base_colors <- rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", "#ffe6b7", 
                        "#aadce0", "#72bcd5", "#528fad", "#376795", "#1e466e"))

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
           )
)
dev.off() # 关闭图形设备，保存 PDF