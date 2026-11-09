#!/usr/bin/env Rscript
# Figure-ordered code from the reviewed manuscript analysis package.
# Each panel contains its original quoted R code.
# Source this chapter, then use mg_list_panels() or mg_show_panel("PANEL_ID").
# Source paths in panel records identify the original analysis scripts.
.mg_panels <- new.env(parent=emptyenv())
.mg_order <- character()
mg_register_panel <- function(id, code) {
  if (!id %in% .mg_order) .mg_order <<- c(.mg_order,id)
  assign(id, code, envir=.mg_panels)
}
mg_list_panels <- function() .mg_order
mg_panel_code <- function(id) {
  if (!exists(id,envir=.mg_panels,inherits=FALSE)) stop("Unknown panel: ",id)
  get(id,envir=.mg_panels,inherits=FALSE)
}
mg_show_panel <- function(id) {
  blocks <- mg_panel_code(id)
  if (is.null(blocks)) {cat(id, ": No R block is included for this panel.\n");return(invisible(NULL))}
  for (block in blocks) {
    cat("# Source: ",block$source," [actual ",block$ranges,"; requested ",block$requested_ranges,"]\n# Role: ",block$role,"\n",sep="")
    cat(paste(deparse(block$code,width.cutoff=120L),collapse="\n"),"\n\n")
  }
  invisible(blocks)
}


# Figure 5

# F5A  Cell-cycle phase on the human microglia UMAP
# Kind: computational
mg_register_panel("F5A", list(
  list(source="source/project/code_upload/Fig.01.r", ranges="930-947", requested_ranges="930-947", role="prepared-upload implementation", code=quote({
# Fig.14: Cell Cycle -----------------------------------------------------------
cell_cycle_path <- file.path(figure_dir, "cellcycle_dimplot.pdf")
phase_colors <- c("G1" = "#1F77B4", "S" = "#FF7F0E", "G2M" = "#2CA02C")

tmp_obj <- qread(input_obj_v2)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

DefaultAssay(tmp_obj) <- "RNA"
tmp_obj <- CellCycleScoring(
  tmp_obj,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)

p_cycle <- DimPlot(tmp_obj, group.by = "Phase", cols = phase_colors)
ggsave(cell_cycle_path, p_cycle, width = 6, height = 5)
  })),
  list(source="source/project/Reference/combine/scripts/04.module.r", ranges="883-895", requested_ranges="883-895", role="original implementation", code=quote({
# MODULE：cell cycle------------------
tmp_obj <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

DefaultAssay(tmp_obj) <- "RNA"
tmp_obj <- CellCycleScoring(tmp_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# 活性递进配色
phase_colors <- c("G1" = "#1F77B4", "S" = "#FF7F0E", "G2M" = "#2CA02C")
p_cycle <- DimPlot(tmp_obj, group.by = 'Phase', cols = phase_colors)
ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/results/figure/cellcycle_dimplot.pdf"),p_cycle,width = 6,height = 5)
  }))
))

# F5B  Monocle3 pseudotime trajectory
# Kind: computational
mg_register_panel("F5B", list(
  list(source="source/project/code_upload/Pseudotime.r", ranges="14-42", requested_ranges="14-42", role="counts to CDS, existing UMAP transfer, trajectory and plotting", code=quote({
sco <- qread('data/processed_data/mg_cca_meta_v2.qs')
data <- GetAssayData(sco, assay = 'RNA', layer  = 'counts')
cell_metadata <- sco@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))

cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, method = "PCA")
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = 'PCA')

int.embed <- Embeddings(sco, reduction = "umap")
int.embed <- int.embed[rownames(cds@int_colData$reducedDims$UMAP), ]
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds, k = 100, cluster_method = 'louvain')
cds <- learn_graph(cds, learn_graph_control = list(euclidean_distance_ratio = 0.8))
cds <- order_cells(cds) 

qsave(cds, 'scripts/Pseudotime/data/processed_data/save_cds.qs')

# Fig.01 ---------------
p_traj <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                     label_leaves = FALSE, label_branch_points = FALSE) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  NoAxes() +
  ggtitle("Trajectory Pseudotime") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("scripts/Pseudotime/results/figure/Trajectory_Pseudotime.pdf", plot = p_traj, width = 5, height = 4)
ggsave("scripts/Pseudotime/results/figure/Trajectory_Pseudotime.png", plot = p_traj, width = 5, height = 4)
  })),
  list(source="source/project/Reference/combine/scripts/Pseudotime/scripts/module_pseudotime.r", ranges="53-105", requested_ranges="53-105", role="original implementation and saved-CDS reload", code=quote({
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
  }))
))

# F5C  Pseudotime heatmap, cell-type bar, and selected GO labels
# Kind: computational_with_manual_layout
mg_register_panel("F5C", list(
  list(source="source/project/code_upload/Pseudotime.r", ranges="44-155", requested_ranges="44-155", role="graph_test, q=0/Moran I >0.15, five clusters, GO, heatmap and ordered-cell bar", code=quote({
Track_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 50)
write.csv(Track_genes, "scripts/Pseudotime/results/table/Track_genes_all.csv", row.names = FALSE)
qsave(Track_genes, 'scripts/Pseudotime/data/processed_data/Track_genes.qs')

genes <- row.names(subset(Track_genes, q_value == 0 & morans_I > 0.15))
mat <- pre_pseudotime_matrix(cds_obj = cds, gene_list = genes)
qsave(mat, 'scripts/Pseudotime/data/processed_data/mat.qs')

ck <- clusterData(exp = mat, cluster.method = "kmeans", cluster.num = 5)
qsave(ck, 'scripts/Pseudotime/data/processed_data/ck.qs')

my_base_colors <- rev(c("#810000", "#CE1212", "#F05454", "#ffd06f",
                        "#ffe6b7", "#aadce0", "#72bcd5", "#528fad",
                        "#376795", "#1e466e"))

my_anno_colors <- c(
  "MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A",
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628"
)

# Fig.02 ---------------
pdf('scripts/Pseudotime/results/figure/monocle3_2.pdf', height = 10, width = 8, onefile = FALSE)
visCluster(object = ck,
           add.sampleanno = FALSE,
           plot.type = "both",
           cluster_rows = FALSE,
           cluster.order = c(1, 2, 5, 4, 3),
           markGenes = sample(rownames(mat), 30, replace = FALSE),
           ht.col.list = list(
             col_range = seq(-2, 2, length.out = 10), 
             col_color = my_base_colors
           ),
           line.side = "left")
dev.off()

enrich <- enrichCluster(object = ck,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 100,
                        seed = 24)

groups <- unique(enrich$group)
for (g in groups) {
  tmp_data <- enrich |> dplyr::filter(group == g)
  file_name <- paste0("scripts/Pseudotime/results/table/GO_Result_", g, ".csv")
  write.csv(tmp_data, file = file_name, row.names = FALSE)
}

palette <- c("Grays", "Light Grays", "Blues 2", "Blues 3", "Purples 2")

gglist <- lapply(seq_along(unique(enrich$group)), function(x){
  tmp <- enrich |> dplyr::filter(group == unique(enrich$group)[x]) |>
    dplyr::arrange(desc(pvalue))
  
  tmp$Description <- factor(tmp$Description, levels = tmp$Description)
  
  p <- ggplot(tmp) +
    geom_col(aes(x = -log10(pvalue), y = Description, fill = -log10(pvalue)), width = 0.75) +
    geom_line(aes(x = log10(ratio), y = as.numeric(Description)), color = "grey50") +
    geom_point(aes(x = log10(ratio), y = Description), size = 3, color = "orange") +
    theme_bw() +
    scale_y_discrete(position = "right", labels = function(x) stringr::str_wrap(x, width = 40)) +
    scale_x_continuous(sec.axis = sec_axis(~., name = "log10(ratio)")) +
    colorspace::scale_fill_binned_sequential(palette = palette[x]) +
    ylab("")
  
  return(p)
})

names(gglist) <- paste("C", 1:5, sep = "")

# Fig.03 ---------------
pdf('scripts/Pseudotime/results/figure/Trajectory_Heatmap_with_GO.pdf', height = 10, width = 14, onefile = FALSE)
visCluster(object = ck,
           add.sampleanno = FALSE,
           plot.type = "both",
           cluster_rows = FALSE,
           cluster.order = c(1, 2, 5, 4, 3),
           markGenes = sample(rownames(mat), 30, replace = FALSE),
           ht.col.list = list(
             col_range = seq(-2, 2, length.out = 10), 
             col_color = my_base_colors
           ),
           line.side = "left",
           ggplot.panel.arg = c(4.5, 0.5, 16, "grey90", NA),
           gglist = gglist)
dev.off()

# Fig.04 ---------------
pt_values <- pseudotime(cds)
pt_values <- pt_values[is.finite(pt_values)]
ordered_barcodes <- names(sort(pt_values))
sorted_cell_types <- cds@colData[ordered_barcodes, "CellType1"]

df_bar <- data.frame(
  pos = 1:length(sorted_cell_types), 
  y = 1,
  CellType = sorted_cell_types
) 

p_bar <- ggplot(df_bar, aes(x = pos, y = y, fill = CellType)) +
  geom_raster() + 
  scale_fill_manual(values = my_anno_colors) + 
  theme_void() + 
  theme(
    legend.position = "none", 
    plot.margin = margin(0, 0, 0, 0) 
  )

ggsave("scripts/Pseudotime/results/figure/cell_type_bar_only.pdf", plot = p_bar, width = 8, height = 0.4, limitsize = FALSE)
  })),
  list(source="source/project/Reference/combine/scripts/Pseudotime/scripts/module_pseudotime.r", ranges="114-215;256-317", requested_ranges="114-212;256-317", role="saved graph_test/matrix/ck and GO tables; heatmap cluster order and bar", code=quote({
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
  }))
))

# F5D  TRIM2, RNF220, MBP, MOBP, ST18 and PLP1 along pseudotime
# Kind: computational_with_manual_layout
mg_register_panel("F5D", list(
  list(source="source/project/code_upload/Pseudotime.r", ranges="157-169", requested_ranges="157-169", role="original gene-pseudotime plot", code=quote({
# Fig.05 ---------------
genes_to_plot <- c('TRIM2','RNF220','ST18','MBP','BCAS1','MOBP','PLP1') 

p_jitter <- monocle3::plot_genes_in_pseudotime(
  cds[genes_to_plot, ], 
  color_cells_by = "CellType1", 
  min_expr = 0.5, 
  ncol = 7,
  panel_order = genes_to_plot 
) + 
  scale_colour_manual(values = my_anno_colors)

ggsave("scripts/Pseudotime/results/figure/Genes_Jitterplot.pdf", plot = p_jitter, width = 16, height = 2)
  })),
  list(source="source/project/Reference/combine/scripts/Pseudotime/scripts/module_pseudotime.r", ranges="321-339", requested_ranges="321-339", role="same plotting function but active gene list changed", code=quote({
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
  }))
))

# F5E  Mouse-myelin debris stimulation of HMC3 experimental schematic
# Kind: experimental_schematic
mg_register_panel("F5E", NULL)

# F5F  Iba1/Trim2/MOG immunofluorescence in control versus myelin-exposed HMC3
# Kind: experimental_micrographs
mg_register_panel("F5F", NULL)

# F5G  TRIM2/MBP/MOG qPCR treatment-versus-control barplots
# Kind: experimental_figure
mg_register_panel("F5G", NULL)

# F5H  OLM gene heatmap in snRNA-seq and six-sample bulk RNA-seq
# Kind: computational
mg_register_panel("F5H", list(
  list(source="source/project/Reference/combine/scripts/Heatmap.r", ranges="1-55;61-108", requested_ranges="1-55;61-108", role="left: AverageExpression RNA by CellType1, gene list and row scaling; right: vsd.rds, same 26-gene list, row scaling", code=quote({
# MODULE 9: Microglia Other Function Gene Heatmap ------------------------------
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/"))
rm(list=ls())

## 9.1 Parameters --------------------------------------------------------------
input_file      <- 'data/processed_data/mg_cca_meta_v2.qs'
output_fig_dir  <- 'results/figure/'
output_filename <- paste0(output_fig_dir, 'bulk_gene_heatmap.pdf')

dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

target_genes <- c("PLLP", "FAM13C", "ITGA2", "NCKAP5", "TTLL11", 
                  "ZNF608", "AK5", "CORO2B", "TRIM9", "DIXDC1", 
                  "GAB1", "ARNT2", "KIAA1958", "TNRC6C", "RNF144A", 
                  "SH3TC2", "TRPC1", "PLEKHB1", "TTLL7", "TRIM2", 
                  "ENOSF1", "PRKD1", "LRRC1", "SEMA4D", "FAXDC2", 
                  "SLC13A3")


## 9.2 Process Data ------------------------------------------------------------
library(Seurat)
library(pheatmap)
library(qs)

if(!file.exists(input_file)) stop("Input file not found.")
tmp_obj <- qread(input_file)

# Calculate Average Expression
avg_exp_list <- AverageExpression(
  tmp_obj,
  group.by = "CellType1",
  assays = "RNA"
)
avg_mat <- as.data.frame(avg_exp_list$RNA)

# Filter and Sort Genes
# Ensure we only define genes that exist in the matrix and maintain input order
valid_genes <- intersect(target_genes, rownames(avg_mat))
plot_data   <- avg_mat[valid_genes, , drop = FALSE]
plot_data   <- plot_data[match(valid_genes, rownames(plot_data)), ]


## 9.3 Generate & Save Plot ----------------------------------------------------
pdf(output_filename, width = 3, height = 3.5)

pheatmap(
  plot_data,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  show_rownames = TRUE,
  legend = TRUE
)

dev.off()

# MODULE 9: Microglia Other Function Gene Heatmap ------------------------------
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/"))
rm(list=ls())

## 9.1 Parameters --------------------------------------------------------------
output_fig_dir  <- 'results/figure/'
output_filename <- paste0(output_fig_dir, 'bulk_heatmap.pdf')

dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

target_genes <- c("PLLP", "FAM13C", "ITGA2", "NCKAP5", "TTLL11", 
                  "ZNF608", "AK5", "CORO2B", "TRIM9", "DIXDC1", 
                  "GAB1", "ARNT2", "KIAA1958", "TNRC6C", "RNF144A", 
                  "SH3TC2", "TRPC1", "PLEKHB1", "TTLL7", "TRIM2", 
                  "ENOSF1", "PRKD1", "LRRC1", "SEMA4D", "FAXDC2", 
                  "SLC13A3")


## 9.2 Process Data ------------------------------------------------------------
library(Seurat)
library(pheatmap)
library(qs)

vsd <- readRDS(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/processed_data/vsd.rds"))
mat <- assay(vsd)[target_genes, ]

# Filter and Sort Genes
# Ensure we only define genes that exist in the matrix and maintain input order
valid_genes <- intersect(target_genes, rownames(mat))
plot_data   <- mat[valid_genes, , drop = FALSE]
plot_data   <- plot_data[match(valid_genes, rownames(plot_data)), ]


## 9.3 Generate & Save Plot ----------------------------------------------------
pdf(output_filename, width = 3, height = 3.5)

pheatmap(
  plot_data,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  show_rownames = TRUE,
  legend = TRUE
)

dev.off()

message("Module 9 Done.")
  })),
  list(source="source/project/Reference/combine/scripts/BULK/scripts/DESeq2.r", ranges="12-30", requested_ranges="12-30", role="upstream bulk count input, sample selection, filtering and VST", code=quote({
# 1. 数据加载---------------------
raw_counts <- read.csv(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/processed_data/All_gene_counts.list"),sep = '\t',header=TRUE,row.names=1)
raw_counts <- raw_counts[,c(-4,-8)]

# 4. 数据预处理与质控----------------
col_data <- read.csv(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/saved_data/sample_info.txt"),sep = '\t',header=TRUE,row.names = 1)
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = col_data,
  design = ~ group)

# 过滤低表达基因
keep <- rowSums(counts(dds) >= 10) >= 4 # 我们的过滤标准（至少4个样本中计数≥10）在统计效能和基因保留之间取得了良好平衡
dds <- dds[keep,]
cat("过滤后保留", nrow(dds), "个基因\n")

# 方差稳定变换
vsd <- vst(dds, blind = FALSE) # vst（方差稳定变换）：计算速度更快，特别适用于大样本量。rlog（正则化对数变换）：计算速度较慢，但对于小样本量可能表现稍好。
saveRDS(vsd,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/processed_data/vsd.rds"))
  }))
))

# F5I  GO barplot of five terms from preselected workbook
# Kind: computational_precomputed_table
mg_register_panel("F5I", list(
  list(source="source/project/Reference/combine/scripts/20260519_go_barplot.r", ranges="1-15", requested_ranges="1-15", role="exact same five terms, Count and p.adjust color as final panel", code=quote({
library(readxl)
library(ggplot2)
library(dplyr)

go_data <- read_excel(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/saved_data/20260519_go.xlsx"))

p <- ggplot(go_data, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_col(width = 0.65) +
  scale_fill_gradient(low = "#2E86AB", high = "#A2BCE0") +
  labs(title = "GO Enrichment Bar Plot", x = "Gene Count", y = "GO Term", fill = "Adjusted p-value") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10), panel.grid = element_blank())

ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/results/figure/20260519_go_barplot.pdf"), plot = p, width = 8, height = 4, dpi = 300)
ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/results/figure/20260519_go_barplot.png"), plot = p, width = 8, height = 4, dpi = 300)
  }))
))
