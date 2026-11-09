# CellChat流程

# 全局设置----------------
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/CellChat"))
rm(list = ls())

# library------------------
library(CellChat)
library(Seurat)
library(qs)
library(ggplot2)
library(patchwork)
library(igraph)
library(ComplexHeatmap)

# 加载数据---------------------
seu <- qread('data/processed_data/stroke_cellchat.qs')
seu.data.input <- GetAssayData(seu, assay = "RNA", slot = "data")
seu.identity <- subset(seu@meta.data, select = "cellchat_type")
cellchat <- createCellChat(object = seu.data.input, meta = seu.identity, group.by = "cellchat_type")

cellchat <- addMeta(cellchat, meta = seu.identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") 
groupSize <- as.numeric(table(cellchat@idents))

# 设置配体-受体相互作用数据库
CellChatDB <- CellChatDB.human  
CellChatDB.use <- subsetDB(CellChatDB, search = c("Cell-Cell Contact"))
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = FALSE,population.size = T)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# weight<-as.data.frame(cellchat@net$weight)
# write.table(weight,file="results/table/weight_con_DLB.xls",sep="\t")
# count<-as.data.frame(cellchat@net$count)
# write.table(count,file="results/table/count_con_DLB.xls",sep="\t")
# netP<-as.data.frame(cellchat@netP[["prob"]])
# write.table(netP,file="results/table/netP_con_DLB.xls",sep="\t")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# 保存-----------
qsave(cellchat,'data/processed_data/con_stroke_cellchat.qs')

# cellchat <- qread('data/processed_data/con_mouse_cellchat.qs')
# netVisual_circle(cellchat@net$weight, vertex.weight = 13, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# netVisual_circle(cellchat@net$weight, vertex.weight = 13, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",sources.use = c("OLIG"),targets.use = c("cluster0","cluster1","cluster2","cluster3"))
# netVisual_circle(cellchat@net$weight, vertex.weight = 13, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",sources.use = c("cluster0","cluster1","cluster2","cluster3"))
# 
# netVisual_bubble(cellchat, sources.use = c("cluster2"), remove.isolate = TRUE,  angle.x = 45)
# netVisual_bubble(cellchat, sources.use = c("OLIG"), targets.use = c("cluster0","cluster1","cluster2","cluster3"),remove.isolate = TRUE,  angle.x = 45)
# netVisual_bubble(cellchat, sources.use = c("cluster0","cluster1","cluster2","cluster3"), targets.use = c("OLIG") ,remove.isolate = TRUE,  angle.x = 45)
# 
# 
# 
# all_cell_types <- levels(cellchat@idents)
# targets_to_use <- setdiff(all_cell_types, c("cluster0", "cluster1", "cluster2", "cluster3"))
# netVisual_bubble(
#   cellchat, 
#   sources.use = "cluster2", 
#   targets.use = targets_to_use, 
#   remove.isolate = TRUE,  
#   angle.x = 45
# )
# 
# 
# MODULE : CellChat Interaction Visualization -------------------

## 1 Parameters --------------------------------------------------------------
# 加载所需包 (请确保已安装)
library(qs)
library(CellChat)
library(ggplot2)
library(patchwork)

# 文件路径设置
input_file <- 'data/processed_data/con_stroke_cellchat.qs'
out_dir <- 'results/figure/'
out_file_circle_1  <- paste0(out_dir, "stroke_cellchat_circle_multi.pdf")
out_file_circle_2  <- paste0(out_dir, "stroke_cellchat_circle_olig.pdf")
out_file_circle_3  <- paste0(out_dir, "stroke_cellchat_circle_subset.pdf")
out_file_bubble    <- paste0(out_dir, "stroke_bubble_cellchat.pdf")

# 定义细胞群组 (方便后续随意更改)
clusters_source_multi <- c("MG0", "MG1", "MG2", "MG3", "MG4", "MG5")
cluster_olig          <- c("OLIG")
cluster_source_single <- c("MG0")
clusters_target_main  <- c("Neuron", "Astrocyte", "OLIG", "OPC",  "Ependymal", "Vascular", "T cell")

# 子集分析群组 (将上述群组结合)
subset_celltypes <- c(cluster_olig, clusters_source_multi)

# 气泡图关注的受体-配体对
lr_pairs_list <- c(
  "NCAM1_NCAM2", "NCAM1_NCAM1", "JAM3_JAM3", "CD6_ALCAM", "ALCAM_CD6",
  "CNTN2_CNTN2_CNTNAP2", "CNTN2_CNTN2", "MAG_MAG",
  "SEMA4D_PLXNB1", "CDH1_ITGA1_ITGB1", "CD86_CD28",
  "CD80_CD28", "HLA-B_CD8A","HLA-A_CD8A", "ICAM1_ITGAX_ITGB2"
)


# 统一的颜色配置 (通过名称匹配，支持自动提取子集颜色)
cell_colors <- c(
  "MG0" = "#4DAF4A", 
  "MG1" = "#A65628",
  "MG2" = "#377EB8", 
  "MG3" = "#984EA3", 
  "MG4" = "#E6AB02", 
  "MG5" = "#8B0000", 
  "Neuron"         = "#E31A1C",
  "Astrocyte"      = "#FF7F00",
  "OLIG"           = "#1F78B4",
  "OPC"            = "#A6CEE3",
  "Ependymal"      = "#6A3D9A",
  "Vascular"       = "#FB9A99",
  "T cell"         = "#FDBF6F"
)


## 2 Process Data ------------------------------------------------------------
# 检查并创建输出目录
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 加载数据
cellchat <- qread(input_file)

# 提取并处理权重矩阵 (用于子集圆圈图)
cell_weight <- cellchat@net$weight

net_weight_no_self <- cellchat@net$weight
diag(net_weight_no_self) <- 0

mat_subset_rows <- cell_weight[rownames(cell_weight) %in% subset_celltypes, ]
mat_subset <- mat_subset_rows[, colnames(mat_subset_rows) %in% subset_celltypes]
diag(mat_subset) <- 0

# 处理气泡图所需的 L-R pairs 数据框
pairLR.use_df <- data.frame(interaction_name = lr_pairs_list, stringsAsFactors = FALSE)

# 提取子集作图专属颜色字典
new_col_subset <- cell_colors[names(cell_colors) %in% subset_celltypes]


## 3 Generate Plots ----------------------------------------------------------
# --- Base R Plots (由于Base R直接写入设备，这里用函数封装作图逻辑) ---
plot_circle_1 <- function() {
  netVisual_circle(net_weight_no_self,
                   arrow.size = 0.5, arrow.width = 0.8, vertex.weight = 13,
                   weight.scale = TRUE, label.edge = FALSE,
                   color.use = cell_colors, title.name = "Interaction weights/strength",
                   sources.use = clusters_source_multi)
}

plot_circle_2 <- function() {
  netVisual_circle(net_weight_no_self,
                   arrow.size = 0.5, arrow.width = 0.8, vertex.weight = 13,
                   weight.scale = TRUE, label.edge = FALSE,
                   color.use = cell_colors, title.name = "Interaction weights/strength",
                   sources.use = cluster_olig)
}

plot_circle_3 <- function() {
  netVisual_circle(mat_subset,
                   arrow.size = 1, arrow.width = 0.8, vertex.weight = 13,
                   weight.scale = TRUE, label.edge = FALSE,
                   color.use = new_col_subset, title.name = "Interaction weights/strength",
                   edge.width.max = 15)
}

# --- ggplot2 Plots (气泡图) ---
p1 <- netVisual_bubble(cellchat,
                       sources.use = cluster_source_single,
                       targets.use = clusters_target_main,
                       pairLR.use = pairLR.use_df,
                       remove.isolate = FALSE, show.legend = FALSE) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

p2 <- netVisual_bubble(cellchat,
                       sources.use = clusters_source_multi,
                       targets.use = cluster_olig,
                       pairLR.use = pairLR.use_df,
                       remove.isolate = FALSE, show.legend = FALSE) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p3 <- netVisual_bubble(cellchat,
                       sources.use = cluster_olig,
                       targets.use = clusters_source_multi,
                       pairLR.use = pairLR.use_df,
                       remove.isolate = FALSE) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

# 拼接气泡图
combined_bubble <- p1 | p2 | p3


## 4 Save Results ------------------------------------------------------------
# 保存圆圈图 1
pdf(out_file_circle_1, width = 5, height = 5)
plot_circle_1()
invisible(dev.off())

# 保存圆圈图 2 (修复了原代码中覆盖同名文件的问题)
pdf(out_file_circle_2, width = 5, height = 5)
plot_circle_2()
invisible(dev.off())

# 保存圆圈图 3 (子集)
pdf(out_file_circle_3, width = 5, height = 5)
plot_circle_3()
invisible(dev.off())

# 保存拼接后的气泡图
pdf(out_file_bubble, height = 3.5, width = 6.5)
print(combined_bubble)
invisible(dev.off())
