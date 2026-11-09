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


# Supplementary Figure 5

# S5A  Iba1/Trim2/DAPI immunofluorescence, needle-like cells and manually marked arrows,20um
# Kind: experimental_or_manual
mg_register_panel("S5A", NULL)

# S5B  Iba1/Trim2/DAPI overview and zoom with 500um/200um scalebars
# Kind: experimental_or_manual
mg_register_panel("S5B", NULL)

# S5C  OLM proportion in five observations
# Kind: experimental_quantification
mg_register_panel("S5C", list(
  list(source="source/recovered/experimental/S3C_OLM_proportion.R", ranges="4-49", requested_ranges="1-49", role="Recovered OLM proportion plot", code=quote({
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[[1L]] else file.path(Sys.getenv("MG_PLOT_INPUT_ROOT", "data/plot_inputs"), "2026-7-29-统计图.xlsx")
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "S3C")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
raw_v <- read_excel(input_file, sheet = "Sheet2", range = "A1:F11", col_names = FALSE)
value <- suppressWarnings(as.numeric(unlist(raw_v[, 6])))
olm <- data.frame(group = "OLM", value = value[!is.na(value)])
stopifnot(nrow(olm) == 5L)
p_olm <- ggplot(olm, aes(x = group, y = value)) +
  geom_boxplot(
    width = 0.38, outlier.shape = NA,
    fill = "#1786C0", alpha = 0.35,
    color = "#1786C0", linewidth = 0.9,
    staplewidth = 0.55
  ) +
  geom_point(
    shape = 21, size = 2.6, stroke = 0.9,
    fill = "#1786C0", color = "black",
    position = position_jitter(width = 0.075, height = 0, seed = 20260729)
  ) +
  scale_y_continuous(
    name = expression(Trim2^"+" * Iba1^"+" * OLM/Iba1^"+" ~ MG),
    breaks = seq(0, 0.12, 0.02),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL) +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(linewidth = 0.9),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black"),
    plot.margin = margin(8, 8, 7, 7)
  ) +
  coord_cartesian(ylim = c(0, 0.12))
ggsave(file.path(output_dir, "S3C_OLM_proportion.png"), p_olm, width = 3.0, height = 4.3, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "S3C_OLM_proportion.pdf"), p_olm, width = 3.0, height = 4.3, device = cairo_pdf)
write.csv(olm, file.path(output_dir, "S3C_plot_values.csv"), row.names = FALSE)
cat("S3C values:\n")
print(olm, digits = 17)
cat("\nSession:\n")
print(sessionInfo())
  }))
))

# S5D  DLB CellChat circle, major Chia source annotations and MG0-MG5
# Kind: computational
mg_register_panel("S5D", list(
  list(source="source/project/Reference/combine/scripts/CellChat/scripts/03.con_DLB_cellchat.r", ranges="1-98", requested_ranges="1-98", role="", code=quote({
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
seu <- qread('data/processed_data/Chia_DLB.qs')

unique(seu$Sample)
# [1] DLB_2 DLB_4 DLB_5 DLB_6 DLB_7 DLB_8
# Levels: ADD_4 ADD_5 ADD_6 ADD_7 CTRL_4 CTRL_5 CTRL_6 CTRL_7 DLB_2 DLB_4 DLB_5 DLB_6 DLB_7 DLB_8 PDD_1 PDD_2 PDD_4 PDD_5 PDD_7 PDD_9
# seu <- subset(seu,subset = Sample %in% c('DLB_2'))

seu.data.input <- GetAssayData(seu, assay = "SCT", slot = "data")
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

cellchat <- computeCommunProb(cellchat, raw.use = FALSE
                              # ,population.size = T
                              )
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
qsave(cellchat,'data/processed_data/con_DLB_cellchat.qs')

cellchat <- qread('data/processed_data/con_DLB_cellchat.qs')
# netVisual_circle(cellchat@net$weight, vertex.weight = 13, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cell_colors <- c(
  "MG0" = "#4DAF4A", 
  "MG1" = "#A65628",
  "MG2" = "#377EB8", 
  "MG3" = "#984EA3", 
  "MG4" = "#E6AB02", 
  "MG5" = "#8B0000", 
  "Excitatory neuron"         = "#E31A1C",
  "Inhibitory neuron"         = "#E32636",
  "Astrocyte"      = "#FF7F00",
  "Oligodendrocyte"           = "#1F78B4",
  "OPC"            = "#A6CEE3",
  "Ependymal cell"      = "#6A3D9A",
  "Endothelial cell"       = "#FB9A99"
)


# 1. 获取当前通讯矩阵中所有细胞类型的准确名称和顺序
group.names <- rownames(cellchat@net$weight)

# 2. 根据该顺序，从你的颜色字典中提取并重新排列颜色
# 确保所有 group.names 都在你的 cell_colors 字典中有定义，否则会产生 NA
ordered_colors <- cell_colors[group.names]

pdf('results/figure/DLB_cellchat.pdf',width = 5,height = 5)
netVisual_circle(
  cellchat@net$weight, 
  vertex.weight = 13, 
  weight.scale = T, 
  label.edge = F, 
  title.name = "Interaction weights/strength",
  sources.use = c("MG0","MG1","MG2","MG3","MG4","MG5"),
  color.use = ordered_colors  # 注意这里使用的是重新对齐顺序的 ordered_colors
)
dev.off()
  }))
))

# S5E  DLB CellChat selected ligand-receptor bubbles, OLM-major cells/MG-OLIG/OLIG-MG
# Kind: computational
mg_register_panel("S5E", list(
  list(source="source/project/Reference/combine/scripts/CellChat/scripts/multi_DLB_cellchat_pop.r", ranges="1-278", requested_ranges="1-278", role="", code=quote({
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

# 加载全局数据---------------------
seu <- qread('data/processed_data/Chia_DLB.qs')

# 获取所有以 "DLB" 开头的样本名
all_samples <- unique(seu$Sample)
dlb_samples <- grep("^DLB", all_samples, value = TRUE)
print(paste("找到以下 DLB 样本将进行分析:", paste(dlb_samples, collapse = ", ")))

# 预设颜色字典 (移到循环外部，避免重复定义)
cell_colors <- c(
  "MG0" = "#4DAF4A", 
  "MG1" = "#A65628",
  "MG2" = "#377EB8", 
  "MG3" = "#984EA3", 
  "MG4" = "#E6AB02", 
  "MG5" = "#8B0000", 
  "Excitatory neuron"         = "#E31A1C",
  "Inhibitory neuron"         = "#E32636",
  "Astrocyte"      = "#FF7F00",
  "Oligodendrocyte"            = "#1F78B4",
  "OPC"             = "#A6CEE3",
  "Ependymal cell"      = "#6A3D9A",
  "Endothelial cell"       = "#FB9A99"
)

# 提前设置好受体-配体相互作用数据库
CellChatDB <- CellChatDB.human  
CellChatDB.use <- subsetDB(CellChatDB, search = c("Cell-Cell Contact"))

# 开始循环---------------------
for (sample_id in dlb_samples) {
  
  cat("\n======================================================\n")
  cat("正在处理样本:", sample_id, "\n")
  cat("======================================================\n")
  
  # 使用 tryCatch 捕获潜在错误，防止循环中断
  tryCatch({
    # 1. 提取当前样本的数据
    seu_sub <- subset(seu, subset = Sample == sample_id)
    
    seu.data.input <- GetAssayData(seu_sub, assay = "SCT", slot = "data")
    seu.identity <- subset(seu_sub@meta.data, select = "cellchat_type")
    
    # 2. 创建并初始化 CellChat 对象
    cellchat <- createCellChat(object = seu.data.input, meta = seu.identity, group.by = "cellchat_type")
    cellchat <- addMeta(cellchat, meta = seu.identity, meta.name = "labels")
    cellchat <- setIdent(cellchat, ident.use = "labels") 
    
    cellchat@DB <- CellChatDB.use
    
    # 3. 核心计算流程
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- smoothData(cellchat, adj = PPI.human)
    
    cellchat <- computeCommunProb(cellchat, raw.use = FALSE,population.size = T)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    
    # 4. 保存 CellChat QS 文件 (动态命名)
    qs_file_path <- paste0('data/processed_data/con_', sample_id, '_pop_cellchat.qs')
    qsave(cellchat, qs_file_path)
    cat("已保存 QS 数据至:", qs_file_path, "\n")
    
    # 5. 获取当前通讯矩阵中的细胞名称和对应的颜色
    group.names <- rownames(cellchat@net$weight)
    ordered_colors <- cell_colors[group.names]
    
    # 6. 绘制并保存 Circle Plot (动态命名)
    circle_pdf_path <- paste0('results/figure/', sample_id, '_pop_cellchat_circle.pdf')
    pdf(circle_pdf_path, width = 5, height = 5)
    netVisual_circle(
      cellchat@net$weight, 
      vertex.weight = 13, 
      weight.scale = T, 
      label.edge = F, 
      title.name = paste("Interaction weights/strength -", sample_id),
      sources.use = c("MG0","MG1","MG2","MG3","MG4","MG5"),
      color.use = ordered_colors 
    )
    dev.off()
    cat("已保存 Circle 图至:", circle_pdf_path, "\n")
    
    # 7. 绘制并保存 Bubble Plot (动态命名)
    bubble_pdf_path <- paste0('results/figure/', sample_id, '_pop_cellchat_bubble.pdf')
    # 确保图赋值给 p1，这样 ggsave 才能正确保存它
    p1 <- netVisual_bubble(cellchat, sources.use = c("MG1"), remove.isolate = TRUE, angle.x = 45) +
      ggtitle(paste("MG1 Interactions -", sample_id))
    ggsave(bubble_pdf_path, plot = p1, width = 5, height = 5)
    cat("已保存 Bubble 图至:", bubble_pdf_path, "\n")
    
  }, error = function(e) {
    # 如果发生错误，打印报错信息并继续下一个样本
    cat("处理样本", sample_id, "时发生错误:\n")
    print(e$message)
    cat("跳过此样本，继续处理下一个...\n")
  })
}

cat("\n所有样本处理完毕！\n")





# MODULE : CellChat Interaction Visualization -------------------

## 1 Parameters --------------------------------------------------------------
# 加载所需包 (请确保已安装)
library(qs)
library(CellChat)
library(ggplot2)
library(patchwork)

# 文件路径设置
input_file <- 'data/processed_data/con_DLB_2_pop_cellchat.qs'
out_dir <- 'results/figure/'
out_file_circle_1  <- paste0(out_dir, "DLB_2_cellchat_circle_multi.pdf")
out_file_circle_2  <- paste0(out_dir, "DLB_2_cellchat_circle_olig.pdf")
out_file_circle_3  <- paste0(out_dir, "DLB_2_cellchat_circle_subset.pdf")
out_file_bubble    <- paste0(out_dir, "DLB_2_bubble_cellchat.pdf")

# 定义细胞群组 (方便后续随意更改)
clusters_source_multi <- c("MG0", "MG1", "MG2", "MG3", "MG4", "MG5")
cluster_olig          <- c("Oligodendrocyte")
cluster_source_single <- c("MG0")
clusters_target_main  <- c("Excitatory neuron","Inhibitory neuron","Astrocyte","Oligodendrocyte","OPC","Ependymal cell","Endothelial cell")

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
  "Excitatory neuron"         = "#E31A1C",
  "Inhibitory neuron"         = "#E32636",
  "Astrocyte"      = "#FF7F00",
  "Oligodendrocyte"            = "#1F78B4",
  "OPC"             = "#A6CEE3",
  "Ependymal cell"      = "#6A3D9A",
  "Endothelial cell"       = "#FB9A99"
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

  }))
))
