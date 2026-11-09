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
    
    cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    
    # 4. 保存 CellChat QS 文件 (动态命名)
    qs_file_path <- paste0('data/processed_data/con_', sample_id, '_cellchat.qs')
    qsave(cellchat, qs_file_path)
    cat("已保存 QS 数据至:", qs_file_path, "\n")
    
    # 5. 获取当前通讯矩阵中的细胞名称和对应的颜色
    group.names <- rownames(cellchat@net$weight)
    ordered_colors <- cell_colors[group.names]
    
    # 6. 绘制并保存 Circle Plot (动态命名)
    circle_pdf_path <- paste0('results/figure/', sample_id, '_cellchat_circle.pdf')
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
    bubble_pdf_path <- paste0('results/figure/', sample_id, '_cellchat_bubble.pdf')
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