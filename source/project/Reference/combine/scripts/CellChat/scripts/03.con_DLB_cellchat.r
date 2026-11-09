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
ggsave('results/figure/DLB_2_cellchat.pdf',p1,width = 5,height = 5)

netVisual_bubble(cellchat, sources.use = c("MG1"), remove.isolate = TRUE,  angle.x = 45)
