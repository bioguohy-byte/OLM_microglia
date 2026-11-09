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
seu <- qread('data/processed_data/mouse_wm_cellchat.qs')
seu.data.input <- GetAssayData(seu, assay = "RNA", slot = "data")
seu.identity <- subset(seu@meta.data, select = "cellchat_type")
cellchat <- createCellChat(object = seu.data.input, meta = seu.identity, group.by = "cellchat_type")

cellchat <- addMeta(cellchat, meta = seu.identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") 
groupSize <- as.numeric(table(cellchat@idents))

# 设置配体-受体相互作用数据库
CellChatDB <- CellChatDB.mouse  
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact" ))
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.mouse)

cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# 保存-----------
qsave(cellchat,'data/processed_data/all_mouse_cellchat.qs')

# all_dlb_cellchat <- qread('data/processed_data/all_mouse_cellchat.qs')
# netVisual_circle(cellchat@net$weight, vertex.weight = 13, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
