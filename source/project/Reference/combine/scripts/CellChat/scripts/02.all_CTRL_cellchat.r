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
seu <- qread('data/processed_data/Chia_CTRL.qs')
seu.data.input <- GetAssayData(seu, assay = "SCT", slot = "data")
seu.identity <- subset(seu@meta.data, select = "cellchat_type")
cellchat <- createCellChat(object = seu.data.input, meta = seu.identity, group.by = "cellchat_type")

cellchat <- addMeta(cellchat, meta = seu.identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") 
groupSize <- as.numeric(table(cellchat@idents))

# 设置配体-受体相互作用数据库
CellChatDB <- CellChatDB.human  
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact" ))
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

weight<-as.data.frame(cellchat@net$weight)
write.table(weight,file="results/table/weight_all_CTRL.xls",sep="\t")
count<-as.data.frame(cellchat@net$count)
write.table(count,file="results/table/count_all_CTRL.xls",sep="\t")
netP<-as.data.frame(cellchat@netP[["prob"]])
write.table(netP,file="results/table/netP_all_CTRL.xls",sep="\t")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# 保存-----------
qsave(cellchat,'data/processed_data/all_CTRL_cellchat.qs')