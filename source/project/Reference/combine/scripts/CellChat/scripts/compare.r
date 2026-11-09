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



DLB<-qread('data/processed_data/con_DLB_cellchat.qs')
DLB <- updateCellChat(DLB)
CTRL<-qread('data/processed_data/con_CTRL_cellchat.qs')
CTRL <- updateCellChat(CTRL)

# 统一细胞类型层级
group.new = levels(CTRL@idents)
DLB <- liftCellChat(DLB, group.new)
object.list <- list(DLB = DLB, CTRL = CTRL)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)


# 可视化差异通讯网络
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
ggsave(filename="results/figure/interaction_circle2.pdf",width = 12, height = 7)

# 对比两个区域的通讯网络
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F,edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# 排名网络对比
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use=c("#B3128A","#74CAFF"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use=c("#B3128A","#74CAFF"))
p3=gg1 + gg2
ggsave(filename="results/figure/flow.pdf",plot=p3, width = 12, height = 8)


levels(CTRL@idents)


# 保存所有通讯信息
g1<-netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13),  comparison = c(1, 2), angle.x = 45,color.text=c("#B3128A","#74CAFF"))
g1data<-as.data.frame(g1$data)
write.table(g1data,file="All_ineration_info.xls",sep="\t")

# 分析信号角色散点图
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
p1=patchwork::wrap_plots(plots = gg)
ggsave(filename="major_scatter.pdf",plot=p1, width = 8, height = 4)

# 分析Microglia(PCDH9high)细胞的通讯变化
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Microglia(PCDH9high)", color.use = c("#7AB000", "#B3128A","#74CAFF"))
gg1data<-as.data.frame(gg1$data)
write.table(gg1data,file="Microglia(PCDH9high)_scatter_info.xls",sep="\t")
ggsave(filename="Microglia(PCDH9high)_scatter.pdf",plot=gg1, width = 6, height = 4)

# 分析特定信号通路的通讯
# 所有细胞向Microglia(PCDH9high)的通讯
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CADM","NEGR"))
netVisual_bubble(cellchat, sources.use = c(1:13), targets.use = 4, pairLR.use = pairLR.use, remove.isolate = TRUE,  comparison = c(1, 2), angle.x = 45,thresh=0.05)
netVisual_bubble(Hippo, sources.use = c(1:13), targets.use = 4, pairLR.use = pairLR.use, remove.isolate = TRUE)

# Microglia(PCDH9high)向其他细胞的通讯
netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:12), pairLR.use = pairLR.use, remove.isolate = TRUE,  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(Hippo, sources.use = 6, targets.use = c(1:12), pairLR.use = pairLR.use, remove.isolate = TRUE)
