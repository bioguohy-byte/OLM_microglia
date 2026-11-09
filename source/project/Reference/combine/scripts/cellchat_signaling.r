library(CellChat)
library(Seurat)
library(qs)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(future)

# 0. 并行计算设置 (根据服务器核心数调整)
plan("multicore", workers = 10) 
options(future.globals.maxSize = 100 * 1024^3) # 设置 100GB 内存上限

# 1. 环境与路径设置
main_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/")
setwd(main_dir)

dir.create("results/figure/", recursive = TRUE, showWarnings = FALSE)
dir.create("data/saved_data/", recursive = TRUE, showWarnings = FALSE)

# 读取预处理好的数据
su_cellchat <- qread("data/processed_data/su_cellchat_ready.qs")

# 2. 初始化 CellChat 对象
Idents(su_cellchat) <- "cellchat_type"
if(!"samples" %in% colnames(su_cellchat@meta.data)){
  su_cellchat$samples <- su_cellchat$Sample 
}

data.input <- GetAssayData(su_cellchat, assay = "RNA", slot = "data")
meta <- su_cellchat@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellchat_type")

# 3. 设置数据库：使用全库 (包含 Secreted, ECM, Cell-Cell Contact)
# 这一步能确保你关心的 NRXN3 和 CNTN2 不会被漏掉
CellChatDB <- CellChatDB.human 
cellchat@DB <- CellChatDB 

# 4. 核心计算流程
cellchat <- subsetData(cellchat) # 提取数据库相关的基因
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 使用 PPI.human 平滑处理 (有助于找回稀疏数据中的信号)
cellchat <- smoothData(cellchat, adj = PPI.human)

# 计算通信概率 (raw.use = FALSE 以应用 PPI 平滑结果)
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, nboot = 100) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# --- 5. 优先保存对象 (安全保险) ---
# 使用 Full_DB 命名以区分之前的版本
save_path <- "data/saved_data/su_cellchat_Full_DB_ppi.qs"
qsave(cellchat, save_path)
message(">>> [SUCCESS] 全库计算完成，对象已保存至: ", save_path)

# --- 6. 自动化图片导出 ---

# A. 整体互作网络图
pdf("results/figure/Full_Overall_Network.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd=TRUE)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Number of Interactions (Full DB)")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Interaction Strength (Full DB)")
dev.off()

# B. 信号角色热图 ( outgoing & incoming )
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 12, height = 16)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 12, height = 16)

pdf("results/figure/Full_Role_Heatmap.pdf", width = 16, height = 14)
draw(ht1 + ht2)
dev.off()

# C. 基础气泡图 (展示所有显著的 L-R pairs)
pb <- netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE)
ggsave("results/figure/Full_L-R_Bubble_All.pdf", plot = pb, width = 14, height = 25, limitsize = FALSE)

# --- 7. 针对 MG0 的定向分析绘图 ---

cellchat <- qread('data/saved_data/su_cellchat_Full_DB_ppi.qs')

# MG0 作为发送者到其他所有细胞
p_mg0_out <- netVisual_bubble(cellchat, 
                              sources.use = "MG0", 
                              targets.use = c("OPC","OLIG","GABA.N","AST","Glut.N","Ependymal","Endothelial","Choroid.plexus","Cajal-Retzius"),
                              remove.isolate = TRUE, angle.x = 45) + 
  ggtitle("MG0 Outgoing Signaling (Full DB)")

# 其他细胞作为发送者到 MG0
p_mg0_in <- netVisual_bubble(cellchat, 
                             sources.use = c("OPC","OLIG","GABA.N","AST","Glut.N","Ependymal","Endothelial","Choroid.plexus","Cajal-Retzius"),
                             targets.use = "MG0",
                             remove.isolate = TRUE, angle.x = 45) + 
  ggtitle("MG0 Incoming Signaling (Full DB)")

ggsave("results/figure/MG0_Specific_Bubble.pdf", plot = (p_mg0_out / p_mg0_in), width = 12, height = 20)

message(">>> [SUCCESS] 所有全库分析图片已保存完成！")