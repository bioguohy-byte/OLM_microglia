library(Seurat)
library(qs)

library(future) # 新增
plan("multicore", workers = 20) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限

# 1. 路径设置与数据读取
main_path <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Su_et_al./data/processed_data/Su_CSC_alldata_v2.qs")
sub_path  <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs")

seu_main <- qread(main_path)
seu_sub  <- qread(sub_path)

# 2. 配置变量
target_col    <- "CellType"     # 主对象原始列名
sub_label_col <- "CellType1"    # 子对象标签列名
result_col    <- "cellchat_type" # 最终用于 CellChat 的列名

# 3. 执行映射逻辑
# A. 初始化：先复制主对象的标签
seu_main[[result_col]] <- as.character(seu_main[[target_col, drop = TRUE]])

# B. 建立映射向量并替换共有细胞标签
mapping_vec <- setNames(as.character(seu_sub[[sub_label_col, drop = TRUE]]), colnames(seu_sub))
common_cells <- intersect(names(mapping_vec), colnames(seu_main))
seu_main@meta.data[common_cells, result_col] <- mapping_vec[common_cells]

# 4. 标准化与保存
su_cellchat <- NormalizeData(seu_main)

# 验证映射结果
message("Final Identities Checklist:")
print(table(su_cellchat[[result_col, drop = TRUE]]))

# 保存处理好的 Seurat 对象
qsave(su_cellchat, file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/su_cellchat_ready.qs"))


library(CellChat)
library(Seurat)
library(qs)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap) # 必须加载，用于处理热图保存

# 1. 环境准备
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/"))
dir.create("results/figure/", recursive = TRUE, showWarnings = FALSE)
dir.create("data/saved_data/", recursive = TRUE, showWarnings = FALSE)

# 读取预处理好的数据
su_cellchat <- qread("data/processed_data/su_cellchat_ready.qs")

# 2. 初始化 CellChat 对象
Idents(su_cellchat) <- "cellchat_type"
# 解决 samples 列警告问题
if(!"samples" %in% colnames(su_cellchat@meta.data)){
  su_cellchat$samples <- su_cellchat$Sample 
}

data.input <- GetAssayData(su_cellchat, assay = "RNA", slot = "data")
meta <- su_cellchat@meta.data

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellchat_type")

# 3. 设置数据库：仅 Cell-Cell Contact
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
cellchat@DB <- CellChatDB.use

# 4. 核心计算逻辑
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 使用 PPI.human 进行平滑处理
cellchat <- smoothData(cellchat, adj = PPI.human)

# 5. 通信概率与网络聚合
cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# --- 6. 先执行对象保存（防崩溃保险） ---
qsave(cellchat, "data/saved_data/su_cellchat_contact_ppi.qs")
message(">>> 数据已安全保存至 data/saved_data/")

# --- 7. 图片自动化导出（修复版） ---

# A. 整体互作 Circle Plot (Base R)
pdf("results/figure/Contact_Overall_Network.pdf", width = 12, height = 6)
par(mfrow = c(1,2), xpd=TRUE)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Number of Interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Interaction Strength")
dev.off()

# B. 信号角色热图 (修复后的保存方式)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 10, height = 14)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 10, height = 14)

# 使用 PDF 设备保存 ComplexHeatmap 列表
pdf("results/figure/Contact_Role_Heatmap.pdf", width = 14, height = 12)
draw(ht1 + ht2) # 使用 draw 函数将两个热图拼在一起
dev.off()

# C. L-R 气泡图 (ggplot2 对象，继续使用 ggsave)
pb <- netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE)
# 因为群组较多，建议增加保存宽度
ggsave("results/figure/Contact_L-R_Bubble.pdf", plot = pb, width = 12, height = 16)

message(">>> 所有图片已存入 results/figure/，流程顺利完成！")


cellchat <- qread('data/saved_data/su_cellchat_contact_ppi.qs')
p1 <- netVisual_circle(cellchat@net$weight, vertex.weight = 20, weight.scale = T, sources.use = c("MG0","MG1","MG2","MG3","MG4","MG5"), targets.use = c("OPC","OLIG","GABA.N","AST","Glut.N","Ependymal","Endothelial","Choroid.plexus","Cajal-Retzius"),
                       label.edge = F, title.name = "Number of Interactions")

p2 <- netVisual_circle(cellchat@net$weight, vertex.weight = 20, weight.scale = T, sources.use = c("OPC","OLIG","GABA.N","AST","Glut.N","Ependymal","Endothelial","Choroid.plexus","Cajal-Retzius"), targets.use = c("MG0","MG1","MG2","MG3","MG4","MG5"),
                       label.edge = F, title.name = "Number of Interactions")

pairLR.use_1 <- data.frame(interaction_name = c("NRXN3_NLGN1","NRXN3_LRRTM4","CD6_ALCAM","CNTN2_L1CAM"))
netVisual_bubble(cellchat,sources.use = c("MG0","MG1","MG2","MG3","MG4","MG5"), targets.use = c("OPC","OLIG","GABA.N","AST","Glut.N","Ependymal","Endothelial","Choroid.plexus","Cajal-Retzius"), pairLR.use = pairLR.use_1,  angle.x = 45,thresh=0.05)
# filter p-value> 0.05
netVisual_bubble(Hippo, sources.use = c(1:12), targets.use = 6, pairLR.use = pairLR.use, remove.isolate = TRUE)

netVisual_bubble(cellchat, 
                 sources.use = c("MG0"), targets.use = c("OPC","OLIG","GABA.N","AST","Glut.N","Ependymal","Endothelial","Choroid.plexus","Cajal-Retzius"),
                 # pairLR.use = pairLR.use_1,
                 remove.isolate = TRUE, 
                 angle.x = 45)

netVisual_bubble(cellchat, 
                 sources.use = c("OPC","OLIG","GABA.N","AST","Glut.N","Ependymal","Endothelial","Choroid.plexus","Cajal-Retzius"), targets.use = c("MG0"),
                 remove.isolate = TRUE, 
                 angle.x = 45)
