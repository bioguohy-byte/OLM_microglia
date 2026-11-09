rm(list=ls())
# ==============================================================================
# Seurat V4 标准无偏差整合 (No Reference) - Human & Mouse Only
# ==============================================================================
# 设置内存限制为 100GB
options(future.globals.maxSize = 100 * 1024^3) 

library(dplyr)
library(Seurat) # V4
library(patchwork)
library(ggplot2)
library(qs) 
library(biomaRt) 
library(future) 

# # 并行设置 (根据服务器负载开启)
plan("multicore", workers = 20)

message("Current Seurat Version: ", packageVersion("Seurat"))

# 设定输出路径 (建议修改文件夹名以区分包含猴子的版本)
out_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)

# ==============================================================================
# 1. BiomaRt 连接 (Dec 2021)
# ==============================================================================
message("--- Connecting to Ensembl Archive (Dec 2021) ---")
host_url <- "https://dec2021.archive.ensembl.org/"

human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host_url)
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host_url)

message("BiomaRt connections established (Human & Mouse).")

# ==============================================================================
# 2. 定义转换函数 (同源基因转换)
# ==============================================================================
Convert_With_Mart <- function(obj, source_mart, target_mart, species_name) {
  DefaultAssay(obj) <- "RNA"
  counts_matrix <- GetAssayData(obj, assay = "RNA", slot = "counts")
  original_genes <- rownames(counts_matrix)
  
  message(paste0("Querying ", length(original_genes), " genes for ", species_name, "..."))
  
  # 获取同源基因映射
  gene_map <- getLDS(attributes = c("external_gene_name"), 
                     filters = "external_gene_name", 
                     values = original_genes, 
                     mart = source_mart, 
                     attributesL = c("hgnc_symbol"), 
                     martL = target_mart, 
                     uniqueRows = TRUE)
  
  colnames(gene_map) <- c("Original", "Human")
  gene_map <- gene_map[gene_map$Human != "" & !is.na(gene_map$Human), ]
  
  # 筛选有效基因
  valid_genes <- intersect(original_genes, gene_map$Original)
  counts_subset <- counts_matrix[valid_genes, ]
  
  # 处理多对一映射 (去重)
  gene_map <- gene_map[!duplicated(gene_map$Original), ]
  rownames(gene_map) <- gene_map$Original
  
  new_names <- gene_map[rownames(counts_subset), "Human"]
  
  # 如果有重复的 Human Symbol，进行合并
  if(any(duplicated(new_names))) {
    rownames(counts_subset) <- new_names
    counts_final <- as(rowsum(as.matrix(counts_subset), group = new_names), "dgCMatrix")
  } else {
    rownames(counts_subset) <- new_names
    counts_final <- counts_subset
  }
  
  # 创建新对象
  new_obj <- CreateSeuratObject(counts = counts_final, meta.data = obj@meta.data)
  return(new_obj)
}

# ==============================================================================
# 3. 加载数据与转换
# ==============================================================================
message("--- Loading and Converting Data ---")

# --- Mouse ---
message("Processing Mouse Data...")
mouse_wm <- qread(Sys.getenv("MG_EXTERNAL_017"))
# 筛选特定的 cluster
mouse_wm_mg <- subset(mouse_wm, subset = CellType1 %in% "Microglia")
rm(mouse_wm); gc()

# 转换为 Human Gene Symbol
mouse_wm_mg <- Convert_With_Mart(mouse_wm_mg, mouse_mart, human_mart, "Mouse")
mouse_wm_mg$Species <- "Mouse"
mouse_wm_mg$integration_batch <- "Mouse"

# --- Human ---
message("Processing Human Data...")
human_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))
human_mg$Species <- "Human"
# 假设 Human 数据集内部包含多个 Study，这里将其作为批次处理
human_mg$integration_batch <- paste0("Human_", make.names(human_mg$Study)) 
DefaultAssay(human_mg) <- "RNA"

###################
message("--- Checking gene numbers before intersection ---")
cat("Mouse genes:", nrow(mouse_wm_mg), "\n")
cat("Human genes:", nrow(human_mg), "\n")

# =============================
# 关键步骤：统一 gene space (仅 Human & Mouse)
# =============================
message("--- Taking shared gene intersection ---")

shared_genes <- intersect(rownames(mouse_wm_mg), rownames(human_mg))

cat("Shared genes:", length(shared_genes), "\n")

# 质量检查
if(length(shared_genes) < 8000){
  warning("Shared genes < 8000. Check ortholog mapping quality!")
}

# 强制二者使用相同 gene space
mouse_wm_mg <- subset(mouse_wm_mg, features = shared_genes)
human_mg    <- subset(human_mg, features = shared_genes)

# ==============================================================================
# 4. Merge & Split
# ==============================================================================
message("--- Merging datasets ---")
# 仅合并 Mouse 和 Human
combined_raw <- merge(
  x = mouse_wm_mg,
  y = human_mg
)
rm(mouse_wm_mg, human_mg); gc()

message("--- Splitting object by integration_batch ---")
# 这里会将 Mouse (1个batch) 和 Human (可能有多个Study batch) 分开
obj.list <- SplitObject(combined_raw, split.by = "integration_batch")
rm(combined_raw); gc()

# ==============================================================================
# 5. Pre-processing & Find Anchors (标准模式)
# ==============================================================================
message("--- Normalizing and finding features ---")
obj.list <- lapply(X = obj.list, FUN = function(x) {
  DefaultAssay(x) <- "RNA" 
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

message("--- Selecting integration features ---")
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)

# 保存一下选定的 features，防止后续复现需要
qsave(features, file = file.path(out_dir, "00.integration_features.qs"))

message("--- Finding Integration Anchors (Method = CCA, All-to-All) ---")
# 寻找 Mouse 和 Human 各个 Batch 之间的 Anchors
anchors <- FindIntegrationAnchors(object.list = obj.list, 
                                  anchor.features = features) 

# 保存中间过程：Anchors (非常重要，耗时最长的一步)
message("--- Saving Anchors ---")
qsave(anchors, file = file.path(out_dir, "01.integration_anchors_noref.qs"))

# ==============================================================================
# 6. Clean Anchors (补丁: 移除不必要的 Assay 节省内存)
# ==============================================================================
message("--- Cleaning residual SCT/Integrated assays from anchor objects ---")
anchors@object.list <- lapply(anchors@object.list, function(x) {
  DefaultAssay(x) <- "RNA"
  if ("SCT" %in% Assays(x)) x[["SCT"]] <- NULL
  if ("integrated" %in% Assays(x)) x[["integrated"]] <- NULL
  return(x)
})

# ==============================================================================
# 7. Integrate & Process
# ==============================================================================
message("--- Integrating Data ---")
integrated_mg <- IntegrateData(anchorset = anchors)

# 清理内存
rm(anchors, obj.list); gc()

message("--- Dimension Reduction ---")
DefaultAssay(integrated_mg) <- "integrated"
integrated_mg <- ScaleData(integrated_mg, verbose = FALSE)
integrated_mg <- RunPCA(integrated_mg, npcs = 30, verbose = FALSE)
integrated_mg <- RunUMAP(integrated_mg, reduction = "pca", dims = 1:5)
integrated_mg <- FindNeighbors(integrated_mg, reduction = "pca", dims = 1:5)
integrated_mg <- FindClusters(integrated_mg, resolution = 0.5)

# ==============================================================================
# 8. Save Final Object
# ==============================================================================
message("--- Saving Final Object ---")
qsave(integrated_mg, file = file.path(out_dir, "02.cross_species_integrated_noref.qs"))

# ==============================================================================
# 9. Plotting (放在最后)
# ==============================================================================
message("--- Plotting ---")

# 图1：按物种着色
p1 <- DimPlot(integrated_mg, reduction = "umap", group.by = "Species", raster=FALSE) + 
  ggtitle("By Species") + 
  theme(plot.title = element_text(hjust = 0.5))

# 图2：按整合批次着色 (Mouse vs Human_Study1 vs Human_Study2 ...)
p2 <- DimPlot(integrated_mg, reduction = "umap", group.by = "integration_batch", raster=FALSE, label = TRUE) + 
  NoLegend() + 
  ggtitle("By Batch") + 
  theme(plot.title = element_text(hjust = 0.5))

# 图3：按聚类着色
p3 <- DimPlot(integrated_mg, reduction = "umap", group.by = "integrated_snn_res.0.5", label = TRUE) + 
  ggtitle("By Cluster (res 0.5)")

# 组合图
combined_plot <- (p1 + p2) / p3

ggsave(file.path(out_dir, "integration_noref_check.pdf"), combined_plot, width = 14, height = 12)

message("Finished! Human-Mouse integration complete.")



#-------------------
tmp <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/02.cross_species_integrated_noref.qs"))
human_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))


human_celltype <- human_mg[["CellType1"]]
tmp <- AddMetaData(object = tmp, metadata = human_celltype)
DimPlot(tmp,group.by = "CellType1")
DimPlot(tmp,group.by = "integrated_snn_res.0.5",label = T)

DefaultAssay(tmp) <- 'integrated'
tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:5)
tmp <- FindClusters(tmp, resolution = 0.1)
DimPlot(tmp,group.by = "integrated_snn_res.0.1")

tmp$CellType2 <- as.character(tmp$integrated_snn_res.0.1)
tmp$CellType2[tmp$CellType2 %in% c('0')] <- 'cluster0'
tmp$CellType2[tmp$CellType2 %in% c('1')] <- 'cluster1'
tmp$CellType2[tmp$CellType2 %in% c('2')] <- 'cluster2'
tmp$CellType2[tmp$CellType2 %in% c('3')] <- 'cluster3'

DimPlot(tmp,group.by = 'CellType2',split.by = "CellType2")


qsave(tmp,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/03.human_mouse_celldefine.qs"))



# 跨物种相关性---------
library(dplyr)
library(Seurat) 
library(patchwork)
library(ggplot2)
library(qs) 
library(biomaRt) 
library(future) 


tmp <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/03.human_mouse_celldefine.qs"))

tmp$CellType1 <- as.character(tmp$CellType1)
tmp$CellType3 <- ifelse(is.na(tmp$CellType1), 
                        tmp$CellType2, 
                        tmp$CellType1)
qsave(tmp,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/03.human_mouse_celldefine_v2.qs"))


cluster.averages <- AverageExpression(tmp, 
                                      assays = "RNA", 
                                      group.by = "CellType3") 

mouse_clusters <- c('cluster0','cluster1','cluster2','cluster3')
human_clusters <- c("MG0","MG1","MG2","MG3","MG4","MG5")

# 2. 计算Spearman 相关性矩阵
ht_df <- cluster.averages[["RNA"]] %>% cor()

# 3. 绘制热图
pheatmap::pheatmap(ht_df[mouse_clusters, human_clusters],
                   show_rownames = T, show_colnames = T,
                   cluster_rows = F, cluster_cols = F,
                   use_raster = F, border_color = "white",
                   display_numbers = TRUE, number_color = "black",
                   fontsize_number = 15, fontsize = 15, legend = F,
                   color = colorRampPalette(c("#72bcd5", "#F7F5F2", "#FD5D5D"))(30))


# 把CellType2映射回小鼠原来的单细胞数据中------------
mouse_wm <- qread(Sys.getenv("MG_EXTERNAL_017"))
mouse_wm_mg <- subset(mouse_wm, subset = CellType1 %in% "Microglia")

tmp <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/03.human_mouse_celldefine.qs"))

mouse_wm_mg$Mapped_CellType <- tmp$CellType2[match(colnames(mouse_wm_mg), colnames(tmp))]

qsave(mouse_wm_mg,Sys.getenv("MG_EXTERNAL_021"))



unique(mouse_wm_mg$Mapped_CellType)
DimPlot(mouse_wm_mg,group.by = 'Mapped_CellType')

Idents(mouse_wm_mg) <- "Mapped_CellType"
mouse_mg_markers <- FindAllMarkers(mouse_wm_mg, 
                                   only.pos = TRUE, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.25)

# 1. 读取所有的 Marker 基因数据

# 2. 设定目标 Cluster 并提取符合阈值的基因
target_cluster <- "cluster2"
cluster_color <- "#D62728" # 你定义的 cluster 专属红色

cluster_markers <- mouse_mg_markers %>%
  dplyr::filter(cluster == target_cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(avg_log2FC > 0.5) %>%
  arrange(desc(avg_log2FC)) # 按 log2FC 降序排列

# 3. 基因 ID 转换 (SYMBOL 转 ENTREZID)
# 注意：这里默认使用人的数据库 org.Hs.eg.db，如果是小鼠请改为 org.Mm.eg.db
g1 <- cluster_markers$gene
gg <- bitr(g1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

# 4. 运行单一群的 GO 富集分析 (BP)
ego_cluster <- enrichGO(
  gene = gg$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# 5. 提取富集结果并准备 Top 10 作图数据
plot_data <- as.data.frame(ego_cluster) %>%
  arrange(p.adjust) %>%  # p.adjust 越小越显著
  head(10) %>%           # 只取 Top 10
  mutate(logP = -log10(p.adjust)) # 计算用于绘图的纵坐标

# 设置因子水平，确保条形图按照显著性从上到下排列
plot_data$Description <- factor(plot_data$Description, levels = rev(plot_data$Description))

# 6. 绘制 cluster 专属融合风格柱状图
p_cluster <- ggplot(plot_data, aes(x = logP, y = Description)) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7, fill = cluster_color) +
  # 将 GO 条目名称直接写在条形图起始位置
  geom_text(aes(x = 0.05, label = Description), size = 4.5, hjust = 0, color = "black") +
  theme_minimal() +
  labs(
    title = paste0("Top 10 GO Terms (BP): Cluster ", target_cluster),
    x = "-log10(p.adjust)",
    y = ""
  ) +
  theme(
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"), 
    axis.text.y = element_blank(), # 隐藏原有的 Y 轴刻度文字
    axis.ticks.x = element_line(color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12)
  ) +
  # 留出右侧 expansion 空间，防止过长的文字标签被切掉
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) 

# 7. 保存结果
dir.create("results/figure/cluster_go/", recursive = TRUE, showWarnings = FALSE)

file_name_base <- paste0('results/figure/cluster_go/Cluster_', target_cluster, '_top10_GO')
ggsave(paste0(file_name_base, ".pdf"), plot = p_cluster, width = 8, height = 5)
ggsave(paste0(file_name_base, ".png"), plot = p_cluster, width = 8, height = 5, dpi = 300)

cat(target_cluster, "的 GO 分析与作图已全部完成！\n")


# MODULE : Calculate and Visualize Module Score for Cluster Marker Genes -----
rm(list=ls())
library(Seurat)
library(ggplot2)
library(qs)
library(dplyr)

## 1 Parameters --------------------------------------------------------------
# 1.1 IO Paths
seurat_path <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/03.human_mouse_celldefine_v2.qs")
marker_path <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/saved_data/CellType1_markers.qs")
output_dir  <- 'results/figure/'

# 1.2 Analysis Variables
assay_name         <- 'RNA'
target_cluster     <- "MG0"
top_n              <- 50
p_val_adj_thresh   <- 0.05      
avg_log2FC_thresh  <- 0.5       


## 2 Process Data ------------------------------------------------------------
# 2.1 Load Data
seurat_obj <- qread(seurat_path)
seurat_obj <- subset(seurat_obj,subset = Species %in% c("Mouse"))

DefaultAssay(seurat_obj) <- assay_name
all_markers <- qread(marker_path)

# 2.2 Filter and Extract Target Genes
target_top_genes <- all_markers %>%
  dplyr::filter(cluster == target_cluster) %>%
  dplyr::filter(p_val_adj < p_val_adj_thresh, avg_log2FC > avg_log2FC_thresh) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::pull(gene) 

if(length(target_top_genes) == 0) {
  stop(sprintf("未找到符合条件的 Marker 基因，请检查过滤条件 (P.adj < %s, log2FC > %s)。", 
               p_val_adj_thresh, avg_log2FC_thresh))
}

# 2.3 Calculate Module Score
score_base_name <- paste0(target_cluster, "_Top", top_n, "_Score")
score_col_name  <- paste0(score_base_name, "1") 

seurat_obj <- AddModuleScore(
  object   = seurat_obj,
  features = list(target_top_genes),
  name     = score_base_name
)


## 3 Generate Plots ----------------------------------------------------------
plot_title <- paste0(target_cluster, " Top ", top_n, " Signature Score")

p_score <- FeaturePlot(
  seurat_obj, 
  features = score_col_name, 
  order    = TRUE,          
  cols     = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(10), 
  raster   = FALSE
) +
  theme_classic() +
  theme(
    axis.line        = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    axis.title       = element_blank(),
    legend.position  = "right",              
    legend.text      = element_text(size = 8),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
    strip.text       = element_text(face = "bold", size = 12)
  ) +
  labs(title = plot_title) + 
  NoAxes()

print(p_score)


## 4 Save Results ------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_formats <- c("pdf", "png")

for (fmt in save_formats) {
  save_name <- sprintf('ModuleScore_Split_%s_Top%d.%s', target_cluster, top_n, fmt)
  save_path <- file.path(output_dir, save_name)
  
  ggsave(
    filename = save_path, 
    plot     = p_score, 
    height   = 5, 
    width    = 6, 
    bg       = '#ffffff'
  )
}
message("Module Done.")
