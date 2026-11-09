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
out_dir <- Sys.getenv("MG_EXTERNAL_035")
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
mouse_wm <- qread(Sys.getenv("MG_EXTERNAL_019"))
# 筛选特定的 cluster
mouse_wm_mg <- subset(mouse_wm, subset = integrated_snn_res.0.2 %in% c('1','5','17'))
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
  y = human_mg,
  add.cell.ids = c("Mouse", "Human") 
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