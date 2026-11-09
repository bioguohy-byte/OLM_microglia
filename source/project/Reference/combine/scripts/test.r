rm(list=ls())
options(future.globals.maxSize = 100 * 1024^3) # 设置内存限制为100GB，防止大数据集整合时报错

# ==============================================================================
# 0. Library Loading
# ==============================================================================
library(dplyr)
library(Seurat) # 确保是 V4 版本
library(patchwork)
library(ggplot2)
library(qs) 
library(homologene) # 用于同源基因转换

library(future) # 新增
plan("multicore", workers = 20) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限


# 检查 Seurat 版本，确保心里有数
message("Current Seurat Version: ", packageVersion("Seurat"))

# 设定输出路径
out_dir <- Sys.getenv("MG_EXTERNAL_030")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)

# # ==============================================================================
# # 1. Load Data & Subsetting
# # ==============================================================================
# message("--- Loading and subsetting data ---")
# 
# # 1.1 Mouse Data (TaxID: 10090)
# mouse_wm <- qread(Sys.getenv("MG_EXTERNAL_019"))
# mouse_wm_mg <- subset(mouse_wm, subset = integrated_snn_res.0.2 %in% c('1','5','17'))
# mouse_wm_mg$Species <- "Mouse"
# mouse_wm_mg$integration_batch <- "Mouse" 
# rm(mouse_wm); gc()
# 
# # 1.2 Monkey Data (TaxID: 9544 Rhesus / 9541 Cynomolgus)
# monkey_hip <- qread(Sys.getenv("MG_EXTERNAL_020"))
# monkey_hip_mg <- subset(monkey_hip, subset = integrated_snn_res.0.1 %in% c('0','9','11'))
# monkey_hip_mg$Species <- "Monkey"
# monkey_hip_mg$integration_batch <- "Monkey"
# rm(monkey_hip); gc()
# 
# # 1.3 Human Data (TaxID: 9606)
# human_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))
# human_mg$Species <- "Human"
# # 确保 Study 列作为批次名时是安全的字符
# human_mg$integration_batch <- paste0("Human_", make.names(human_mg$Study)) 
# 
# # ==============================================================================
# # 2. Gene Homology Handling (Seurat V4 Style)
# # ==============================================================================
# message("--- Converting gene names using Homologene for Seurat V4 ---")
# 
# # 定义 V4 专用的转换函数
# Convert_Using_Homologene_V4 <- function(obj, species_taxid) {
#   # [Seurat V4 Standard] 使用 GetAssayData 提取 count 矩阵
#   counts_matrix <- GetAssayData(obj, assay = "RNA", slot = "counts")
#   original_genes <- rownames(counts_matrix)
#   
#   message(paste0("Querying homologs for TaxID: ", species_taxid, "..."))
#   
#   # 使用 homologene 数据库查询 (9606 = Human)
#   homolog_table <- homologene(original_genes, inTax = species_taxid, outTax = 9606)
#   
#   if(nrow(homolog_table) == 0) {
#     stop("No homologs found! Check network connection or TaxID.")
#   }
#   
#   # 提取输入ID (如 '10090') 和 输出ID ('9606') 的列名
#   input_col <- as.character(species_taxid)
#   output_col <- "9606"
#   
#   # 1. 筛选矩阵：只保留有同源信息的基因
#   # 注意：homologene 可能会把多个 mouse 基因映射到一个 human 基因，也可能反之
#   # 我们这里以 homologene 的结果为准
#   valid_genes <- intersect(original_genes, homolog_table[[input_col]])
#   counts_subset <- counts_matrix[valid_genes, ]
#   
#   # 2. 确保顺序一致
#   match_idx <- match(rownames(counts_subset), homolog_table[[input_col]])
#   ordered_homologs <- homolog_table[match_idx, ]
#   
#   # 3. 替换基因名
#   new_gene_names <- ordered_homologs[[output_col]]
#   
#   # 4. 处理多对一映射 (Duplicate Gene Names)
#   # 在 Seurat V4 中，创建对象不允许有重复行名，必须聚合
#   if(any(duplicated(new_gene_names))) {
#     message(paste0("  Found ", sum(duplicated(new_gene_names)), " duplicated human gene symbols. Aggregating expression..."))
#     
#     # 将稀疏矩阵转为普通矩阵进行 rowsum (虽然耗内存，但是最稳妥的内置方法)
#     # 如果内存不够，可以使用 Matrix.utils::aggregate.Matrix (但在新版R中常有问题)
#     # 这里使用基础 R 的 rowsum 配合 sparse matrix 转换
#     
#     # 临时给行名赋值（允许重复）
#     rownames(counts_subset) <- new_gene_names
#     
#     # 按行名求和 (这一步是关键，解决多对一问题)
#     # as.matrix 可能会爆内存，如果内存足够大则没问题。
#     # 优化写法：直接利用 Matrix 包特性，避免转 dense matrix
#     # 但为了代码稳健性，只要细胞数不是百万级，rowsum(as.matrix) 是可以的。
#     # 既然是 subset 之后的小胶质细胞，数据量应该可控。
#     counts_aggregated <- rowsum(as.matrix(counts_subset), group = new_gene_names)
#     
#     counts_final <- as(counts_aggregated, "dgCMatrix")
#   } else {
#     rownames(counts_subset) <- new_gene_names
#     counts_final <- counts_subset
#   }
#   
#   message(paste0("  Converted: ", length(original_genes), " -> ", nrow(counts_final), " Human Genes"))
#   
#   # [Seurat V4 Standard] 重建对象
#   new_obj <- CreateSeuratObject(counts = counts_final, meta.data = obj@meta.data)
#   return(new_obj)
# }
# 
# # --- 执行转换 ---
# # Mouse (10090)
# mouse_wm_mg <- Convert_Using_Homologene_V4(mouse_wm_mg, species_taxid = 10090)
# 
# # Monkey (9544 - Rhesus)
# monkey_hip_mg <- Convert_Using_Homologene_V4(monkey_hip_mg, species_taxid = 9544)
# 
# # Human 不需要转换，但建议确保也是唯一 gene names
# # 如果 Human 数据集里有重复基因名，CreateSeuratObject 会报错，前面 qread 读进来的应该是处理好的
# 
# # ==============================================================================
# # 3. Merge & Split (Seurat V4 Workflow)
# # ==============================================================================
# message("--- Merging datasets ---")
# 
# # 先合并
# combined_raw <- merge(
#   x = mouse_wm_mg,
#   y = list(monkey_hip_mg, human_mg),
#   add.cell.ids = c("Mouse", "Monkey", "Human") 
# )
# 
# # 清理内存
# rm(mouse_wm_mg, monkey_hip_mg, human_mg); gc()
# 
# message("--- Splitting object by integration_batch ---")
# # [Seurat V4 Standard] 必须拆分成 List 才能进行 FindIntegrationAnchors
# obj.list <- SplitObject(combined_raw, split.by = "integration_batch")
# 
# # 释放大对象
# rm(combined_raw); gc()
# 
# # ==============================================================================
# # 4. Integration (Seurat V4: CCA or RPCA)
# # ==============================================================================
# message("--- Pre-processing list for integration ---")
# 
# # 标准预处理
# obj.list <- lapply(X = obj.list, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# 
# message("--- Selecting integration features ---")
# features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)
# 
# message("--- Finding Integration Anchors (Method = CCA) ---")
# # [Seurat V4 Standard] 寻找锚点
# # 如果数据量巨大，可以在这里加 reduction = "rpca"
# anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
# 
# # *** 保存中间结果 (V4 对象) ***
# qsave(anchors, file = file.path(out_dir, "01.integration_anchors_v4.qs"))
# 
# message("--- Integrating Data ---")
# # [Seurat V4 Standard] 执行整合
# integrated_mg <- IntegrateData(anchorset = anchors)
# 
# # 清理
# rm(anchors, obj.list); gc()
# 
# # ==============================================================================
# # 5. Post-Processing
# # ==============================================================================
# message("--- Dimension Reduction ---")
# 
# # [Seurat V4 Standard] 切换到 integrated assay
# DefaultAssay(integrated_mg) <- "integrated"
# 
# integrated_mg <- ScaleData(integrated_mg, verbose = FALSE)
# integrated_mg <- RunPCA(integrated_mg, npcs = 30, verbose = FALSE)
# integrated_mg <- RunUMAP(integrated_mg, reduction = "pca", dims = 1:30)
# 
# integrated_mg <- FindNeighbors(integrated_mg, reduction = "pca", dims = 1:30)
# integrated_mg <- FindClusters(integrated_mg, resolution = 0.5)
# 
# # ==============================================================================
# # 6. Visualization & Save
# # ==============================================================================
# message("--- Plotting ---")
# 
# p1 <- DimPlot(integrated_mg, reduction = "umap", group.by = "Species", raster=FALSE) + ggtitle("By Species")
# p2 <- DimPlot(integrated_mg, reduction = "umap", group.by = "integration_batch", raster=FALSE, label = TRUE) + ggtitle("By Batch") + NoLegend()
# 
# ggsave(file.path(out_dir, "integration_v4_check.pdf"), p1 + p2, width = 16, height = 7)
# 
# message("--- Saving Final Seurat V4 Object ---")
# qsave(integrated_mg, file = file.path(out_dir, "02.cross_species_integrated_seurat_v4.qs"))
# 
# message("Finished successfully!")


# 设定路径
out_dir <- Sys.getenv("MG_EXTERNAL_030")

# 1. 读取已经算好的 Anchors (救命稻草)
message("Loading calculated anchors...")
anchors <- qread(file.path(out_dir, "01.integration_anchors_v4.qs"))

# 2. 【关键步骤】清理 Anchor 对象里的 SCT 残留
message("Cleaning residual SCT assays from anchors...")

# 遍历 anchor 对象里的每一个 dataset，把 SCT 和 integrated 删干净
anchors@object.list <- lapply(anchors@object.list, function(x) {
  # 强制设为 RNA
  DefaultAssay(x) <- "RNA"
  
  # 如果有 SCT，删掉
  if ("SCT" %in% Assays(x)) {
    x[["SCT"]] <- NULL
    message(paste("Removed SCT from:", unique(x$integration_batch)))
  }
  
  # 如果有旧的 integrated，删掉
  if ("integrated" %in% Assays(x)) {
    x[["integrated"]] <- NULL
  }
  
  return(x)
})

# 3. 继续运行 IntegrateData
message("--- Resuming Integration (Cleaned) ---")
# 这里的 normalization.method 默认是 LogNormalize，这正是我们想要的
integrated_mg <- IntegrateData(anchorset = anchors)

# 4. 后续流程 (PCA, UMAP, Save)
message("--- Post-processing ---")
DefaultAssay(integrated_mg) <- "integrated"

integrated_mg <- ScaleData(integrated_mg, verbose = FALSE)
integrated_mg <- RunPCA(integrated_mg, npcs = 30, verbose = FALSE)
integrated_mg <- RunUMAP(integrated_mg, reduction = "pca", dims = 1:30)

# 聚类
integrated_mg <- FindNeighbors(integrated_mg, reduction = "pca", dims = 1:30)
integrated_mg <- FindClusters(integrated_mg, resolution = 0.5)

# 5. 画图检查
library(ggplot2)
p1 <- DimPlot(integrated_mg, reduction = "umap", group.by = "Species", raster=FALSE) + ggtitle("Species")
p2 <- DimPlot(integrated_mg, reduction = "umap", group.by = "integration_batch", raster=FALSE, label=TRUE) + ggtitle("Batch") + NoLegend()

ggsave(file.path(out_dir, "integration_v4_fixed.pdf"), p1 + p2, width = 16, height = 7)

# 6. 保存最终结果
message("--- Saving Final Clean Object ---")
qsave(integrated_mg, file = file.path(out_dir, "02.cross_species_integrated_seurat_v4.qs"))

message("Fixed and Finished!")