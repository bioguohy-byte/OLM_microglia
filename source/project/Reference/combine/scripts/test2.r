rm(list=ls())
# 设置内存限制为100GB，防止大数据集整合时报错
options(future.globals.maxSize = 100 * 1024^3) 

# ==============================================================================
# 0. Library Loading & Setup
# ==============================================================================
library(dplyr)
library(Seurat) # 确保是 V4 版本
library(patchwork)
library(ggplot2)
library(qs) 
library(biomaRt) # 用于在线同源转换
library(future) 

# 并行设置
plan("multicore", workers = 10) # 建议先设10，设太高容易爆内存，根据实际情况调
message("Current Seurat Version: ", packageVersion("Seurat"))

# 设定输出路径
out_dir <- Sys.getenv("MG_EXTERNAL_031")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)

# ==============================================================================
# 1. BiomaRt Connection (Dec 2021 Archive)
# ==============================================================================
message("--- Connecting to Ensembl Archive (Dec 2021) ---")
host_url <- "https://dec2021.archive.ensembl.org/"

# 1.1 人类 (目标)
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host_url)

# 1.2 猴子 (源: 恒河猴)
# 如果你是食蟹猴，请改为 'mfascicularis_gene_ensembl'
monkey_mart <- useMart("ensembl", dataset = "mmulatta_gene_ensembl", host = host_url)

# 1.3 小鼠 (源)
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host_url)

message("BiomaRt connections established.")

# ==============================================================================
# 2. Define Conversion Function
# ==============================================================================
Convert_With_Mart <- function(obj, source_mart, target_mart, species_name) {
  
  # 强制指定使用 RNA assay
  DefaultAssay(obj) <- "RNA"
  
  # 获取原始矩阵
  counts_matrix <- GetAssayData(obj, assay = "RNA", slot = "counts")
  original_genes <- rownames(counts_matrix)
  
  message(paste0("Querying ", length(original_genes), " genes for ", species_name, "..."))
  
  # 执行查询
  # source_symbol 通常是 'external_gene_name' (猴/鼠在 Ensembl 都是这个)
  gene_map <- getLDS(attributes = c("external_gene_name"), 
                     filters = "external_gene_name", 
                     values = original_genes, 
                     mart = source_mart, 
                     attributesL = c("hgnc_symbol"), 
                     martL = target_mart, 
                     uniqueRows = TRUE)
  
  colnames(gene_map) <- c("Original", "Human")
  
  # 过滤空名
  gene_map <- gene_map[gene_map$Human != "" & !is.na(gene_map$Human), ]
  
  # 取交集保留矩阵
  valid_genes <- intersect(original_genes, gene_map$Original)
  counts_subset <- counts_matrix[valid_genes, ]
  
  # 构建映射字典 (去重)
  gene_map <- gene_map[!duplicated(gene_map$Original), ]
  rownames(gene_map) <- gene_map$Original
  
  # 获取新名字
  new_names <- gene_map[rownames(counts_subset), "Human"]
  
  # 聚合重复基因 (RowSum) - 处理多对一映射
  if(any(duplicated(new_names))) {
    message(paste0("  Aggregating ", sum(duplicated(new_names)), " duplicated human genes..."))
    rownames(counts_subset) <- new_names
    # 转矩阵聚合再转回稀疏矩阵 (处理同名基因相加)
    # 注意：如果内存不够，这里可能会卡一下，但通常子集数据没问题
    counts_final <- as(rowsum(as.matrix(counts_subset), group = new_names), "dgCMatrix")
  } else {
    rownames(counts_subset) <- new_names
    counts_final <- counts_subset
  }
  
  message(paste0("  Final features: ", nrow(counts_final)))
  
  # 重建对象
  new_obj <- CreateSeuratObject(counts = counts_final, meta.data = obj@meta.data)
  return(new_obj)
}

# ==============================================================================
# 3. Load Data & Convert
# ==============================================================================
message("--- Loading and Converting Data ---")

# --- 3.1 Mouse ---
message("Processing Mouse Data...")
mouse_wm <- qread(Sys.getenv("MG_EXTERNAL_019"))
mouse_wm_mg <- subset(mouse_wm, subset = integrated_snn_res.0.2 %in% c('1','5','17'))
rm(mouse_wm); gc()

# 转换 Mouse
mouse_wm_mg <- Convert_With_Mart(mouse_wm_mg, mouse_mart, human_mart, "Mouse")
mouse_wm_mg$Species <- "Mouse"
mouse_wm_mg$integration_batch <- "Mouse"

# --- 3.2 Monkey ---
message("Processing Monkey Data...")
monkey_hip <- qread(Sys.getenv("MG_EXTERNAL_020"))
monkey_hip_mg <- subset(monkey_hip, subset = integrated_snn_res.0.1 %in% c('0','9','11'))
rm(monkey_hip); gc()

# 转换 Monkey
monkey_hip_mg <- Convert_With_Mart(monkey_hip_mg, monkey_mart, human_mart, "Monkey")
monkey_hip_mg$Species <- "Monkey"
monkey_hip_mg$integration_batch <- "Monkey"

# --- 验证 PTPRC 是否还在 (Crucial Check) ---
if("PTPRC" %in% rownames(monkey_hip_mg)) {
  message("✅ SUCCESS: PTPRC is present in the converted Monkey object!")
} else {
  warning("❌ WARNING: PTPRC is missing! Check gene mapping.")
}

# --- 3.3 Human ---
message("Processing Human Data...")
human_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))
human_mg$Species <- "Human"
# Human 不需要转换，直接用
human_mg$integration_batch <- paste0("Human_", make.names(human_mg$Study)) 
# 确保使用 RNA assay
DefaultAssay(human_mg) <- "RNA"

# ==============================================================================
# 4. Prepare Integration (Merge & Split)
# ==============================================================================
message("--- Merging datasets ---")

combined_raw <- merge(
  x = mouse_wm_mg,
  y = list(monkey_hip_mg, human_mg),
  add.cell.ids = c("Mouse", "Monkey", "Human") 
)
rm(mouse_wm_mg, monkey_hip_mg, human_mg); gc()

message("--- Splitting object by integration_batch ---")
obj.list <- SplitObject(combined_raw, split.by = "integration_batch")
rm(combined_raw); gc()

# ==============================================================================
# 5. Pre-processing & Find Anchors
# ==============================================================================
message("--- Normalizing and finding features ---")

# 标准预处理
obj.list <- lapply(X = obj.list, FUN = function(x) {
  DefaultAssay(x) <- "RNA" # 再次确保是 RNA
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

message("--- Selecting integration features ---")
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)

message("--- Finding Integration Anchors (Method = CCA) ---")

# [优化] 使用人类数据作为 Reference，防止猴子游离
# 自动寻找所有人类样本的索引
ref_idx <- grep("Human", names(obj.list))
message(paste("Using Human datasets as reference indices:", paste(ref_idx, collapse=",")))

anchors <- FindIntegrationAnchors(object.list = obj.list, 
                                  anchor.features = features, 
                                  reference = ref_idx) # 指定参考系

# *** 保存中间结果 ***
qsave(anchors, file = file.path(out_dir, "01.integration_anchors_v4.qs"))

# ==============================================================================
# 6. Clean Anchors (修复 SCT/Integrated 报错的关键步骤)
# ==============================================================================
message("--- Cleaning residual SCT/Integrated assays from anchors ---")

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

# ==============================================================================
# 7. Integrate & Post-process
# ==============================================================================
message("--- Integrating Data ---")
integrated_mg <- IntegrateData(anchorset = anchors)
rm(anchors, obj.list); gc()

message("--- Dimension Reduction ---")
DefaultAssay(integrated_mg) <- "integrated"

integrated_mg <- ScaleData(integrated_mg, verbose = FALSE)
integrated_mg <- RunPCA(integrated_mg, npcs = 30, verbose = FALSE)
integrated_mg <- RunUMAP(integrated_mg, reduction = "pca", dims = 1:30)

integrated_mg <- FindNeighbors(integrated_mg, reduction = "pca", dims = 1:30)
integrated_mg <- FindClusters(integrated_mg, resolution = 0.5)

# ==============================================================================
# 8. Visualization & Save
# ==============================================================================
message("--- Plotting ---")

p1 <- DimPlot(integrated_mg, reduction = "umap", group.by = "Species", raster=FALSE) + ggtitle("By Species")
p2 <- DimPlot(integrated_mg, reduction = "umap", group.by = "integration_batch", raster=FALSE, label = TRUE) + ggtitle("By Batch") + NoLegend()

ggsave(file.path(out_dir, "integration_biomart_fixed.pdf"), p1 + p2, width = 16, height = 7)

message("--- Saving Final Object ---")
qsave(integrated_mg, file = file.path(out_dir, "02.cross_species_integrated_seurat_v4_biomart.qs"))

message("Finished successfully! PTPRC should now be correct.")