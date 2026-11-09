rm(list=ls())
# ==============================================================================
# Seurat V4 标准无偏差整合 (No Reference)
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

# # 并行设置
# plan("multicore", workers = 10) 

message("Current Seurat Version: ", packageVersion("Seurat"))

# 设定输出路径 (注意区分文件夹，以免覆盖 Reference 的结果)
out_dir <- Sys.getenv("MG_EXTERNAL_032")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)

# ==============================================================================
# 1. BiomaRt 连接 (Dec 2021)
# ==============================================================================
message("--- Connecting to Ensembl Archive (Dec 2021) ---")
host_url <- "https://dec2021.archive.ensembl.org/"

human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host_url)
monkey_mart <- useMart("ensembl", dataset = "mmulatta_gene_ensembl", host = host_url)
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host_url)

message("BiomaRt connections established.")

# ==============================================================================
# 2. 定义转换函数
# ==============================================================================
Convert_With_Mart <- function(obj, source_mart, target_mart, species_name) {
  DefaultAssay(obj) <- "RNA"
  counts_matrix <- GetAssayData(obj, assay = "RNA", slot = "counts")
  original_genes <- rownames(counts_matrix)
  
  message(paste0("Querying ", length(original_genes), " genes for ", species_name, "..."))
  
  gene_map <- getLDS(attributes = c("external_gene_name"), 
                     filters = "external_gene_name", 
                     values = original_genes, 
                     mart = source_mart, 
                     attributesL = c("hgnc_symbol"), 
                     martL = target_mart, 
                     uniqueRows = TRUE)
  
  colnames(gene_map) <- c("Original", "Human")
  gene_map <- gene_map[gene_map$Human != "" & !is.na(gene_map$Human), ]
  
  valid_genes <- intersect(original_genes, gene_map$Original)
  counts_subset <- counts_matrix[valid_genes, ]
  
  gene_map <- gene_map[!duplicated(gene_map$Original), ]
  rownames(gene_map) <- gene_map$Original
  
  new_names <- gene_map[rownames(counts_subset), "Human"]
  
  if(any(duplicated(new_names))) {
    rownames(counts_subset) <- new_names
    counts_final <- as(rowsum(as.matrix(counts_subset), group = new_names), "dgCMatrix")
  } else {
    rownames(counts_subset) <- new_names
    counts_final <- counts_subset
  }
  
  new_obj <- CreateSeuratObject(counts = counts_final, meta.data = obj@meta.data)
  return(new_obj)
}

# ==============================================================================
# 3. 加载数据与转换
# ==============================================================================
message("--- Loading and Converting Data ---")

# --- Mouse ---
mouse_wm <- qread(Sys.getenv("MG_EXTERNAL_019"))
mouse_wm_mg <- subset(mouse_wm, subset = integrated_snn_res.0.2 %in% c('1','5','17'))
rm(mouse_wm); gc()

mouse_wm_mg <- Convert_With_Mart(mouse_wm_mg, mouse_mart, human_mart, "Mouse")
mouse_wm_mg$Species <- "Mouse"
mouse_wm_mg$integration_batch <- "Mouse"

# --- Monkey ---
monkey_hip <- qread(Sys.getenv("MG_EXTERNAL_020"))
monkey_hip_mg <- subset(monkey_hip, subset = integrated_snn_res.0.1 %in% c('0','9','11'))
rm(monkey_hip); gc()

monkey_hip_mg <- Convert_With_Mart(monkey_hip_mg, monkey_mart, human_mart, "Monkey")
monkey_hip_mg$Species <- "Monkey"
monkey_hip_mg$integration_batch <- "Monkey"

# --- Human ---
human_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))
human_mg$Species <- "Human"
human_mg$integration_batch <- paste0("Human_", make.names(human_mg$Study)) 
DefaultAssay(human_mg) <- "RNA"

###################
message("--- Checking gene numbers before intersection ---")
cat("Mouse genes:", nrow(mouse_wm_mg), "\n")
cat("Monkey genes:", nrow(monkey_hip_mg), "\n")
cat("Human genes:", nrow(human_mg), "\n")

# =============================
# 关键步骤：统一 gene space
# =============================
message("--- Taking shared gene intersection ---")

shared_genes <- Reduce(intersect, list(
  rownames(mouse_wm_mg),
  rownames(monkey_hip_mg),
  rownames(human_mg)
))

cat("Shared genes:", length(shared_genes), "\n")

# 如果 shared_genes < 10000，说明 mapping 丢失严重
if(length(shared_genes) < 8000){
  warning("Shared genes < 8000. Check ortholog mapping!")
}

# 强制三者使用相同 gene space
mouse_wm_mg   <- subset(mouse_wm_mg, features = shared_genes)
monkey_hip_mg <- subset(monkey_hip_mg, features = shared_genes)
human_mg      <- subset(human_mg, features = shared_genes)


###################



# ==============================================================================
# 4. Merge & Split
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

message("--- Finding Integration Anchors (Method = CCA, All-to-All) ---")

# 【核心修改】这里去掉了 reference 参数
# Seurat 将会在所有数据集之间两两寻找锚点 (All pairwise anchors)
anchors <- FindIntegrationAnchors(object.list = obj.list, 
                                  anchor.features = features) 
# reference = NULL (默认)

qsave(anchors, file = file.path(out_dir, "01.integration_anchors_noref.qs"))

# ==============================================================================
# 6. Clean Anchors (补丁)
# ==============================================================================
message("--- Cleaning residual SCT/Integrated assays ---")
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
rm(anchors, obj.list); gc()

message("--- Dimension Reduction ---")
DefaultAssay(integrated_mg) <- "integrated"
integrated_mg <- ScaleData(integrated_mg, verbose = FALSE)
integrated_mg <- RunPCA(integrated_mg, npcs = 30, verbose = FALSE)
integrated_mg <- RunUMAP(integrated_mg, reduction = "pca", dims = 1:30)
integrated_mg <- FindNeighbors(integrated_mg, reduction = "pca", dims = 1:30)
integrated_mg <- FindClusters(integrated_mg, resolution = 0.5)

# ==============================================================================
# 8. Save
# ==============================================================================
message("--- Plotting ---")
p1 <- DimPlot(integrated_mg, reduction = "umap", group.by = "Species", raster=FALSE) + ggtitle("By Species")
p2 <- DimPlot(integrated_mg, reduction = "umap", group.by = "integration_batch", raster=FALSE, label = TRUE) + NoLegend()

ggsave(file.path(out_dir, "integration_noref_check.pdf"), p1 + p2, width = 16, height = 7)

message("--- Saving Final Object ---")
qsave(integrated_mg, file = file.path(out_dir, "02.cross_species_integrated_noref.qs"))

message("Finished! Unbiased integration complete.")