rm(list=ls())
# ==============================================================================
# Seurat V4 RPCA 整合流程 (Cross-species: Human/Monkey/Mouse)
# ==============================================================================

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(qs) 
library(biomaRt) 
library(future) 

# --- 并行与内存设置 ---
# 解决 "UNRELIABLE VALUE" 报错的核心设置
options(future.rng.onMisuse = "ignore") 
options(future.globals.maxSize = 80 * 1024^3) # 80GB

plan("multicore", workers = 20) 

message("Current Seurat Version: ", packageVersion("Seurat"))

# 设定输出路径 (改为 rpca 目录)
out_dir <- Sys.getenv("MG_EXTERNAL_034")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)

# ==============================================================================
# 1. BiomaRt 连接 (Dec 2021)
# ==============================================================================
message("--- Connecting to Ensembl Archive (Dec 2021) ---")
host_url <- "https://dec2021.archive.ensembl.org/"

human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host_url)
monkey_mart <- useMart("ensembl", dataset = "mmulatta_gene_ensembl", host = host_url)
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host_url)

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
# 3. 加载数据与同源转换
# ==============================================================================
message("--- Loading and Converting Data ---")

# --- Mouse ---
mouse_wm <- qread(Sys.getenv("MG_EXTERNAL_019"))
mouse_wm_mg <- subset(mouse_wm, subset = integrated_snn_res.0.2 %in% c('1','5','17'))
rm(mouse_wm); gc()

mouse_wm_mg <- Convert_With_Mart(mouse_wm_mg, mouse_mart, human_mart, "Mouse")
mouse_wm_mg$Species <- "Mouse"
mouse_wm_mg$integration_batch <- "Mouse"
DefaultAssay(mouse_wm_mg) <- "RNA"

# --- Monkey ---
monkey_hip <- qread(Sys.getenv("MG_EXTERNAL_020"))
monkey_hip_mg <- subset(monkey_hip, subset = integrated_snn_res.0.1 %in% c('0','9','11'))
rm(monkey_hip); gc()

monkey_hip_mg <- Convert_With_Mart(monkey_hip_mg, monkey_mart, human_mart, "Monkey")
monkey_hip_mg$Species <- "Monkey"
monkey_hip_mg$integration_batch <- "Monkey"
DefaultAssay(monkey_hip_mg) <- "RNA"

# --- Human ---
human_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))
human_mg$Species <- "Human"
human_mg$integration_batch <- paste0("Human_", make.names(human_mg$Study)) 
DefaultAssay(human_mg) <- "RNA"

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
if(length(shared_genes) < 8000) warning("Shared genes < 8000. Check ortholog mapping!")

# 强制使用相同 gene space
mouse_wm_mg   <- subset(mouse_wm_mg, features = shared_genes)
monkey_hip_mg <- subset(monkey_hip_mg, features = shared_genes)
human_mg      <- subset(human_mg, features = shared_genes)

# ==============================================================================
# 4. RPCA 预处理 (Split -> SelectFeatures -> Scale/PCA each)
# ==============================================================================
message("--- Preparing Object List for RPCA ---")

# 1. 构建列表
obj.list <- list(Mouse = mouse_wm_mg, Monkey = monkey_hip_mg, Human = human_mg)
rm(mouse_wm_mg, monkey_hip_mg, human_mg); gc()

# 2. 标准化 & 找高变基因
obj.list <- lapply(X = obj.list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})

# 3. 选择整合用的 Feature
message("--- Selecting Integration Features ---")
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)

# 4. 关键：RPCA 要求对 List 中每个对象单独做 PCA
message("--- Running PCA on each object (Required for RPCA) ---")
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  return(x)
})

# ==============================================================================
# 5. Find Anchors (RPCA) & Integrate
# ==============================================================================
message("--- Finding Integration Anchors (Method = RPCA) ---")

# k.anchor: 默认是5。如果整合不够强，可以适当增加 (例如 20)
# reduction: 必须指定为 "rpca"
anchors <- FindIntegrationAnchors(object.list = obj.list, 
                                  anchor.features = features,
                                  reduction = "rpca", 
                                  dims = 1:30)

message("--- Integrating Data ---")
integrated_mg <- IntegrateData(anchorset = anchors, dims = 1:30)

# 清理内存
rm(anchors, obj.list); gc()

# ==============================================================================
# 6. Downstream Analysis (Integrated Assay)
# ==============================================================================
message("--- Post-integration Processing ---")

DefaultAssay(integrated_mg) <- "integrated"

# 标准流程：Scale -> PCA -> UMAP -> Cluster
integrated_mg <- ScaleData(integrated_mg, verbose = FALSE)
integrated_mg <- RunPCA(integrated_mg, npcs = 30, verbose = FALSE)
integrated_mg <- RunUMAP(integrated_mg, reduction = "pca", dims = 1:30)
integrated_mg <- FindNeighbors(integrated_mg, reduction = "pca", dims = 1:30)
integrated_mg <- FindClusters(integrated_mg, resolution = 0.5)

# ==============================================================================
# 7. Save FIRST (安全第一)
# ==============================================================================
message("--- Saving Final Object (Before Plotting) ---")
# 先保存，确保计算成果不丢失
qsave(integrated_mg, file = file.path(out_dir, "02.cross_species_integrated_rpca.qs"))

# ==============================================================================
# 8. Plot
# ==============================================================================
message("--- Plotting ---")

# 修正了之前的语法错误：raster 和 label 必须在 DimPlot 括号内
p1 <- DimPlot(integrated_mg, reduction = "umap", group.by = "Species", raster = FALSE) + ggtitle("By Species")

p2 <- DimPlot(integrated_mg, reduction = "umap", group.by = "integration_batch", 
              raster = FALSE, label = FALSE) + NoLegend() + ggtitle("By Batch")

p3 <- DimPlot(integrated_mg, reduction = "umap", 
              raster = FALSE, label = TRUE) + NoLegend() + ggtitle("Clusters")

ggsave(file.path(out_dir, "integration_rpca_check.pdf"), p1 + p2 + p3, width = 18, height = 6)

message("Finished! RPCA integration complete.")