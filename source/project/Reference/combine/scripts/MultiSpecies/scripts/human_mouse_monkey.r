#!/usr/bin/env Rscript

# MODULE : Multi-Species (Human, Mouse, Monkey) snRNA-seq Integration --------
# Description: 无偏差 CCA 整合人、鼠、猴三个物种的单细胞转录组数据 (Seurat V4).
#              基于同源基因 (Ensembl Dec 2021) 统一转换为 Human Gene Symbol.
#              并按模块化、非硬编码要求重构，增加中间处理数据落盘。

rm(list=ls())
gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(qs)
  library(biomaRt)
  library(future)
})

message("Current Seurat Version: ", packageVersion("Seurat"))

## 1 Parameters --------------------------------------------------------------

# --- Global & Environment Settings ---
MAX_MEMORY_GB <- 100
WORKERS <- 20
options(future.globals.maxSize = MAX_MEMORY_GB * 1024^3)
plan("multicore", workers = WORKERS)

# --- I/O Paths ---
# 输出路径
OUT_DIR <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data")

# 中间结果保存路径 (新增)
MOUSE_PROCESSED_QS <- file.path(OUT_DIR, "00.mouse_processed.qs")
MONKEY_PROCESSED_QS <- file.path(OUT_DIR, "00.monkey_processed.qs")

# 输入文件路径
HUMAN_DATA_PATH  <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs")
MOUSE_DATA_PATH  <- Sys.getenv("MG_EXTERNAL_019")
MONKEY_DATA_PATH <- Sys.getenv("MG_EXTERNAL_020")

# --- Filtering Parameters ---
# Mouse 过滤条件
MOUSE_CLUSTERS <- c('1', '5', '17')
MOUSE_CLUSTER_COL <- "integrated_snn_res.0.2"

# Monkey 过滤条件
MONKEY_CLUSTERS <- c('0', '9', '11')
MONKEY_CLUSTER_COL <- "integrated_snn_res.0.1"
MONKEY_SAMPLES <- c('MFA_OF1', 'MFA_OF4', 'MFA_OM2', 'MFA_YF4', 'MFA_YM3')

# --- BiomaRt Settings ---
BIOMART_HOST <- "https://dec2021.archive.ensembl.org/"

# --- Seurat Analysis Parameters ---
N_FEATURES <- 2000
PCA_DIMS <- 1:30
UMAP_DIMS <- 1:10
CLUSTER_RES <- 0.1

# --- Plotting Colors ---
SPECIES_COLORS <- c(
  "Human"  = "#377EB8",
  "Mouse"  = "#E41A1C",
  "Monkey" = "#4DAF4A"
)

## 2 Process Data ------------------------------------------------------------

# 确保输出文件夹存在
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# -- 2.1 定义同源基因转换函数 --
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

# -- 2.2 连接 Ensembl BiomaRt --
message("--- Connecting to Ensembl Archive (Dec 2021) ---")
human_mart  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = BIOMART_HOST)
mouse_mart  <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = BIOMART_HOST)
monkey_mart <- useMart("ensembl", dataset = "mmulatta_gene_ensembl", host = BIOMART_HOST)
message("BiomaRt connections established.")

# -- 2.3 加载并处理 Human 数据 --
message("--- Processing Human Data ---")
human_mg <- qread(HUMAN_DATA_PATH)
DefaultAssay(human_mg) <- "RNA"
human_mg$Species <- "Human"
human_mg$integration_batch <- paste0("Human_", make.names(human_mg$Study)) 

# -- 2.4 加载并处理 Mouse 数据 --
message("--- Processing Mouse Data ---")
mouse_wm <- qread(MOUSE_DATA_PATH)
mouse_meta_clusters <- mouse_wm[[MOUSE_CLUSTER_COL]][,1]
mouse_wm_mg <- subset(mouse_wm, cells = colnames(mouse_wm)[mouse_meta_clusters %in% MOUSE_CLUSTERS])
rm(mouse_wm); gc()

mouse_wm_mg <- Convert_With_Mart(mouse_wm_mg, mouse_mart, human_mart, "Mouse")
mouse_wm_mg$Species <- "Mouse"
mouse_wm_mg$integration_batch <- "Mouse"

message(paste0("--- Saving Processed Mouse Data to: ", MOUSE_PROCESSED_QS, " ---"))
qsave(mouse_wm_mg, MOUSE_PROCESSED_QS)

# -- 2.5 加载并处理 Monkey 数据 --
message("--- Processing Monkey Data ---")
monkey_hip <- qread(MONKEY_DATA_PATH)
monkey_meta_clusters <- monkey_hip[[MONKEY_CLUSTER_COL]][,1]
monkey_hip_mg <- subset(monkey_hip, cells = colnames(monkey_hip)[monkey_meta_clusters %in% MONKEY_CLUSTERS])
monkey_hip_mg <- subset(monkey_hip_mg, subset = orig.ident %in% MONKEY_SAMPLES)
rm(monkey_hip); gc()

monkey_hip_mg <- Convert_With_Mart(monkey_hip_mg, monkey_mart, human_mart, "Monkey")
monkey_hip_mg$Species <- "Monkey"
monkey_hip_mg$integration_batch <- "Monkey"

message(paste0("--- Saving Processed Monkey Data to: ", MONKEY_PROCESSED_QS, " ---"))
qsave(monkey_hip_mg, MONKEY_PROCESSED_QS)

# # -- 2.6 获取共同基因并合并 --
# message("--- Checking gene numbers before intersection ---")
# cat("Human genes:", nrow(human_mg), "\n")
# cat("Mouse genes:", nrow(mouse_wm_mg), "\n")
# cat("Monkey genes:", nrow(monkey_hip_mg), "\n")
# 
# shared_genes <- intersect(intersect(rownames(human_mg), rownames(mouse_wm_mg)), rownames(monkey_hip_mg))
# cat("Shared genes across all 3 species:", length(shared_genes), "\n")
# 
# if(length(shared_genes) < 8000){
#   warning("Shared genes < 8000. Please double check ortholog mapping quality!")
# }
# 
# # 强制统一基因空间并合并
# human_mg      <- subset(human_mg, features = shared_genes)
# mouse_wm_mg   <- subset(mouse_wm_mg, features = shared_genes)
# monkey_hip_mg <- subset(monkey_hip_mg, features = shared_genes)
# 
# message("--- Merging 3 Species Datasets ---")
# combined_raw <- merge(
#   x = human_mg,
#   y = list(mouse_wm_mg, monkey_hip_mg),
#   add.cell.ids = c("Human", "Mouse", "Monkey")
# )
# rm(human_mg, mouse_wm_mg, monkey_hip_mg); gc()
# 
# # -- 2.7 拆分与预处理 (Split & Pre-process) --
# message("--- Splitting object by integration_batch ---")
# obj.list <- SplitObject(combined_raw, split.by = "integration_batch")
# rm(combined_raw); gc()
# 
# message("--- Normalizing and finding features ---")
# obj.list <- lapply(X = obj.list, FUN = function(x) {
#   DefaultAssay(x) <- "RNA" 
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = N_FEATURES)
# })
# 
# # -- 2.8 CCA 整合核心步骤 --
# message("--- Selecting integration features ---")
# features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = N_FEATURES)
# # 及时保存特征基因作为分析检查点
# qsave(features, file = file.path(OUT_DIR, "00.integration_features.qs"))
# 
# message("--- Finding Integration Anchors (Method = CCA) ---")
# anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features) 
# # 及时保存 Anchors (耗时极长，作为关键检查点)
# qsave(anchors, file = file.path(OUT_DIR, "01.integration_anchors_noref.qs"))
# 
# message("--- Cleaning residual assays to save memory ---")
# anchors@object.list <- lapply(anchors@object.list, function(x) {
#   DefaultAssay(x) <- "RNA"
#   if ("SCT" %in% Assays(x)) x[["SCT"]] <- NULL
#   if ("integrated" %in% Assays(x)) x[["integrated"]] <- NULL
#   return(x)
# })
# gc()
# 
# message("--- Integrating Data ---")
# integrated_mg <- IntegrateData(anchorset = anchors)
# rm(anchors, obj.list); gc()
# 
# # -- 2.9 降维与聚类 (DR & Clustering) --
# message("--- Dimension Reduction and Clustering ---")
# DefaultAssay(integrated_mg) <- "integrated"
# integrated_mg <- ScaleData(integrated_mg, verbose = FALSE)
# integrated_mg <- RunPCA(integrated_mg, npcs = max(PCA_DIMS), verbose = FALSE)
# integrated_mg <- RunUMAP(integrated_mg, reduction = "pca", dims = UMAP_DIMS)
# integrated_mg <- FindNeighbors(integrated_mg, reduction = "pca", dims = UMAP_DIMS)
# integrated_mg <- FindClusters(integrated_mg, resolution = CLUSTER_RES)
# 
# 
# # 3 Generate Plots ----------------------------------------------------------
# 
# message("--- Generating QA Plots ---")
# 
# p1 <- DimPlot(integrated_mg, reduction = "umap", group.by = "Species", raster = FALSE) + 
#   ggtitle("By Species") + 
#   theme(plot.title = element_text(hjust = 0.5))
# 
# p2 <- DimPlot(integrated_mg, reduction = "umap", group.by = "integration_batch", raster = FALSE, label = TRUE) + 
#   NoLegend() + 
#   ggtitle("By Batch") + 
#   theme(plot.title = element_text(hjust = 0.5))
# 
# cluster_col_name <- paste0("integrated_snn_res.", CLUSTER_RES)
# p3 <- DimPlot(integrated_mg, reduction = "umap", group.by = cluster_col_name, label = TRUE) + 
#   ggtitle(paste0("By Cluster (res ", CLUSTER_RES, ")"))
# 
# combined_plot <- (p1 + p2) / p3
# 
# 
# ## 4 Save Results ------------------------------------------------------------
# 
# message("--- Saving Final Outputs ---")
# 
# # 1. 保存最终的 Seurat 对象
# final_obj_path <- file.path(OUT_DIR, "02.cross_species_integrated_noref_Human_Mouse_Monkey.qs")
# qsave(integrated_mg, file = final_obj_path)
# message("Saved integrated object to: ", final_obj_path)
# 
# # 2. 保存质控图纸
# pdf_path <- file.path(OUT_DIR, "integration_noref_Human_Mouse_Monkey_check.pdf")
# ggsave(pdf_path, combined_plot, width = 14, height = 12)
# message("Saved check plots to: ", pdf_path)
# 
# message("=====================================================================")
# message("SUCCESS: Multi-Species (Human, Mouse, Monkey) integration complete!")
# message("=====================================================================")