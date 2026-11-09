rm(list=ls())

# 0. 环境设置与配置 (Setup) -------------------------------------------------------
# 设置内存限制为 100GB
options(future.globals.maxSize = 100 * 1024^3) 

library(dplyr)
library(Seurat) # V4
library(patchwork)
library(ggplot2)
library(qs) 
library(biomaRt) 
library(future) 

# 并行设置
plan("multicore", workers = 20)

message("Current Seurat Version: ", packageVersion("Seurat"))

# 设定输出路径
out_dir <- Sys.getenv("MG_EXTERNAL_036")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)


# 1. 工具函数定义 (Functions) -----------------------------------------------------

# 1.1 定义同源基因转换函数 ----
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


# 2. 数据加载与预处理 (Data Loading) ----------------------------------------------

# 2.1 连接 Ensembl BiomaRt ----
message("--- Connecting to Ensembl Archive (Dec 2021) ---")
host_url <- "https://dec2021.archive.ensembl.org/"

human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host_url)
monkey_mart <- useMart("ensembl", dataset = "mmulatta_gene_ensembl", host = host_url)
message("BiomaRt connections established.")

# 2.2 处理 Monkey 数据 ----
message("Processing Monkey Data...")
monkey_hip <- qread(Sys.getenv("MG_EXTERNAL_020"))

# 筛选 Monkey 特定的 cluster
monkey_hip_mg <- subset(monkey_hip, subset = integrated_snn_res.0.1 %in% c('0','9','11'))

## 选两个样本的猴子
monkey_hip_mg <- subset(monkey_hip_mg,subset = orig.ident %in% c('MFA_OF1','MFA_OF4','MFA_OM2','MFA_YF4','MFA_YM3') )

rm(monkey_hip); gc()

# 转换为 Human Gene Symbol
monkey_hip_mg <- Convert_With_Mart(monkey_hip_mg, monkey_mart, human_mart, "Monkey")
monkey_hip_mg$Species <- "Monkey"
monkey_hip_mg$integration_batch <- "Monkey"

# 2.3 处理 Human 数据 ----
message("Processing Human Data...")
human_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))

human_mg$Species <- "Human"
human_mg$integration_batch <- paste0("Human_", make.names(human_mg$Study)) 
DefaultAssay(human_mg) <- "RNA"

# 2.4 基因交集检查 (Intersection) ----
message("--- Checking gene numbers before intersection ---")
cat("Monkey genes:", nrow(monkey_hip_mg), "\n")
cat("Human genes:", nrow(human_mg), "\n")

message("--- Taking shared gene intersection ---")
shared_genes <- intersect(rownames(monkey_hip_mg), rownames(human_mg))
cat("Shared genes:", length(shared_genes), "\n")

if(length(shared_genes) < 8000){
  warning("Shared genes < 8000. Check ortholog mapping quality!")
}

# 强制统一基因空间
monkey_hip_mg <- subset(monkey_hip_mg, features = shared_genes)
human_mg      <- subset(human_mg, features = shared_genes)


# 3. 合并与拆分 (Merge & Split) ---------------------------------------------------

message("--- Merging datasets ---")
combined_raw <- merge(
  x = human_mg,
  y = monkey_hip_mg,
  add.cell.ids = c("Human", "Monkey") 
)
rm(monkey_hip_mg, human_mg); gc()

message("--- Splitting object by integration_batch ---")
obj.list <- SplitObject(combined_raw, split.by = "integration_batch")
rm(combined_raw); gc()


# 4. 寻找整合锚点 (Find Anchors) --------------------------------------------------

# 4.1 标准化与高变基因 (Norm & HVG) ----
message("--- Normalizing and finding features ---")
obj.list <- lapply(X = obj.list, FUN = function(x) {
  DefaultAssay(x) <- "RNA" 
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# 4.2 选择整合特征 (Select Features) ----
message("--- Selecting integration features ---")
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)
qsave(features, file = file.path(out_dir, "00.integration_features.qs"))

# 4.3 运行 CCA (FindIntegrationAnchors) ----
message("--- Finding Integration Anchors (Method = CCA) ---")
anchors <- FindIntegrationAnchors(object.list = obj.list, 
                                  anchor.features = features) 
# 保存 Anchors
qsave(anchors, file = file.path(out_dir, "01.integration_anchors_noref.qs"))

# 4.4 清理 Anchors 对象 (Clean) ----
message("--- Cleaning residual assays ---")
anchors@object.list <- lapply(anchors@object.list, function(x) {
  DefaultAssay(x) <- "RNA"
  if ("SCT" %in% Assays(x)) x[["SCT"]] <- NULL
  if ("integrated" %in% Assays(x)) x[["integrated"]] <- NULL
  return(x)
})


# 5. 整合与降维 (Integration & DR) ------------------------------------------------

# 5.1 运行整合 (IntegrateData) ----
message("--- Integrating Data ---")
integrated_mg <- IntegrateData(anchorset = anchors)
rm(anchors, obj.list); gc()

# 5.2 降维与聚类 (PCA/UMAP/Cluster) ----
message("--- Dimension Reduction ---")
DefaultAssay(integrated_mg) <- "integrated"
integrated_mg <- ScaleData(integrated_mg, verbose = FALSE)
integrated_mg <- RunPCA(integrated_mg, npcs = 30, verbose = FALSE)

# UMAP
integrated_mg <- RunUMAP(integrated_mg, reduction = "pca", dims = 1:10)

# Clustering
integrated_mg <- FindNeighbors(integrated_mg, reduction = "pca", dims = 1:10)
integrated_mg <- FindClusters(integrated_mg, resolution = 0.5)


# 6. 结果输出 (Output) ------------------------------------------------------------

# 6.1 保存对象 (Save Object) ----
message("--- Saving Final Object ---")
qsave(integrated_mg, file = file.path(out_dir, "02.cross_species_integrated_noref.qs"))

# 6.2 绘图检查 (Plotting) ----
message("--- Plotting ---")
# 图1：按物种
p1 <- DimPlot(integrated_mg, reduction = "umap", group.by = "Species", raster=FALSE) + 
  ggtitle("By Species") + 
  theme(plot.title = element_text(hjust = 0.5))

# 图2：按批次
p2 <- DimPlot(integrated_mg, reduction = "umap", group.by = "integration_batch", raster=FALSE, label = TRUE) + 
  NoLegend() + 
  ggtitle("By Batch") + 
  theme(plot.title = element_text(hjust = 0.5))

# 图3：按聚类
p3 <- DimPlot(integrated_mg, reduction = "umap", group.by = "integrated_snn_res.0.5", label = TRUE) + 
  ggtitle("By Cluster (res 0.5)")

combined_plot <- (p1 + p2) / p3
ggsave(file.path(out_dir, "integration_noref_check.pdf"), combined_plot, width = 14, height = 12)

message("Finished! Human-Monkey integration complete.")