rm(list=ls())
# ==============================================================================
# Seurat + Harmony 整合流程 (Cross-species)
# ==============================================================================

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(qs) 
library(biomaRt) 
library(future) 
library(harmony) # 必需加载 Harmony

library(future) # 新增
plan("multicore", workers = 20) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限


message("Current Seurat Version: ", packageVersion("Seurat"))

# 设定输出路径 (改为 harmony 目录)
out_dir <- Sys.getenv("MG_EXTERNAL_033")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)

# ==============================================================================
# 1. BiomaRt 连接 (Dec 2021) - 保持不变
# ==============================================================================
message("--- Connecting to Ensembl Archive (Dec 2021) ---")
host_url <- "https://dec2021.archive.ensembl.org/"

human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host_url)
monkey_mart <- useMart("ensembl", dataset = "mmulatta_gene_ensembl", host = host_url)
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host_url)

# ==============================================================================
# 2. 定义转换函数 - 保持不变
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
# 3. 加载数据与转换 - 保持不变
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
# 4. 预处理 (标准化与高变基因选择) - 为 Harmony 做准备
# ==============================================================================
message("--- Creating Object List for Feature Selection ---")
# 为了更准确地去除批次效应，我们先单独对每个数据集找高变基因，然后取并集或使用 SelectIntegrationFeatures
obj.list <- list(Mouse = mouse_wm_mg, Monkey = monkey_hip_mg, Human = human_mg)
rm(mouse_wm_mg, monkey_hip_mg, human_mg); gc()

obj.list <- lapply(X = obj.list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})

message("--- Selecting Integration Features ---")
# 这步确保我们选出的基因在不同批次中都是高变的
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)

message("--- Merging Data ---")
# Harmony 需要一个合并的大对象 (unintegrated counts)
combined_obj <- merge(
  x = obj.list[[1]],
  y = obj.list[-1],
  add.cell.ids = c("Mouse", "Monkey", "Human") # 注意这里顺序要和 obj.list 一致
)
rm(obj.list); gc()

DefaultAssay(combined_obj) <- "RNA"
VariableFeatures(combined_obj) <- features

# ==============================================================================
# 5. Run PCA & Harmony
# ==============================================================================
message("--- Scaling and PCA ---")
# ScaleData 只需要对选定的 features 进行，节省内存
combined_obj <- ScaleData(combined_obj, features = features, verbose = FALSE)
combined_obj <- RunPCA(combined_obj, npcs = 30, verbose = FALSE)

message("--- Running Harmony ---")
# group.by.vars: 你想去除的批次效应列名
# plot_convergence: 是否画图看收敛过程
combined_obj <- RunHarmony(combined_obj, 
                           group.by.vars = "integration_batch", 
                           plot_convergence = TRUE)

# 保存 Harmony 迭代图（可选）
ggsave(file.path(out_dir, "harmony_convergence.pdf"))

# ==============================================================================
# 6. Downstream Analysis (基于 Harmony reduction)
# ==============================================================================
message("--- Dimension Reduction (UMAP & Clustering) ---")

# 注意：reduction 必须设为 'harmony'
combined_obj <- RunUMAP(combined_obj, reduction = "harmony", dims = 1:30)
combined_obj <- FindNeighbors(combined_obj, reduction = "harmony", dims = 1:30)
combined_obj <- FindClusters(combined_obj, resolution = 0.5)

# ==============================================================================
# 7. Save & Plot
# ==============================================================================

message("--- Saving Final Object ---")
qsave(combined_obj, file = file.path(out_dir, "02.cross_species_integrated_harmony.qs"))
message("--- Plotting ---")
# group.by 使用默认的 identity (即 clusters) 和 metadata
p1 <- DimPlot(combined_obj, reduction = "umap", group.by = "Species", raster=FALSE)
p2 <- DimPlot(combined_obj, reduction = "umap", group.by = "integration_batch", raster=FALSE, label = FALSE) + NoLegend() 
p3 <- DimPlot(combined_obj, reduction = "umap", label = TRUE) + NoLegend() 

ggsave(file.path(out_dir, "integration_harmony_check.pdf"), p1 + p2 + p3, width = 18, height = 6)



message("Finished! Harmony integration complete.")