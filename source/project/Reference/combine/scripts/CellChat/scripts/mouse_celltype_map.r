# 全局设置----------------
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/CellChat"))
rm(list = ls())

# library------------------
library(CellChat)
library(Seurat)
library(qs)
library(ggplot2)
library(patchwork)
library(igraph)
library(ComplexHeatmap)

# 处理数据-------------------
seurat_obj <- qread(Sys.getenv("MG_EXTERNAL_017"))
mg_meta <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/mouse_wm_mg_v2.qs"))

## 亚群细胞映射回大群FUN---------------
transfer_sub_labels <- function(seu_main, 
                                seu_sub, 
                                main_col = "CellType1",    # 主对象原始列名
                                sub_col = "CellType2",    # 子对象标签列名
                                new_col = "cellchat_type" # 新生成的列名
) {
  
  # 1. 检查列名是否存在
  if (!main_col %in% colnames(seu_main@meta.data)) {
    stop(paste("错误: 主对象 (seu_main) 中找不到列名:", main_col))
  }
  if (!sub_col %in% colnames(seu_sub@meta.data)) {
    stop(paste("错误: 子对象 (seu_sub) 中找不到列名:", sub_col))
  }
  
  # 2. 初始化：先复制主对象的标签到新列
  # 使用 as.character 确保转换为字符型，避免 factor level 对应错误的风险
  seu_main[[new_col]] <- as.character(seu_main[[main_col, drop = TRUE]])
  
  # 3. 建立映射向量
  # 提取子对象标签，并以细胞名为 names
  mapping_vec <- setNames(as.character(seu_sub[[sub_col, drop = TRUE]]), colnames(seu_sub))
  
  # 4. 找到共有细胞
  common_cells <- intersect(names(mapping_vec), colnames(seu_main))
  
  # 5. 执行替换
  if (length(common_cells) > 0) {
    seu_main@meta.data[common_cells, new_col] <- mapping_vec[common_cells]
    message(sprintf("成功: 已更新 %d 个共有细胞的标签。", length(common_cells)))
  } else {
    warning("警告: 两个对象之间没有检测到共有的细胞 ID (colnames)，未进行替换。")
  }
  
  # 6. 返回修改后的主对象
  return(seu_main)
}


seurat_obj <- transfer_sub_labels(
  seurat_obj,
  mg_meta,
  main_col = "combine_celltype",
  sub_col = "new_celltype",
  new_col = "cellchat_type"
)

# 处理数据2---------------------
Idents(seurat_obj) <- "cellchat_type"
desired_order <- c("Mouse_MG0_like","Other_MG","Neuron","OLIG","OPC","Astrocyte","Vascular","T cell","Ependymal")
seurat_obj$cellchat_type <- factor(seurat_obj$cellchat_type, levels = desired_order)
unique(seurat_obj$cellchat_type)

qsave(seurat_obj,file = 'data/processed_data/mouse_wm_cellchat.qs')

