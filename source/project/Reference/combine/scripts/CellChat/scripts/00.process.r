# 处理Chia文章中的DLB路易痴呆症的单细胞数据
# 因为MAM在DLB中的细胞比例是最多的，我们准备单独做DLB的细胞互作，与正常组的细胞互作，两者对比一下，看哪些通路是有差异的

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
seurat_obj <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Chia_et_al./data/processed_data/Chia_v2.qs"))
sce <- SplitObject(seurat_obj, split.by = "Condition")

mg_meta <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))

## 亚群细胞映射回大群FUN---------------
transfer_sub_labels <- function(seu_main, 
                                seu_sub, 
                                main_col = "CellType",    # 主对象原始列名
                                sub_col = "CellType1",    # 子对象标签列名
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


Chia_DLB <- sce[["DLB"]]
Chia_DLB <- transfer_sub_labels(
  Chia_DLB,
  mg_meta,
  main_col = "CellType",
  sub_col = "CellType1",
  new_col = "cellchat_type"
)

Chia_CTRL <- sce[["CTRL"]]
Chia_CTRL <- transfer_sub_labels(
  Chia_CTRL,
  mg_meta,
  main_col = "CellType",
  sub_col = "CellType1",
  new_col = "cellchat_type"
)

# 处理数据2---------------------
Idents(Chia_DLB) <- "cellchat_type"
Idents(Chia_CTRL) <- "cellchat_type"

qsave(Chia_DLB,file = 'data/processed_data/Chia_DLB.qs')
qsave(Chia_CTRL,file = 'data/processed_data/Chia_CTRL.qs')


