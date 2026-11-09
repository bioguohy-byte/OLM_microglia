setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies"))


library(qs)
library(Seurat)

manual_colors <- c(
  "MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A", 
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628"
)


combine <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/02.cross_species_integrated_noref_Human_Mouse_Monkey.qs"))
DimPlot(combine,group.by = 'CellType1', cols = manual_colors, split.by = 'Species')
DimPlot(combine,group.by = 'integrated_snn_res.0.1')


integrated_mg$new_celltype <- as.character(integrated_mg$CellType1)
target_cells <- integrated_mg$integrated_snn_res.0.1 == '1' & integrated_mg$Species == 'Mouse'
integrated_mg$new_celltype[target_cells] <- "Mouse_MG0_like"
Idents(integrated_mg) <- "new_celltype"
table(integrated_mg$new_celltype)

target_cells <- integrated_mg$integrated_snn_res.0.1 == '1' & integrated_mg$Species == 'Monkey'
integrated_mg$new_celltype[target_cells] <- "Monkey_MG0_like"
Idents(integrated_mg) <- "new_celltype"
table(integrated_mg$new_celltype)


#-------------
keep_types <- c("MG0", "Mouse_MG0_like", "Monkey_MG0_like")
integrated_mg$new_celltype <- as.character(integrated_mg$new_celltype)
integrated_mg$new_celltype <- ifelse(
  integrated_mg$new_celltype %in% keep_types,
  integrated_mg$new_celltype, 
  "Other_MG"
)
Idents(integrated_mg) <- "new_celltype"
table(integrated_mg$new_celltype, useNA = "ifany")

DimPlot(integrated_mg,group.by = 'new_celltype',split.by = 'Species')


#-----------------
# 1. 创建命名颜色向量 (Named Vector)
custom_colors <- c(
  "MG0"             = "#4DAF4A",
  "Mouse_MG0_like"  = "#33A02C",  
  "Monkey_MG0_like" = "#B2DF8A",  
  "Other_MG"        = "#D9D9D9"   
)

library(ggplot2)

(p1 <- DimPlot(integrated_mg, 
        reduction = "umap", 
        group.by = "new_celltype", 
        cols = custom_colors,
        split.by = 'Species',
        pt.size = 0.3,       # 由于 Other_MG 细胞很多，适当调整点的大小
        raster = FALSE) + 
  ggtitle("MG0 Cross-Species Subsets") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoAxes()
)
ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/figure/species_dimplot.pdf"),p1,width = 9,height = 3)
ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/figure/species_dimplot.png"),p1,width = 9,height = 3)

# 1. 为三个物种分配 3 种区分度高的经典颜色
species_colors <- c(
  "Human"  = "#E69F00",  
  "Mouse"  = "#56B4E9",  
  "Monkey" = "#009E73"   
)

# 2. 画在同一个 UMAP 当中
(p2 <- DimPlot(integrated_mg, 
              reduction = "umap", 
              group.by = "Species",     # 确保这里的拼写是 Species
              cols = species_colors,    # 应用上面定义的 3 个颜色
              pt.size = 0.2,            # 调整点的大小
              raster = FALSE) + 
  ggtitle("Species") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoAxes()                              # 取消横纵坐标及刻度
)

ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/figure/species_celltype.pdf"),p2,width = 4,height = 3)
ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/figure/species_celltype.png"),p2,width = 4,height = 3)


# process data ---------------
# mouse_wm <- qread(Sys.getenv("MG_EXTERNAL_019"))
# mouse_wm_mg <- subset(mouse_wm, subset = integrated_snn_res.0.2 %in% c('1','5','17'))
# qsave(mouse_wm_mg,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/mouse_wm_mg.qs"))
mouse_wm_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/mouse_wm_mg.qs"))


# monkey_hip <- qread(Sys.getenv("MG_EXTERNAL_020"))
# monkey_hip_mg <- subset(monkey_hip, subset = integrated_snn_res.0.1 %in% c('0','9','11'))
# monkey_hip_mg <- subset(monkey_hip_mg,subset = orig.ident %in% c('MFA_OF1','MFA_OF4','MFA_OM2','MFA_YF4','MFA_YM3') )
# qsave(monkey_hip_mg,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/monkey_hip_mg.qs"))
monkey_hip_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/monkey_hip_mg.qs"))


# ==============================================================================
# 1. 将注释映射回原始 Monkey 对象 (monkey_hip_mg)
# ==============================================================================
message("Mapping annotations back to Monkey object...")

# 提取整合对象中属于 Monkey 的细胞名和对应的注释
monkey_integrated_barcodes <- colnames(integrated_mg)[integrated_mg$Species == "Monkey"]
monkey_celltypes <- integrated_mg$new_celltype[integrated_mg$Species == "Monkey"]

# 去除 "Monkey_" 前缀，恢复原始 barcode
orig_monkey_barcodes <- gsub("^Monkey_", "", monkey_integrated_barcodes)

# 构建 Dataframe (行名必须是原始 barcode)
monkey_anno_df <- data.frame(
  new_celltype = monkey_celltypes,
  row.names = orig_monkey_barcodes
)

# 安全映射：AddMetaData 会自动根据行名(barcode)严格对齐
monkey_hip_mg <- AddMetaData(monkey_hip_mg, metadata = monkey_anno_df)
qsave(monkey_hip_mg,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/monkey_hip_mg_v2.qs"))
# 检查结果
print("Monkey Annotation Summary:")
print(table(monkey_hip_mg$new_celltype, useNA = "ifany"))


# ==============================================================================
# 2. 将注释映射回原始 Mouse 对象 (mouse_wm_mg)
# ==============================================================================
message("Mapping annotations back to Mouse object...")

# 提取整合对象中属于 Mouse 的细胞名和对应的注释
mouse_integrated_barcodes <- colnames(integrated_mg)[integrated_mg$Species == "Mouse"]
mouse_celltypes <- integrated_mg$new_celltype[integrated_mg$Species == "Mouse"]

# 去除 "Mouse_" 前缀，恢复原始 barcode
orig_mouse_barcodes <- gsub("^Mouse_", "", mouse_integrated_barcodes)

# 构建 Dataframe
mouse_anno_df <- data.frame(
  new_celltype = mouse_celltypes,
  row.names = orig_mouse_barcodes
)

# 安全映射
mouse_wm_mg <- AddMetaData(mouse_wm_mg, metadata = mouse_anno_df)
qsave(mouse_wm_mg,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/mouse_wm_mg_v2.qs"))
# 检查结果
print("Mouse Annotation Summary:")
print(table(mouse_wm_mg$new_celltype, useNA = "ifany"))

# 绘制小提琴图----------------------
#!/usr/bin/env Rscript

# MODULE : FeatureStatPlot Generation ----------------------------------------
# Description: 基于 Seurat 对象，动态生成带有显著性标记的基因表达小提琴图，
#              支持自定义分组顺序、指定对照组（Reference）并自动配对比较。


# 载入依赖包 (假设您的 FeatureStatPlot 来自 SCP 或类似可视化包)
suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(SCP)
})

## 1 Parameters --------------------------------------------------------------

# --- I/O Paths (输入与输出路径) ---
IN_QS_PATH   <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/mouse_wm_mg_v2.qs")
OUT_PDF_PATH <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/figure/Mouse_Violin_Stat_Markers.pdf")

# --- Grouping & Comparisons (细胞分组与统计比较策略) ---
GROUP_COL    <- "new_celltype"
TARGET_ORDER <- c("Mouse_MG0_like", "Other_MG")

# 设定对照组 (Reference Group)，其余组将自动与该组进行显著性检验计算
REF_GROUP    <- "Mouse_MG0_like" 

# --- Plotting Variables (绘图特征与视觉参数) ---
TARGET_GENES <- c("Mog", "Mbp", "Rnf220", "Trim2")

custom_colors <- c(
  "MG0"             = "#4DAF4A",
  "Mouse_MG0_like"  = "#33A02C",  
  "Monkey_MG0_like" = "#B2DF8A",  
  "Other_MG"        = "#D9D9D9"   
)
PALETTE_COLORS <- custom_colors 

# 图片保存尺寸
PLOT_WIDTH   <- 16
PLOT_HEIGHT  <- 4



## 2 Process Data ------------------------------------------------------------
message("--- Step 2: Loading and Processing Data ---")

# 1. 读取 Seurat 对象
seu <- qread(IN_QS_PATH)
DefaultAssay(seu) <- 'RNA'

# 2. 软编码：安全地重置因子的 level 顺序
seu[[GROUP_COL]] <- factor(seu[[GROUP_COL]][, 1], levels = TARGET_ORDER)
message("Current cluster levels: ", paste(levels(seu[[GROUP_COL]][, 1]), collapse = ", "))

# 3. 软编码：动态生成比较组合 (消灭了原来的 celltype[-3] 和 celltype[3] 硬索引)
all_clusters   <- levels(seu[[GROUP_COL]][, 1])
compare_groups <- setdiff(all_clusters, REF_GROUP) # 自动挑出所有非对照组的 cluster

# 将每一个非对照组与对照组(REF_GROUP)绑定，生成 list 格式的 comparison 参数
my_comparisons <- lapply(compare_groups, function(x) c(x, REF_GROUP))
message("Comparisons generated against Reference: ", REF_GROUP)


## 3 Generate Plots ----------------------------------------------------------
message("--- Step 3: Generating Violin Plots ---")

p_violin <- FeatureStatPlot(
  srt         = seu, 
  group.by    = GROUP_COL,
  stat.by     = TARGET_GENES, 
  add_box     = FALSE, 
  comparisons = my_comparisons,
  bg.by       = GROUP_COL,
  palcolor    = PALETTE_COLORS,
  bg_palcolor = PALETTE_COLORS,
  ncol = 4
)

# 在 Rstudio / 终端打印预览 (可选)
# print(p_violin)


## 4 Save Results ------------------------------------------------------------
message("--- Step 4: Saving Plot Results ---")

# 检查输出目录是否存在，不存在则创建
out_dir <- dirname(OUT_PDF_PATH)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 保存为 PDF
ggsave(
  filename = OUT_PDF_PATH, 
  plot     = p_violin, 
  width    = PLOT_WIDTH, 
  height   = PLOT_HEIGHT
)

message("=====================================================================")
message("SUCCESS: Plot saved to -> ", OUT_PDF_PATH)
message("=====================================================================")


# 绘制小提琴图----------------------
#!/usr/bin/env Rscript

# MODULE : FeatureStatPlot Generation ----------------------------------------
# Description: 基于 Seurat 对象，动态生成带有显著性标记的基因表达小提琴图，
#              支持自定义分组顺序、指定对照组（Reference）并自动配对比较。

rm(list=ls())
gc()

# 载入依赖包 (假设您的 FeatureStatPlot 来自 SCP 或类似可视化包)
suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(SCP)
})

## 1 Parameters --------------------------------------------------------------

# --- I/O Paths (输入与输出路径) ---
IN_QS_PATH   <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/monkey_hip_mg_v2.qs")
OUT_PDF_PATH <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/figure/Monkey_Violin_Stat_Markers.pdf")

# --- Grouping & Comparisons (细胞分组与统计比较策略) ---
GROUP_COL    <- "new_celltype"
TARGET_ORDER <- c("Monkey_MG0_like", "Other_MG")

# 设定对照组 (Reference Group)，其余组将自动与该组进行显著性检验计算
REF_GROUP    <- "Monkey_MG0_like" 

# --- Plotting Variables (绘图特征与视觉参数) ---
TARGET_GENES <- c("MOG", "MBP", "RNF220", "TRIM2")

custom_colors <- c(
  "MG0"             = "#4DAF4A",
  "Mouse_MG0_like"  = "#33A02C",  
  "Monkey_MG0_like" = "#B2DF8A",  
  "Other_MG"        = "#D9D9D9"   
)
PALETTE_COLORS <- custom_colors 

# 图片保存尺寸
PLOT_WIDTH   <- 16
PLOT_HEIGHT  <- 4



## 2 Process Data ------------------------------------------------------------
message("--- Step 2: Loading and Processing Data ---")

# 1. 读取 Seurat 对象
seu <- qread(IN_QS_PATH)
DefaultAssay(seu) <- 'RNA'

# 2. 软编码：安全地重置因子的 level 顺序
seu[[GROUP_COL]] <- factor(seu[[GROUP_COL]][, 1], levels = TARGET_ORDER)
message("Current cluster levels: ", paste(levels(seu[[GROUP_COL]][, 1]), collapse = ", "))

# 3. 软编码：动态生成比较组合 (消灭了原来的 celltype[-3] 和 celltype[3] 硬索引)
all_clusters   <- levels(seu[[GROUP_COL]][, 1])
compare_groups <- setdiff(all_clusters, REF_GROUP) # 自动挑出所有非对照组的 cluster

# 将每一个非对照组与对照组(REF_GROUP)绑定，生成 list 格式的 comparison 参数
my_comparisons <- lapply(compare_groups, function(x) c(x, REF_GROUP))
message("Comparisons generated against Reference: ", REF_GROUP)


## 3 Generate Plots ----------------------------------------------------------
message("--- Step 3: Generating Violin Plots ---")

p_violin <- FeatureStatPlot(
  srt         = seu, 
  group.by    = GROUP_COL,
  stat.by     = TARGET_GENES, 
  add_box     = FALSE, 
  comparisons = my_comparisons,
  bg.by       = GROUP_COL,
  palcolor    = PALETTE_COLORS,
  bg_palcolor = PALETTE_COLORS,
  ncol = 4
)

# 在 Rstudio / 终端打印预览 (可选)
# print(p_violin)


## 4 Save Results ------------------------------------------------------------
message("--- Step 4: Saving Plot Results ---")

# 检查输出目录是否存在，不存在则创建
out_dir <- dirname(OUT_PDF_PATH)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 保存为 PDF
ggsave(
  filename = OUT_PDF_PATH, 
  plot     = p_violin, 
  width    = PLOT_WIDTH, 
  height   = PLOT_HEIGHT
)

message("=====================================================================")
message("SUCCESS: Plot saved to -> ", OUT_PDF_PATH)
message("=====================================================================")

