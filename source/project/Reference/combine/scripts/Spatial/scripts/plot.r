# MODULE : Spatial Mapping, Subsetting, and Visualization -------------------

## 1 Parameters --------------------------------------------------------------
# 加载依赖包
library(Seurat)
library(dplyr)
library(ggplot2)

setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/"))

# 主数据路径与列名配置 (更新为 v2 版本及 cell.type_manual2)
file_master_rds <- Sys.getenv("MG_EXTERNAL_022")
col_master_celltype <- 'cell.type_manual2'

# 图片输出目录
dir_out_fig <- 'results/figure/'

# 细胞类型颜色配置
target_colors <- c(
  "ST_MG0" = "red",
  "OLG"    = "blue",
  "OPC"    = "yellow"
)
color_other <- "gray80"

# 绘图基础参数
plot_pt_size <- 0.5
plot_width   <- 6
plot_height  <- 6

# 样本级配置字典 (输出文件统一更新为 _v2 命名)
sample_configs <- list(
  "s2r1" = list(
    in_rds  = Sys.getenv("MG_EXTERNAL_023"),
    out_rds = "data/processed_data/merfish_s2r1_v2.rds",
    coord   = c(x_min = -3000, x_max = -1500, y_min = -8000, y_max = -5300)
  ),
  "s2r2" = list(
    in_rds  = Sys.getenv("MG_EXTERNAL_024"),
    out_rds = "data/processed_data/merfish_s2r2_v2.rds",
    coord   = c(x_min = 800, x_max = 4000, y_min = -10000, y_max = -7500)
  ),
  "s2r3" = list(
    in_rds  = Sys.getenv("MG_EXTERNAL_025"),
    out_rds = "data/processed_data/merfish_s2r3_v2.rds",
    coord   = c(x_min = 100, x_max = 2500, y_min = -8000, y_max = -4500)
  )
)


## 2 Process Data ------------------------------------------------------------
# 读取全局数据
master_merfish <- readRDS(file_master_rds)

# 严格采用新逻辑：提取为主数据框并重命名列
master_celltype_df <- as.data.frame(master_merfish[[col_master_celltype]])
colnames(master_celltype_df) <- 'CellType'

# 释放主数据内存
rm(master_merfish)
gc()

# 初始化存储列表
processed_objects <- list()
subset_objects <- list()

# 批量处理数据 (映射 CellType + 截取坐标子集)
for (smpl in names(sample_configs)) {
  cfg <- sample_configs[[smpl]]
  
  # 1. 读取单样本数据
  obj <- readRDS(cfg$in_rds)
  
  # 2. 映射 CellType (严格采用新逻辑：保留 dataframe 格式赋值)
  smpl_celltype <- master_celltype_df[rownames(master_celltype_df) %in% colnames(obj), , drop = FALSE]
  obj$CellType <- smpl_celltype
  processed_objects[[smpl]] <- obj # 存入列表，留待第4步保存
  
  # 3. 筛选坐标子集
  df_meta <- obj@meta.data %>% select(rotated_x, rotated_y, CellType)
  cells_keep <- df_meta %>% 
    filter(rotated_x >= cfg$coord["x_min"] & rotated_x <= cfg$coord["x_max"] & 
             rotated_y >= cfg$coord["y_min"] & rotated_y <= cfg$coord["y_max"]) %>% 
    rownames()
  
  subset_objects[[smpl]] <- subset(obj, cells = cells_keep)
}


## 3 Generate Plots ----------------------------------------------------------
# 检查并创建图片输出目录
if(!dir.exists(dir_out_fig)) dir.create(dir_out_fig, recursive = TRUE)

plot_list <- list()

# 批量生成图片
for (smpl in names(subset_objects)) {
  obj_sub <- subset_objects[[smpl]]
  
  # 动态生成颜色向量
  current_celltypes <- unique(obj_sub$CellType)
  colors_mapped <- sapply(current_celltypes, function(x) {
    if (x %in% names(target_colors)) return(target_colors[[x]])
    else return(color_other)
  })
  
  # 绘制 DimPlot
  p <- DimPlot(obj_sub, reduction = 'spatial', group.by = 'CellType', pt.size = plot_pt_size) + 
    scale_color_manual(values = colors_mapped) +
    ggtitle(paste("Spatial Plot -", smpl))
  
  plot_list[[smpl]] <- p
}


## 4 Save Results ------------------------------------------------------------
# 批量保存数据和图片
for (smpl in names(sample_configs)) {
  cfg <- sample_configs[[smpl]]
  
  # 保存映射好 CellType 的完整 RDS 对象 (此时已经是 _v2.rds)
  saveRDS(processed_objects[[smpl]], file = cfg$out_rds)
  
  # 保存图片 (PDF & PNG)
  file_base <- paste0(dir_out_fig, smpl, "_plot_v2")
  ggsave(paste0(file_base, ".pdf"), plot = plot_list[[smpl]], width = plot_width, height = plot_height)
  ggsave(paste0(file_base, ".png"), plot = plot_list[[smpl]], width = plot_width, height = plot_height)
}