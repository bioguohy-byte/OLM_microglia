# MODULE : Cell Type Annotation and Module Scoring -------------------

## 1 Parameters --------------------------------------------------------------
# 加载依赖包
library(Seurat)
library(qs)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SCP)

# 文件路径
file_in_merfish  <- Sys.getenv("MG_EXTERNAL_026")
file_in_markers  <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/saved_data/mouse_mg_markers.qs")
file_out_merfish <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/merfish_immune.combined_with4.rds")
file_out_plot    <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/results/figure/module_score_myelin.pdf")

# 细胞类型映射字典 (左边是旧名字，右边是新名字)
col_old_celltype <- 'cell.type_manual'
col_new_celltype <- 'plot_celltype'
celltype_dict <- c(
  'Myelin_MG' = 'ST_MG0',
  'MG_0'      = 'ST_MG1',
  'MG_1'      = 'ST_MG2',
  'MG_2'      = 'ST_MG3',
  'MG_3'      = 'ST_MG4',
  'MG_5'      = 'ST_MG5'
)

# Marker 基因提取与 Module Score 参数
marker_cluster_col <- 'cluster'
target_cluster     <- 'cluster2'
marker_gene_col    <- 'gene'

score_module_name  <- "myelin"
score_feature_name <- paste0(score_module_name, "1") # Seurat AddModuleScore 默认会自动加 '1'
score_ctrl_size    <- 10

# 绘图参数
plot_pt_size <- 0


## 2 Process Data ------------------------------------------------------------
# 2.1 读取数据
merfish <- readRDS(file_in_merfish)
mouse_mg_markers <- qread(file_in_markers)

# 2.2 替换细胞类型名称
# 初始化新列
merfish@meta.data[[col_new_celltype]] <- as.character(merfish@meta.data[[col_old_celltype]])
# 遍历字典进行替换
for (old_name in names(celltype_dict)) {
  new_name <- celltype_dict[[old_name]]
  merfish@meta.data[[col_new_celltype]][merfish@meta.data[[col_new_celltype]] %in% old_name] <- new_name
}

# 2.3 提取目标 cluster 基因并与 MERFISH 数据取交集
filtered_data <- mouse_mg_markers %>% 
  filter(.data[[marker_cluster_col]] == target_cluster)
int_gene <- intersect(filtered_data[[marker_gene_col]], rownames(merfish))

# 2.4 计算 Module Score
merfish <- AddModuleScore(object = merfish,
                          features = list(int_gene),
                          ctrl = score_ctrl_size,
                          name = score_module_name)


## 3 Generate Plots ----------------------------------------------------------
merfish <- readRDS(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/merfish_immune.combined_with4.rds"))

st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35", 
  "ST_MG1" = "#4DBBD5", 
  "ST_MG2" = "#00A087", 
  "ST_MG3" = "#3C5488", 
  "ST_MG4" = "#F39B7F", 
  "ST_MG5" = "#7E6148"
)

(p1 <- DimPlot(merfish, group.by = col_new_celltype,cols = st_mg_colors_npg))

# 生成 VlnPlot
target_order <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")
merfish$plot_celltype <- factor(merfish$plot_celltype, levels = target_order)

celltype <- levels(merfish$plot_celltype)
my_comparisons <- lapply(celltype[-1], function(x)
  c(x,celltype[1]))

p2 <- FeatureStatPlot(
  srt = merfish, group.by = "plot_celltype",
  stat.by = c("myelin1"), add_box = F,
  comparisons = my_comparisons,
  bg.by = 'plot_celltype',
  palcolor = st_mg_colors_npg,
  bg_palcolor = st_mg_colors_npg
)


# 拼接图形
combined_plot <- p1 | p2


## 4 Save Results ------------------------------------------------------------
# 保存更新后的 Seurat 对象 (包含新添加的 celltype 和 module score)
saveRDS(merfish, file = file_out_merfish)

# 保存可视化结果
ggsave(filename = file_out_plot, plot = combined_plot, width = 10, height = 5)
