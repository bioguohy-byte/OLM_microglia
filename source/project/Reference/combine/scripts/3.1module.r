# MODULE : UMAP/tSNE combined with Nuclei Statistics Barplot -------------------

library(dplyr)
library(ggplot2)
library(patchwork)
library(forcats)
library(ggh4x)

# ## 1 Parameters(Study) --------------------------------------------------------------
# # 数据与分组参数
# seurat_obj  <- tmp_obj                   # 输入的 Seurat 对象
# group_col   <- "Study"                   # 用于分组的列名 (e.g., "Study", "Sample")
# type_col    <- "CellType1"               # 细胞类型列名 (预留，以防后续扩展)
# 
# # 颜色参数 (NPG 风格)
# manual_colors <- c(
#   "Chia et al., 2025"    = "#E64B35",
#   "Fullard et al., 2021" = "#4DBBD5",
#   "Grubman et al., 2019" = "#00A087",
#   "Jakel et al., 2019"   = "#3C5488",
#   "Su et al., 2022"      = "#F39B7F",
#   "Yao et al., 2024"     = "#8491B4"
# )
# 
# # 图像文本与布局参数
# plot_title    <- "Number of nuclei (mean genes per nucleus)"
# text_nudge_y  <- 200                     # 柱状图上方文本的偏移量 (可根据最大细胞数动态微调)
# bar_width     <- 0.8                     # 柱子的宽度
# 
# # 输出参数
# output_dir    <- "./results/final_fig/"  # 输出文件夹
# file_prefix   <- "DimPlot_and_Bar_"      # 输出文件前缀
# save_width    <- 12                      # 画布宽度
# save_height   <- 4                       # 画布高度

# ## 1 Parameters(Region) --------------------------------------------------------------
# # 数据与分组参数
# seurat_obj  <- tmp_obj                   # 输入的 Seurat 对象
# group_col   <- "Region"                  # 【修改点 1】: 将分组列名改为 Region
# type_col    <- "CellType1"               # 细胞类型列名 (预留，以防后续扩展)
# 
# # 颜色参数 (高对比度 D3 风格，与之前的颜色完全不同)
# # 【修改点 2】: 替换为对应脑区名称的新颜色向量
# manual_colors <- c(
#   "Choroid Plexus"    = "#FF7F0E",  # 鲜亮橙色
#   "Entorhinal Cortex" = "#2CA02C",  # 森林绿色
#   "Hippocampus"       = "#8B2F5D",  # 酒红色
#   "Medulla"           = "#8C564B",  # 栗棕色
#   "Prefrontal Cortex" = "#E377C2",  # 柔和洋红/粉色
#   "White Matter"      = "#BCBD22"   # 橄榄黄色
# )
# 
# # 图像文本与布局参数
# plot_title    <- "Number of nuclei (mean genes per nucleus)"
# text_nudge_y  <- 200                     # 柱状图上方文本的偏移量
# bar_width     <- 0.8                     # 柱子的宽度
# 
# # 输出参数
# output_dir    <- "./results/final_fig/"  # 输出文件夹
# file_prefix   <- "DimPlot_and_Bar_"      # 输出文件前缀
# save_width    <- 12                      # 画布宽度
# save_height   <- 4                       # 画布高度
# 
# 

## 1 Parameters(Condition) --------------------------------------------------------------
# 数据与分组参数
seurat_obj  <- tmp_obj                   # 输入的 Seurat 对象
group_col   <- "Condition"               # 【修改点 1】: 将分组列名改为 Condition
type_col    <- "CellType1"               # 细胞类型列名 (预留，以防后续扩展)

# 颜色参数 (疾病状态定制配色)
# 【修改点 2】: 包含 7 种颜色，Control 采用经典灰色作为基准底色
manual_colors <- c(
  "Control" = "#B3B3B3",  # 中性灰 (作为健康对照组的视觉基准)
  "AD"      = "#D53E4F",  # 绯红色 (阿尔茨海默病)
  "PDD"     = "#3288BD",  # 海蓝色 (帕金森病痴呆)
  "DLB"     = "#5E4FA2",  # 靛紫色 (路易体痴呆)
  "MS"      = "#F46D43",  # 亮橙色 (多发性硬化症)
  "Stroke"  = "#8C510A",  # 棕褐色 (中风)
  "COVID"   = "#C51B7D"   # 品红色 (新冠)
)

# 图像文本与布局参数
plot_title    <- "Number of nuclei (mean genes per nucleus)"
text_nudge_y  <- 200                     # 柱状图上方文本的偏移量 (若某疾病细胞数极大，可适当调大)
bar_width     <- 0.8                     # 柱子的宽度

# 输出参数
output_dir    <- "./results/final_fig/"  # 输出文件夹
file_prefix   <- "DimPlot_and_Bar_"      # 输出文件前缀
save_width    <- 12                      # 画布宽度
save_height   <- 4                       # 画布高度
## 2 Process Data ------------------------------------------------------------
# 自动创建输出目录（如果不存在）
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 统计每个组的细胞数和平均表达基因数
plot_df_num <- seurat_obj@meta.data %>%
  group_by(Group = .data[[group_col]]) %>%
  summarise(
    Count = n(), 
    Mean_Genes = round(mean(nFeature_RNA, na.rm = TRUE))
  ) %>%
  mutate(Label = paste0(Count, " (", Mean_Genes, ")")) %>%
  # 按数量排序
  mutate(Group = reorder(Group, Count)) %>%
  # 强制将 "Unknown" 移动到第一个 Level (如果有的话)
  mutate(Group = fct_relevel(Group, intersect("Unknown", levels(Group))))


## 3 Generate Plots ----------------------------------------------------------
# 3.1 降维图 (UMAP/tSNE) - 常规
p_dim <- DimPlot(seurat_obj, group.by = group_col, cols = manual_colors, pt.size = 0.1)

# 3.2 柱状统计图
p_stat <- ggplot(plot_df_num, aes(x = Group, y = Count, fill = Group)) +
  geom_col(color = "black", width = bar_width) +
  coord_flip() +
  geom_text(aes(label = Label), hjust = 0, nudge_y = text_nudge_y, size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
  scale_fill_manual(values = manual_colors) +
  labs(x = NULL, y = NULL, title = plot_title) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, size = 14, face = "plain"),
    axis.text       = element_text(size = 12, color = "black")
  )

# 3.3 拼图 (常规降维图 + 柱状图)
p_combined <- (p_dim + p_stat) + 
  plot_layout(widths = c(1, 1))   # 强制两张图的宽度比例为 1:1

# 3.4 降维图 (UMAP/tSNE) - 分面 (Split)
# 移除图例使画面更清爽，因为每个子图的标题已经标明了分组
num_cols <- length(unique(seurat_obj@meta.data[[group_col]]))

# 画图代码
p_split <- DimPlot(seurat_obj, group.by = group_col, cols = manual_colors, split.by = group_col, pt.size = 0.3) +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position = "right", 
    axis.text.y     = element_blank(), 
    axis.text.x     = element_blank(), 
    axis.ticks      = element_blank(),
    axis.line       = element_blank(),
    strip.text      = element_text(color = "white", face = "bold", size = 10) 
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  facet_wrap2(
    vars(!!sym(group_col)),
    nrow = 1,          # 依然强制保持 1 行排版
    ncol = num_cols,   # 【核心修改点】：传入上面动态计算出来的列数变量
    strip = strip_themed(
      background_x = elem_list_rect(fill = manual_colors) 
    )
  )

## 4 Save Results ------------------------------------------------------------
# 4.1 保存拼图 (DimPlot + Barplot)
output_file_combined <- paste0(output_dir, file_prefix, group_col, "_Combined.pdf")
ggsave(
  filename = output_file_combined, 
  plot     = p_combined, 
  width    = save_width, 
  height   = save_height
)
message("Combined Plot saved to: ", output_file_combined)

# 4.2 保存分面图 (Split DimPlot)
output_file_split <- paste0(output_dir, file_prefix, group_col, "_Split.pdf")
# 动态计算分面图需要的画布宽度 (每个分组给大约 3 到 4 的宽度，确保不拥挤)
num_groups <- length(unique(plot_df_num$Group))
save_width_split <- max(12, num_groups * 3.5) 

ggsave(
  filename = output_file_split, 
  plot     = p_split, 
  width    = save_width_split, 
  height   = save_height 
)
message("Split Plot saved to: ", output_file_split)
