setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/"))
rm(list=ls())

# 0. Libraries -----------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Seurat)
library(Matrix)
library(cowplot)
library(reshape2)
library(SingleCellExperiment)
library(umap)
library(ggthemes)
library(dplyr)
library(patchwork)
library(ggsci)
library(scales)
library(tidyverse)
library(magick)
library(magrittr)
library(devtools)
library(ComplexHeatmap)
library(ggplotify)
library(circlize)
library(pheatmap)
library(harmony)
library(qs)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)

# ==============================
# 1. 设置通用变量 (只需修改这里)
# ==============================
target_meta <- "Study"  # 想要分析的 metadata 列名
output_dir  <- Sys.getenv("MG_EXTERNAL_002") # 输出路径

# 定义颜色 (保持不变)
cluster_colors <- c(
  "0"  = "#4E79A7", "1"  = "#807DBA", "2"  = "#B3DE69", "3"  = "#A0CBE8",
  "4"  = "#FDB462", "5"  = "#FCCDE5", "6"  = "#FB8072", "7"  = "#BEBADA",
  "8"  = "#D9D9D9", "9"  = "#FCBBA1", "10" = "#BC80BD", "11" = "#E41A1C",
  "12" = "#FFED6F"
)

# ==============================
# 2. DimPlot 绘图 (p & p2)
# ==============================

# p: 整体分布图
p <- DimPlot(tmp_obj, 
             group.by = target_meta,     # 这里的 Seurat 函数直接接受字符串变量
             reduction = "umap", 
             label = FALSE, 
             repel = TRUE) + 
  labs(title = target_meta,              # 标题自动跟随变量
       x = "UMAP_1", y = "UMAP_2") + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(), 
    axis.text = element_text(color = "black"), 
    panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"), 
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    legend.text = element_text(size = 10)
  )
print(p)

# p2: 拆分视图
p2 <- DimPlot(tmp_obj, 
              group.by = 'integrated_snn_res.0.25', 
              split.by = target_meta,    # split.by 也接受字符串变量
              cols = cluster_colors,
              reduction = "umap", 
              label = FALSE, 
              repel = TRUE) + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(), 
    axis.text = element_text(color = "black"), 
    panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"), 
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    legend.text = element_text(size = 10)
  )
print(p2)

# ==============================
# 3. 统计绘图 (p3)
# ==============================

# 准备数据
plot_data_stats <- tmp_obj@meta.data %>%
  group_by(!!sym(target_meta)) %>%  # 使用 !!sym() 将字符串转为列名进行分组
  summarise(
    Count = n(),
    MeanGenes = round(mean(nFeature_RNA), 0)
  ) %>%
  arrange(desc(Count)) %>%
  mutate(Label = paste0(Count, " (", MeanGenes, ")"))

# 动态设置 Factor Level 以实现排序 (从上到下：多到少)
# 注意：这里使用了 base R 的 [[]] 写法来动态访问列
plot_data_stats[[target_meta]] <- factor(plot_data_stats[[target_meta]], 
                                         levels = rev(plot_data_stats[[target_meta]]))

# 绘图
p3 <- ggplot(plot_data_stats, aes(x = Count, 
                                  y = !!sym(target_meta), # y轴使用动态变量
                                  fill = !!sym(target_meta))) + # 颜色使用动态变量
  geom_col(color = "black", width = 0.8, show.legend = FALSE) +
  geom_text(aes(label = Label, x = 0), 
            hjust = 0, 
            nudge_x = 10, 
            size = 5, 
            color = "black") +
  theme_classic() +
  labs(x = "Number of nuclei (mean genes per nucleus)", 
       y = NULL,
       title = paste("Statistics by", target_meta)) + # 动态标题
  theme(
    axis.line.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    axis.text.y = element_text(size = 14, color = "black", face = "plain"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3)))

# 动态生成文件名
save_path <- paste0(output_dir, target_meta, '_stats.png')
ggsave(p3, file = save_path)
print(p3)

# ==============================
# 4. 比例组成图 (p4)
# ==============================

# 准备数据
plot_data_comp <- tmp_obj@meta.data %>%
  # 两个分组变量：一个是动态的，一个是 Cluster
  group_by(!!sym(target_meta), integrated_snn_res.0.25) %>% 
  summarise(Count = n(), .groups = 'drop') %>%
  # 计算百分比
  group_by(!!sym(target_meta)) %>% 
  mutate(Percentage = Count / sum(Count) * 100)

# 绘图 p4 (堆叠柱状图)
p4 <- ggplot(plot_data_comp, aes(x = !!sym(target_meta), # x轴动态
                                 y = Percentage, 
                                 fill = integrated_snn_res.0.25)) +
  geom_bar(stat = "identity", color = "black", width = 0.7, position = "stack") +
  scale_fill_manual(values = cluster_colors, name = "Clusters") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  labs(y = "Composition Percentage (%)", 
       x = target_meta, 
       title = paste("Cell Composition by", target_meta))

print(p4)

# 绘图 (Facet 分面图)
# 注意：facet_wrap 里的公式需要用 vars() 或者 as.formula()
p5 <- ggplot(plot_data_comp, aes(x = integrated_snn_res.0.25, 
                                 y = Percentage, 
                                 fill = integrated_snn_res.0.25)) +
  geom_bar(stat = "identity", color = "black") +
  # 动态分面：使用 vars(!!sym(target_meta))
  facet_wrap(vars(!!sym(target_meta)), scales = "fixed") + 
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  labs(x = "Cluster", y = paste("Percentage within", target_meta, "(%)"))

print(p5)



#---------------------------
p <- DimPlot(tmp_obj, 
             group.by = 'Region', 
             reduction = "umap",    # 使用 UMAP 降维结果
             label = TRUE,          # 在图上显示亚群名称
             label.size = 0,        # 标签字体大小
             repel = TRUE) +        # 自动排布标签防止重叠
  labs(title = "Region", # 图片标题
       x = "UMAP_1", y = "UMAP_2") +                 # 坐标轴名称
  theme_bw() +                                       # 使用白底黑框主题
  theme(
    panel.grid = element_blank(),                    # 移除背景网格线，使图面整洁
    axis.text = element_text(color = "black"),       # 坐标轴文字设为黑色
    panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"), # 加深外框
    plot.title = element_text(hjust = 0.5, face = "bold"), # 标题居中并加粗
    legend.text = element_text(size = 10)            # 图例文本大小
  )
p


cluster_colors <- c(
  "0"  = "#4E79A7", # HOM-1 (稳态 1)
  "1"  = "#807DBA", # ARM-1 (激活 1)
  "2"  = "#B3DE69", # MAM (髓鞘相关)
  "3"  = "#A0CBE8", # HOM-2 (稳态 2)
  "4"  = "#FDB462", # PIM (促炎)
  "5"  = "#FCCDE5", # SYM (突触调节)
  "6"  = "#FB8072", # DAM-1 (疾病相关 1)
  "7"  = "#BEBADA", # ARM-2 (激活 2)
  "8"  = "#D9D9D9", # SRM (应激反应)
  "9"  = "#FCBBA1", # DAM-2 (疾病相关 2)
  "10" = "#BC80BD", # BAM (边界相关巨噬细胞)
  "11" = "#E41A1C", # CRM (趋化因子富集)
  "12" = "#FFED6F"  # NAM (神经元相关)
)


p2 <- DimPlot(tmp_obj, 
              group.by = 'integrated_snn_res.0.25', 
              split.by = 'Region',
              cols = cluster_colors,
              reduction = "umap",    # 使用 UMAP 降维结果
              label = FALSE,          # 在图上显示亚群名称
              repel = TRUE) +        # 自动排布标签防止重叠
  theme_bw() +                                       # 使用白底黑框主题
  theme(
    panel.grid = element_blank(),                    # 移除背景网格线，使图面整洁
    axis.text = element_text(color = "black"),       # 坐标轴文字设为黑色
    panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"), # 加深外框
    plot.title = element_text(hjust = 0.5, face = "bold"), # 标题居中并加粗
    legend.text = element_text(size = 10)            # 图例文本大小
  )
p2


#----------------------
library(ggplot2)
library(dplyr)

# 1. 准备数据
plot_data <- tmp_obj@meta.data %>%
  # 关键修改：统计每个 Region 下每个 Cluster 的数量
  group_by(Region, integrated_snn_res.0.25) %>% 
  summarise(Count = n(), .groups = 'drop') %>%
  # 关键修改：按 Region 分组，确保同一个 Region 内的比例相加为 100%
  group_by(Region) %>% 
  mutate(Percentage = Count / sum(Count) * 100)

# 2. 绘图（推荐使用堆叠柱状图 Stacked Bar Chart）
p3 <- ggplot(plot_data, aes(x = Region, y = Percentage, fill = integrated_snn_res.0.25)) +
  # position = "stack" 确保柱子堆叠起来，总高度为 100
  geom_bar(stat = "identity", color = "black", width = 0.7, position = "stack") +
  
  # 使用你自定义的颜色
  scale_fill_manual(values = cluster_colors, name = "Clusters") +
  
  # 样式调整
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right" # 在右侧显示图例
  ) +
  labs(y = "Composition Percentage (%)", x = "Brain Region", title = "Cell Composition by Region")

print(p3)


ggplot(plot_data, aes(x = integrated_snn_res.0.25, y = Percentage, fill = integrated_snn_res.0.25)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~ Region, scales = "fixed") + 
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  labs(x = "Cluster", y = "Percentage within Region (%)")




library(ggplot2)
library(dplyr)

# 1. 准备统计数据
# 假设你的基因数存放在 nFeature_RNA 列中
plot_data <- tmp_obj@meta.data %>%
  group_by(Region) %>%
  summarise(
    Count = n(),
    MeanGenes = round(mean(nFeature_RNA), 0)
  ) %>%
  # 按细胞数量从大到小排序
  arrange(desc(Count)) %>%
  # 创建标签文字：例如 "32301 (1171)"
  mutate(Label = paste0(Count, " (", MeanGenes, ")"))

# 为了让图中“从上到下”是“从多到少”，需要对 Factor Level 进行反向排序
plot_data$Region <- factor(plot_data$Region, levels = rev(plot_data$Region))

# 2. 绘图
p4 <- ggplot(plot_data, aes(x = Count, y = Region, fill = Region)) +
  # 绘制横向柱状图
  geom_col(color = "black", width = 0.8, show.legend = FALSE) +
  
  # 添加柱子上的文字标注
  # hjust = 0 表示文字左对齐，nudge_x 可以微调文字与柱子起始位置的距离
  geom_text(aes(label = Label, x=0), 
            hjust = 0, 
            nudge_x = 10, # 根据你的数据量级调整这个值，确保文字不重叠
            size = 5, 
            color = "black") +
  
  # 使用自定义颜色（如果有的话）或默认配色
  # scale_fill_manual(values = region_colors) + 
  theme_classic() +
  
  # 调整坐标轴和布局
  labs(x = "Number of nuclei (mean genes per nucleus)", 
       y = NULL) +
  theme(
    axis.line.y = element_blank(), # 隐藏 Y 轴线
    axis.ticks.y = element_blank(), # 隐藏 Y 轴刻度
    axis.text.y = element_text(size = 14, color = "black", face = "plain"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15))
  ) +
  # 稍微拉长 X 轴范围，防止标签文字被切掉
  scale_x_continuous(expand = expansion(mult = c(0, 0.3)))

ggsave(p4,file = Sys.getenv("MG_EXTERNAL_003"))
print(p4)




#---------------------------
p <- DimPlot(tmp_obj, 
             group.by = 'Condition', 
             reduction = "umap",    # 使用 UMAP 降维结果
             label = TRUE,          # 在图上显示亚群名称
             label.size = 0,        # 标签字体大小
             repel = TRUE) +        # 自动排布标签防止重叠
  labs(title = "Condition", # 图片标题
       x = "UMAP_1", y = "UMAP_2") +                 # 坐标轴名称
  theme_bw() +                                       # 使用白底黑框主题
  theme(
    panel.grid = element_blank(),                    # 移除背景网格线，使图面整洁
    axis.text = element_text(color = "black"),       # 坐标轴文字设为黑色
    panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"), # 加深外框
    plot.title = element_text(hjust = 0.5, face = "bold"), # 标题居中并加粗
    legend.text = element_text(size = 10)            # 图例文本大小
  )
p


cluster_colors <- c(
  "0"  = "#4E79A7", # HOM-1 (稳态 1)
  "1"  = "#807DBA", # ARM-1 (激活 1)
  "2"  = "#B3DE69", # MAM (髓鞘相关)
  "3"  = "#A0CBE8", # HOM-2 (稳态 2)
  "4"  = "#FDB462", # PIM (促炎)
  "5"  = "#FCCDE5", # SYM (突触调节)
  "6"  = "#FB8072", # DAM-1 (疾病相关 1)
  "7"  = "#BEBADA", # ARM-2 (激活 2)
  "8"  = "#D9D9D9", # SRM (应激反应)
  "9"  = "#FCBBA1", # DAM-2 (疾病相关 2)
  "10" = "#BC80BD", # BAM (边界相关巨噬细胞)
  "11" = "#E41A1C", # CRM (趋化因子富集)
  "12" = "#FFED6F"  # NAM (神经元相关)
)


p2 <- DimPlot(tmp_obj, 
              group.by = 'integrated_snn_res.0.25', 
              split.by = 'Condition',
              cols = cluster_colors,
              reduction = "umap",    # 使用 UMAP 降维结果
              label = FALSE,          # 在图上显示亚群名称
              repel = TRUE) +        # 自动排布标签防止重叠
  theme_bw() +                                       # 使用白底黑框主题
  theme(
    panel.grid = element_blank(),                    # 移除背景网格线，使图面整洁
    axis.text = element_text(color = "black"),       # 坐标轴文字设为黑色
    panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"), # 加深外框
    plot.title = element_text(hjust = 0.5, face = "bold"), # 标题居中并加粗
    legend.text = element_text(size = 10)            # 图例文本大小
  )
p2

library(ggplot2)
library(dplyr)

# 1. 准备统计数据
# 假设你的基因数存放在 nFeature_RNA 列中
plot_data <- tmp_obj@meta.data %>%
  group_by(Condition) %>%
  summarise(
    Count = n(),
    MeanGenes = round(mean(nFeature_RNA), 0)
  ) %>%
  # 按细胞数量从大到小排序
  arrange(desc(Count)) %>%
  # 创建标签文字：例如 "32301 (1171)"
  mutate(Label = paste0(Count, " (", MeanGenes, ")"))

# 为了让图中“从上到下”是“从多到少”，需要对 Factor Level 进行反向排序
plot_data$Condition <- factor(plot_data$Condition, levels = rev(plot_data$Condition))

# 2. 绘图
p4 <- ggplot(plot_data, aes(x = Count, y = Condition, fill = Condition)) +
  # 绘制横向柱状图
  geom_col(color = "black", width = 0.8, show.legend = FALSE) +
  
  # 添加柱子上的文字标注
  # hjust = 0 表示文字左对齐，nudge_x 可以微调文字与柱子起始位置的距离
  geom_text(aes(label = Label, x=0), 
            hjust = 0, 
            nudge_x = 10, # 根据你的数据量级调整这个值，确保文字不重叠
            size = 5, 
            color = "black") +
  
  # 使用自定义颜色（如果有的话）或默认配色
  # scale_fill_manual(values = Condition_colors) + 
  theme_classic() +
  
  # 调整坐标轴和布局
  labs(x = "Number of nuclei (mean genes per nucleus)", 
       y = NULL) +
  theme(
    axis.line.y = element_blank(), # 隐藏 Y 轴线
    axis.ticks.y = element_blank(), # 隐藏 Y 轴刻度
    axis.text.y = element_text(size = 14, color = "black", face = "plain"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 12, margin = margin(t = 15))
  ) +
  # 稍微拉长 X 轴范围，防止标签文字被切掉
  scale_x_continuous(expand = expansion(mult = c(0, 0.3)))

ggsave(p4,file = Sys.getenv("MG_EXTERNAL_004"))
print(p4)


library(ggplot2)
library(dplyr)

# 1. 准备数据
plot_data <- tmp_obj@meta.data %>%
  # 关键修改：统计每个 Condition 下每个 Cluster 的数量
  group_by(Condition, integrated_snn_res.0.25) %>% 
  summarise(Count = n(), .groups = 'drop') %>%
  # 关键修改：按 Condition 分组，确保同一个 Condition 内的比例相加为 100%
  group_by(Condition) %>% 
  mutate(Percentage = Count / sum(Count) * 100)

# 2. 绘图（推荐使用堆叠柱状图 Stacked Bar Chart）
p3 <- ggplot(plot_data, aes(x = Condition, y = Percentage, fill = integrated_snn_res.0.25)) +
  # position = "stack" 确保柱子堆叠起来，总高度为 100
  geom_bar(stat = "identity", color = "black", width = 0.7, position = "stack") +
  
  # 使用你自定义的颜色
  scale_fill_manual(values = cluster_colors, name = "Clusters") +
  
  # 样式调整
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right" # 在右侧显示图例
  ) +
  labs(y = "Composition Percentage (%)", x = "Condition", title = "Cell Composition by Condition")

print(p3)


ggplot(plot_data, aes(x = integrated_snn_res.0.25, y = Percentage, fill = integrated_snn_res.0.25)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~ Condition, scales = "fixed") + 
  scale_fill_manual(values = cluster_colors) +
  theme_classic() +
  labs(x = "Cluster", y = "Percentage within Condition (%)")

