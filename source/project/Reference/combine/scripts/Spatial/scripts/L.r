# MODULE : Spatial L-Function (Co-localization) Analysis -------------------

## 1 Parameters --------------------------------------------------------------
# 加载依赖包
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(spatstat)

# 输入输出路径
file_in_rds  <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/merfish_s2r1_v2.rds")
dir_out_fig  <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/results/figure/")
file_out_pdf <- paste0(dir_out_fig, "s2r1_OLG_vs_other_MG_subclusters_L_function_plot_single_row.pdf")

# 空间坐标截取范围
coord_filter <- c(x_min = -3000, x_max = -1500, y_min = -8000, y_max = -5300)

# 分析的细胞类型配置
center_cell     <- "OLG"
target_cells    <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")
min_cell_count  <- 10  # 纳入分析的最小细胞数量阈值

# 空间统计学模拟参数 (L-function)
sim_r_max <- 250   # 最大半径 (μm)
sim_step  <- 1     # 步长
sim_nsim  <- 200   # 随机置换模拟次数

# 绘图配置
plot_r_min <- 5    # 画图时截去的最小半径
plot_ncol  <- length(target_cells) # 自动根据目标细胞数量设定列数 (保证单行)
plot_width <- 12
plot_height<- 3.5

# 细胞类型颜色字典 (包含所有目标细胞)
color_dict <- c(
  "ST_MG0" = "#E64B35", 
  "ST_MG1" = "#4DBBD5", 
  "ST_MG2" = "#00A087", 
  "ST_MG3" = "#3C5488", 
  "ST_MG4" = "#F39B7F", 
  "ST_MG5" = "#7E6148"
)


## 2 Process Data ------------------------------------------------------------
# 2.1 读取与预处理数据
obj <- readRDS(file_in_rds)
df <- obj@meta.data %>% 
  select(rotated_x, rotated_y, CellType) %>%
  filter(rotated_x >= coord_filter["x_min"] & rotated_x <= coord_filter["x_max"] & 
           rotated_y >= coord_filter["y_min"] & rotated_y <= coord_filter["y_max"])

# 确保坐标是数值型
df$rotated_x <- as.numeric(df$rotated_x)
df$rotated_y <- as.numeric(df$rotated_y)

# 2.2 创建 ppp 对象
pp <- ppp(df$rotated_x, df$rotated_y, 
          window = owin(range(df$rotated_x), range(df$rotated_y)), 
          marks = factor(df$CellType))

# 2.3 定义优化后的 L 函数计算模块
centered_L <- function(pp, i, j, r_max, step, nsim) {
  r_seq <- seq(0, r_max, by = step)
  
  # envelope 内部会自动计算观测值($obs)和模拟值($mmean, $lo, $hi)
  L_sim <- envelope(pp, Lcross, i = i, j = j, nsim = nsim, r = r_seq, 
                    correction = "border", simulate = expression(rlabel(pp)),
                    verbose = FALSE) 
  
  if (is.null(L_sim$lo) || is.null(L_sim$hi) || is.null(L_sim$mmean)) {
    warning(paste("Unable to generate simulation envelope for", i, "vs", j))
    return(NULL)
  }
  
  # 直接利用 envelope 返回的 $obs 减去模拟均值 $mmean，完美对齐向量
  data.frame(
    r  = r_seq, 
    L  = L_sim$obs - L_sim$mmean, 
    lo = L_sim$lo - L_sim$mmean, 
    hi = L_sim$hi - L_sim$mmean
  )
}

# 2.4 批量计算细胞对的空间分布
results <- list()
for (target in target_cells) {
  count_center <- sum(df$CellType == center_cell)
  count_target <- sum(df$CellType == target)
  
  pair_name <- paste(center_cell, target, sep = "_vs_")
  print(paste(sprintf("[%s]", pair_name), "- Center Count:", count_center, "| Target Count:", count_target))
  
  if (count_center >= min_cell_count && count_target >= min_cell_count) {
    tryCatch({
      res <- centered_L(pp, center_cell, target, r_max = sim_r_max, step = sim_step, nsim = sim_nsim)
      if (!is.null(res)) {
        res$cell_type1 <- center_cell
        res$cell_type2 <- target
        results[[pair_name]] <- res
        print(paste("   -> Successfully analyzed:", pair_name))
      }
    }, error = function(e) {
      print(paste("   -> Error analyzing:", pair_name, "-", e$message))
    })
  } else {
    print(paste("   -> Skipped due to insufficient cell count"))
  }
}

# 2.5 整合结果为 DataFrame
df_combined <- bind_rows(results)


## 3 Generate Plots ----------------------------------------------------------
# 过滤掉极小半径的伪影噪音
df_filtered <- df_combined %>% filter(r >= plot_r_min)

# 绘制中心化 L 函数包络图 (采用宽格式绘图逻辑)
p_l_func <- ggplot(df_filtered, aes(x = r)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = cell_type2), alpha = 0.3) +
  geom_line(aes(y = L, color = cell_type2), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ cell_type2, scales = "fixed", ncol = plot_ncol) + # <--- 仅将 "free_y" 修改为 "fixed"
  labs(x = "Distance (μm)", y = "Centered L(r)") +
  theme_minimal() +
  theme(
    legend.position = "bottom", # 保持原样不变
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  scale_color_manual(values = color_dict) +
  scale_fill_manual(values = color_dict)


## 4 Save Results ------------------------------------------------------------
# 检查并创建图片输出目录
if(!dir.exists(dir_out_fig)) dir.create(dir_out_fig, recursive = TRUE)

# 保存最终图形
ggsave(filename = file_out_pdf, plot = p_l_func, 
       width = plot_width, height = plot_height, units = "in",  
       device = cairo_pdf, dpi = 300)

print(paste("Plot saved to:", file_out_pdf))