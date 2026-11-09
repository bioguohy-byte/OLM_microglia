library(openxlsx)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rstatix)

# MODULE : Group Comparisons Boxplot -------------------

## 1 Parameters  ---------------------------------------------------
input_filepath  <- "data/processed_data/sheet6.xlsx"
plot_name       <- "sheet6"
output_filepath <- paste0("results/figure/", plot_name, "_boxplot.pdf")

var_group_name  <- "Group"
var_value_name  <- "Value"

mycolor      <- c("#0077c1", "#00a99e", "#6bc72b")
label_x      <- ""
label_y      <- "value"
plot_margins <- c(0.2, 0.2, 0.2, 0.2) 

stat_method    <- "t.test"
# ref_group      <- "SP8-aging"
# my_comparisons <- list(
#   c("SP8-aging", "C57-young"),
#   c("SP8-aging", "samp8-young")
# )
my_symnum      <- list(
  cutpoints = c(0, 0.001, 0.01, 0.05, 1),
  symbols   = c("***", "**", "*", "ns")
)

out_width  <- 6
out_height <- 2

## 2 Process Data ------------------------------------------------------------
df <- read.xlsx(input_filepath)

dt <- melt(df,
           id.vars       = NULL,
           measure.vars  = colnames(df),
           variable.name = var_group_name,
           value.name    = var_value_name)

### 2.5 Calculate Statistics  ---------------------------------------
# 使用 rstatix 提前计算 t.test 和 y 轴坐标位置
stat.test <- dt %>%
  t_test(Value ~ Group, 
         # comparisons = my_comparisons
  ) %>%
  add_xy_position(x = "Group")

# 将 p 值转换为你自定义的星号
stat.test$p.signif <- symnum(stat.test$p, 
                             corr = FALSE, na = FALSE, 
                             cutpoints = my_symnum$cutpoints, 
                             symbols = my_symnum$symbols)


## 3 Generate Plots ----------------------------------------------------------
mytheme2 <- theme_classic() +
  theme(
    text              = element_text(family = "sans", colour = "black", size = 12),
    axis.line         = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks.length = unit(1.5, units = "mm"),
    legend.position   = "none",
    plot.margin       = unit(x = plot_margins, units = "inches")
  )

# 注意：这里去掉了 orientation = "horizontal"
p_base <- ggboxplot(dt, 
                    x     = var_group_name, 
                    y     = var_value_name,
                    color = var_group_name,
                    width = 0.4,
                    fill  = var_group_name,
                    add   = "jitter",
                    add.params = list(size = 0.8)) + 
  scale_colour_manual(values = alpha(mycolor, 1)) +
  scale_fill_manual(values = alpha(mycolor, 0.3)) +
  labs(x = label_x, y = label_y) +
  coord_flip() + # <--- 使用 coord_flip() 翻转整个画布
  mytheme2

# 使用 stat_pvalue_manual 替代 stat_compare_means
p_final <- p_base + stat_pvalue_manual(
  stat.test,
  label        = "p.signif",
  tip.length   = 0.02,
  bracket.size = 0.5,
  angle        = -90,  # 让星号顺时针旋转90度，与括号平行
  coord.flip   = TRUE  # <--- 关键参数：告诉它画布被翻转了，星号和括号会自动对齐
)

print(p_final)

## 4 Save Results ------------------------------------------------------------
out_dir <- dirname(output_filepath)
if(!dir.exists(out_dir) && out_dir != ".") {
  dir.create(out_dir, recursive = TRUE)
}

ggsave(filename = output_filepath, 
       plot     = p_final, 
       width    = out_width, 
       height   = out_height)

