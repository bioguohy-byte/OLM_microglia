library(openxlsx)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rstatix)

# MODULE : Group Comparisons Violin Plot -------------------

## 1 Parameters  ---------------------------------------------------
input_filepath  <- "data/processed_data/sheet6.xlsx"
plot_name       <- "sheet6_violin"
output_filepath <- paste0("results/figure/", plot_name, ".pdf")

var_group_name  <- "Group"
var_value_name  <- "Value"

mycolor      <- c("#0077c1", "#00a99e", "#6bc72b")
label_x      <- ""
label_y      <- "value"
plot_margins <- c(0.2, 0.2, 0.2, 0.2) 

stat_method    <- "t.test"
my_symnum      <- list(
  cutpoints = c(0, 0.001, 0.01, 0.05, 1),
  symbols   = c("***", "**", "*", "ns")
)

# 【修改1：调整画板宽高】竖着的图需要更高一点的空间放显著性括号
out_width  <- 2   
out_height <- 5   

## 2 Process Data ------------------------------------------------------------
df <- read.xlsx(input_filepath)

dt <- melt(df,
           id.vars       = NULL,
           measure.vars  = colnames(df),
           variable.name = var_group_name,
           value.name    = var_value_name)

### 2.5 Calculate Statistics  ---------------------------------------
stat.test <- dt %>%
  t_test(Value ~ Group) %>%
  add_xy_position(x = "Group")

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
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
    legend.position   = "none",
    plot.margin       = unit(x = plot_margins, units = "inches")
  )

p_base <- ggviolin(dt, 
                   x     = var_group_name, 
                   y     = var_value_name,
                   color = var_group_name,
                   width = 0.6,               # 小提琴图稍微放宽一点更好看
                   fill  = var_group_name,
                   add   = "jitter",          # 依然保留散点
                   add.params = list(
                     shape = 21,            # 设置点的形状为 21 号（可描边的圆）
                     size  = 1.5,           # 点的大小
                     color = "black",       # 点的描边颜色（黑色）
                     fill  = var_group_name # 点的填充颜色（映射到分组，与小提琴图一致）
                   )) + 
  scale_colour_manual(values = alpha(mycolor, 1)) +
  scale_fill_manual(values = alpha(mycolor, 0.3)) +
  labs(x = label_x, y = label_y) +
  mytheme2


# 【修改3：清理掉用于翻转的参数】
p_final <- p_base + stat_pvalue_manual(
  stat.test,
  label        = "p.signif",
  tip.length   = 0.02,
  bracket.size = 0.5
  # 删除了 angle = -90 和 coord.flip = TRUE
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