#!/usr/bin/env Rscript
# Figure-ordered code from the reviewed manuscript analysis package.
# Each panel contains its original quoted R code.
# Source this chapter, then use mg_list_panels() or mg_show_panel("PANEL_ID").
# Source paths in panel records identify the original analysis scripts.
.mg_panels <- new.env(parent=emptyenv())
.mg_order <- character()
mg_register_panel <- function(id, code) {
  if (!id %in% .mg_order) .mg_order <<- c(.mg_order,id)
  assign(id, code, envir=.mg_panels)
}
mg_list_panels <- function() .mg_order
mg_panel_code <- function(id) {
  if (!exists(id,envir=.mg_panels,inherits=FALSE)) stop("Unknown panel: ",id)
  get(id,envir=.mg_panels,inherits=FALSE)
}
mg_show_panel <- function(id) {
  blocks <- mg_panel_code(id)
  if (is.null(blocks)) {cat(id, ": No R block is included for this panel.\n");return(invisible(NULL))}
  for (block in blocks) {
    cat("# Source: ",block$source," [actual ",block$ranges,"; requested ",block$requested_ranges,"]\n# Role: ",block$role,"\n",sep="")
    cat(paste(deparse(block$code,width.cutoff=120L),collapse="\n"),"\n\n")
  }
  invisible(blocks)
}


# Figure 6

# F6A  GFP-HMC3 myelin conditioning and corpus-callosum transplantation schematic
# Kind: experimental_schematic
mg_register_panel("F6A", NULL)

# F6B  GFP-HMC3 fluorescence and viral transduction violin plot
# Kind: experimental_quantification
mg_register_panel("F6B", list(
  list(source="source/project/Reference/combine/scripts/Boxplot/scripts/sheet5_vlnplot.r", ranges="1-105", requested_ranges="1-105", role="HMC3 GFP-positive fraction; working directory: Reference/combine/scripts/Boxplot", code=quote({
library(openxlsx)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rstatix)

# MODULE : Group Comparisons Violin Plot -------------------

## 1 Parameters  ---------------------------------------------------
input_filepath  <- "data/processed_data/sheet5.xlsx"
plot_name       <- "sheet5_violin"
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

# ### 2.5 Calculate Statistics  ---------------------------------------
# stat.test <- dt %>%
#   t_test(Value ~ Group) %>%
#   add_xy_position(x = "Group")
# 
# stat.test$p.signif <- symnum(stat.test$p, 
#                              corr = FALSE, na = FALSE, 
#                              cutpoints = my_symnum$cutpoints, 
#                              symbols = my_symnum$symbols)


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


p_final <- p_base

# # 【修改3：清理掉用于翻转的参数】
# p_final <- p_base + stat_pvalue_manual(
#   stat.test,
#   label        = "p.signif",
#   tip.length   = 0.02,
#   bracket.size = 0.5
#   # 删除了 angle = -90 和 coord.flip = TRUE
# )

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
  }))
))

# F6C  Olig2/GFP/DAPI corpus-callosum images after HMC3-control or myelin-exposed HMC3 graft
# Kind: experimental_micrographs
mg_register_panel("F6C", NULL)

# F6D  Olig2-positive OLG fluorescent area fraction after graft
# Kind: experimental_quantification
mg_register_panel("F6D", list(
  list(source="source/recovered/experimental/F6D_OLG_area_fraction.R", ranges="1-5;10-12;14-15;17-20;22-26;29-30;33-33;35-39;42-44;46-49;53-62;64-76;80-86;88-88;91-94;96-99", requested_ranges="1-99", role="Original two-group statistical plot", code=quote({
library(openxlsx)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rstatix)

input_filepath <- file.path(Sys.getenv("MG_PLOT_INPUT_ROOT", "data/plot_inputs"), "sp8_m0.xlsx")
plot_name       <- "sp8_m0"
output_filepath <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "F6D", paste0(plot_name, ".pdf"))

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

out_width  <- 2   
out_height <- 3   

df <- read.xlsx(input_filepath)

dt <- melt(df,
           id.vars       = NULL,
           measure.vars  = colnames(df),
           variable.name = var_group_name,
           value.name    = var_value_name)

stat.test <- dt %>%
  t_test(Value ~ Group) %>%
  add_xy_position(x = "Group")

stat.test$p.signif <- symnum(stat.test$p,
                             corr = FALSE, na = FALSE,
                             cutpoints = my_symnum$cutpoints,
                             symbols = my_symnum$symbols)

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
  # coord_flip() + # <--- 使用 coord_flip() 翻转整个画布
  mytheme2

p_final <- p_base + stat_pvalue_manual(
  stat.test,
  label        = "p.signif",
  tip.length   = 0.02,
  bracket.size = 0.5
  # 删除了 angle = -90 和 coord.flip = TRUE
)

print(p_final)

out_dir <- dirname(output_filepath)
if(!dir.exists(out_dir) && out_dir != ".") {
  dir.create(out_dir, recursive = TRUE)
}

ggsave(filename = output_filepath, 
       plot     = p_final, 
       width    = out_width, 
       height   = out_height)
  }))
))

# F6E  Young C57/young SAMP8/aging SAMP8 schematic and TEM images
# Kind: experimental_schematic_and_micrographs
mg_register_panel("F6E", NULL)

# F6F  G-ratios comparing young C57, young SAMP8 and aging SAMP8
# Kind: experimental_quantification
mg_register_panel("F6F", list(
  list(source="source/recovered/experimental/F6F_G_ratio_age_groups.R", ranges="1-5;10-12;14-15;17-20;22-22;28-31;33-34;37-37;39-43;47-51;54-57;61-69;72-84;87-94;96-96;99-102;104-107", requested_ranges="1-108", role="Original three-group G-ratio plot", code=quote({
library(openxlsx)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rstatix)

input_filepath <- file.path(Sys.getenv("MG_PLOT_INPUT_ROOT", "data/plot_inputs"), "sheet1.xlsx")
plot_name       <- "sheet1"
output_filepath <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "F6F", paste0(plot_name, ".pdf"))

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

out_width  <- 6
out_height <- 2

df <- read.xlsx(input_filepath)

dt <- melt(df,
           id.vars       = NULL,
           measure.vars  = colnames(df),
           variable.name = var_group_name,
           value.name    = var_value_name)

stat.test <- dt %>%
  t_test(Value ~ Group, 
         # comparisons = my_comparisons
  ) %>%
  add_xy_position(x = "Group")

stat.test$p.signif <- symnum(stat.test$p, 
                             corr = FALSE, na = FALSE, 
                             cutpoints = my_symnum$cutpoints, 
                             symbols = my_symnum$symbols)

mytheme2 <- theme_classic() +
  theme(
    text              = element_text(family = "sans", colour = "black", size = 12),
    axis.line         = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks.length = unit(1.5, units = "mm"),
    legend.position   = "none",
    plot.margin       = unit(x = plot_margins, units = "inches")
  )

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

p_final <- p_base + stat_pvalue_manual(
  stat.test,
  label        = "p.signif",
  tip.length   = 0.02,
  bracket.size = 0.5,
  angle        = -90,  # 让星号顺时针旋转90度，与括号平行
  coord.flip   = TRUE  # <--- 关键参数：告诉它画布被翻转了，星号和括号会自动对齐
)

print(p_final)

out_dir <- dirname(output_filepath)
if(!dir.exists(out_dir) && out_dir != ".") {
  dir.create(out_dir, recursive = TRUE)
}

ggsave(filename = output_filepath, 
       plot     = p_final, 
       width    = out_width, 
       height   = out_height)
  }))
))

# F6G  OLM transplantation/remyelination schematic and CC-saline/CC-HMC3/CC-OLM TEM
# Kind: experimental_schematic_and_micrographs
mg_register_panel("F6G", NULL)

# F6H  G-ratios for control, HMC3 and OLM transplantation groups
# Kind: experimental_quantification
mg_register_panel("F6H", list(
  list(source="source/recovered/experimental/F6H_G_ratio_transplant.R", ranges="3-10;12-14;17-19;21-29;32-52;54-54;56-66;68-69;71-108;113-130;132-133;182-185;193-194;198-199", requested_ranges="1-134;182-185;193-194;198-199", role="Original G-ratio branch with Welch tests and BH-adjusted plotted labels", code=quote({
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggbeeswarm)
  library(patchwork)
})

input_file <- file.path(Sys.getenv("MG_PLOT_INPUT_ROOT", "data/plot_inputs"), "2026-7-29-统计图1.xlsx")
output_dir <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "F6H")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw_g <- read_excel(input_file, sheet = "Sheet1", range = "A1:J146", col_names = FALSE)
group_cols <- c(3, 7, 10)
group_names <- as.character(unlist(raw_g[6, group_cols]))

g_ratio <- lapply(seq_along(group_cols), function(i) {
  tibble(
    group = group_names[i],
    value = suppressWarnings(as.numeric(unlist(raw_g[8:nrow(raw_g), group_cols[i]])))
  )
}) |>
  bind_rows() |>
  filter(!is.na(value)) |>
  mutate(group = factor(group, levels = group_names))

pair_index <- combn(group_names, 2, simplify = FALSE)
pairwise_tests <- lapply(pair_index, function(z) {
  test <- t.test(value ~ group, data = filter(g_ratio, group %in% z))
  tibble(group1 = z[1], group2 = z[2], p_value = test$p.value)
}) |>
  bind_rows() |>
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adjusted <= 0.001 ~ "***",
      p_adjusted <= 0.01 ~ "**",
      p_adjusted <= 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    significance_raw = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

print(pairwise_tests |> select(group1, group2, p_value, p_adjusted))

summary_g <- g_ratio |>
  group_by(group) |>
  summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    q1 = quantile(value, 0.25),
    q3 = quantile(value, 0.75),
    .groups = "drop"
  )

palette_g <- setNames(c("#0077c1", "#00a99e", "#6bc72b"), group_names)
set.seed(20260729)

p_g <- ggplot(g_ratio, aes(x = value, y = group, color = group)) +
  geom_boxplot(
    aes(fill = group),
    width = 0.40, outlier.shape = NA,
    linewidth = 0.60, staplewidth = 0
  ) +
  geom_point(
    size = 0.80,
    position = position_jitter(
      width = 0, height = 0.15, seed = 20260729
    )
  ) +
  scale_color_manual(values = palette_g, guide = "none") +
  scale_fill_manual(
    values = scales::alpha(palette_g, 0.30),
    guide = "none"
  ) +
  scale_x_continuous(
    name = "G-ratio",
    limits = c(0.48, 1.025),
    breaks = seq(0.5, 0.9, 0.1),
    expand = expansion(mult = c(0.01, 0))
  ) +
  labs(y = NULL) +
  theme_classic() +
  theme(
    text = element_text(
      family = "sans", colour = "black", size = 12
    ),
    axis.line = element_line(linewidth = 0.60, colour = "black"),
    axis.ticks = element_line(linewidth = 0.60, colour = "black"),
    axis.ticks.length = grid::unit(1.5, units = "mm"),
    legend.position = "none",
    plot.margin = grid::unit(
      c(0.2, 0.2, 0.2, 0.2), units = "inches"
    )
  ) +
  coord_cartesian(clip = "off")

bracket_x <- c(0.972, 0.998, 0.950)
add_comparison_brackets <- function(plot_object, labels) {
  for (i in seq_len(nrow(pairwise_tests))) {
    y1 <- match(pairwise_tests$group1[i], group_names)
    y2 <- match(pairwise_tests$group2[i], group_names)
    x <- bracket_x[i]
    plot_object <- plot_object +
      annotate("segment", x = x, xend = x, y = y1, yend = y2, linewidth = 0.5) +
      annotate("segment", x = x - 0.007, xend = x, y = y1, yend = y1, linewidth = 0.5) +
      annotate("segment", x = x - 0.007, xend = x, y = y2, yend = y2, linewidth = 0.5) +
      annotate(
        "text", x = x + 0.009, y = (y1 + y2) / 2,
        label = labels[i], angle = -90,
        size = 3.5
      )
  }
  plot_object
}

p_g_raw <- add_comparison_brackets(p_g, pairwise_tests$significance_raw)
p_g <- add_comparison_brackets(p_g, pairwise_tests$significance)

ggsave(file.path(output_dir, "图1_G-ratio_散点箱线图.png"), p_g, width = 6, height = 2, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "图1_G-ratio_散点箱线图.pdf"), p_g, width = 6, height = 2, device = cairo_pdf)
ggsave(file.path(output_dir, "图1B_G-ratio_未校正P值_散点箱线图.png"), p_g_raw, width = 6, height = 2, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "图1B_G-ratio_未校正P值_散点箱线图.pdf"), p_g_raw, width = 6, height = 2, device = cairo_pdf)

write.csv(summary_g, file.path(output_dir, "G-ratio_描述统计.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(pairwise_tests, file.path(output_dir, "G-ratio_Welch_t检验_BH校正.csv"), row.names = FALSE, fileEncoding = "UTF-8")

print(summary_g)
print(pairwise_tests)
  }))
))
