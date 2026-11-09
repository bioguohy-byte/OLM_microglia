#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

output_dir <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "S4D")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 输入的百分数换算为小数：8.244%、10.91%、8.24%。
olm_data <- data.frame(
  group = "OLM",
  value = c(0.08244, 0.10910, 0.08240),
  point_x = c(0.94, 1.00, 1.06)
)

plot_olm <- ggplot(olm_data, aes(x = 1, y = value)) +
  geom_boxplot(
    aes(group = 1),
    width = 0.38,
    outlier.shape = NA,
    fill = "#1786C0",
    alpha = 0.35,
    color = "#1786C0",
    linewidth = 0.9,
    staplewidth = 0.55
  ) +
  geom_point(
    aes(x = point_x),
    shape = 21,
    size = 2.6,
    stroke = 0.9,
    fill = "#1786C0",
    color = "black"
  ) +
  scale_x_continuous(
    breaks = 1,
    labels = "OLM",
    limits = c(0.5, 1.5),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    name = expression(Trim2^"+" * Iba1^"+" * OLM/Iba1^"+" ~ MG),
    breaks = seq(0, 0.12, 0.02),
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL) +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(linewidth = 0.9),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    plot.margin = margin(8, 8, 7, 7)
  ) +
  coord_cartesian(ylim = c(0, 0.12))

print(plot_olm)

ggsave(
  file.path(output_dir, "OLM_三值箱线图.png"),
  plot_olm,
  width = 3.0,
  height = 4.3,
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(output_dir, "OLM_三值箱线图.pdf"),
  plot_olm,
  width = 3.0,
  height = 4.3,
  device = cairo_pdf
)

cat("输入值（%）：", paste(olm_data$value * 100, collapse = ", "), "\n")
cat("输出目录：", output_dir, "\n", sep = "")
