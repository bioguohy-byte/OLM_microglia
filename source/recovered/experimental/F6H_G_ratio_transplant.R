#!/usr/bin/env Rscript

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

# Sheet1: group names are in row 6; observations begin in row 8.
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

# Pairwise Welch t tests; the plotted labels use BH-adjusted P values.
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

# Add three neatly nested pairwise brackets at the right of the data.
# Adjacent-group comparisons are inside; the comparison spanning all
# three rows is placed outside, matching the supplied reference style.
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

# Sheet2: five observations are stored in column F.
raw_v <- read_excel(input_file, sheet = "Sheet2", range = "A1:F11", col_names = FALSE)
viral <- tibble(
  group = "HMC3-GFP",
  value = suppressWarnings(as.numeric(unlist(raw_v[, 6])))
) |>
  filter(!is.na(value))

summary_v <- viral |>
  summarise(
    group = first(group),
    n = n(),
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    q1 = quantile(value, 0.25),
    q3 = quantile(value, 0.75)
  )

p_v <- ggplot(viral, aes(x = group, y = value)) +
  geom_boxplot(
    width = 0.38, outlier.shape = NA,
    fill = "#1786C0", alpha = 0.35,
    color = "#1786C0", linewidth = 0.9,
    staplewidth = 0.55
  ) +
  geom_point(
    shape = 21, size = 2.6, stroke = 0.9,
    fill = "#1786C0", color = "black",
    position = position_jitter(width = 0.075, height = 0, seed = 20260729)
  ) +
  scale_y_continuous(
    name = "Viral transduction",
    breaks = seq(0, 0.12, 0.02),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL) +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(linewidth = 0.9),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black"),
    plot.margin = margin(8, 8, 7, 7)
  ) +
  coord_cartesian(ylim = c(0, 0.12))

ggsave(file.path(output_dir, "图1_G-ratio_散点箱线图.png"), p_g, width = 6, height = 2, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "图1_G-ratio_散点箱线图.pdf"), p_g, width = 6, height = 2, device = cairo_pdf)
ggsave(file.path(output_dir, "图1B_G-ratio_未校正P值_散点箱线图.png"), p_g_raw, width = 6, height = 2, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "图1B_G-ratio_未校正P值_散点箱线图.pdf"), p_g_raw, width = 6, height = 2, device = cairo_pdf)
ggsave(file.path(output_dir, "图2_病毒转导率_箱线图.png"), p_v, width = 3.0, height = 4.3, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "图2_病毒转导率_箱线图.pdf"), p_v, width = 3.0, height = 4.3, device = cairo_pdf)

combined <- p_g + p_v + plot_layout(widths = c(2.3, 1))
ggsave(file.path(output_dir, "图1图2_组合版.png"), combined, width = 9.2, height = 4.0, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "图1图2_组合版.pdf"), combined, width = 9.2, height = 4.0, device = cairo_pdf)

write.csv(summary_g, file.path(output_dir, "G-ratio_描述统计.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(pairwise_tests, file.path(output_dir, "G-ratio_Welch_t检验_BH校正.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(summary_v, file.path(output_dir, "病毒转导率_描述统计.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("Finished. Outputs:", output_dir, "\n")
print(summary_g)
print(pairwise_tests)
print(summary_v)
