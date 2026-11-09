#!/usr/bin/env Rscript
# Supplementary Figure 3C: OLM proportion.
# Plot parameters are retained from plot_统计图.R, lines 124-130, 143-169, 175-176.
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[[1L]] else file.path(Sys.getenv("MG_PLOT_INPUT_ROOT", "data/plot_inputs"), "2026-7-29-统计图.xlsx")
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "S3C")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
raw_v <- read_excel(input_file, sheet = "Sheet2", range = "A1:F11", col_names = FALSE)
value <- suppressWarnings(as.numeric(unlist(raw_v[, 6])))
olm <- data.frame(group = "OLM", value = value[!is.na(value)])
stopifnot(nrow(olm) == 5L)
p_olm <- ggplot(olm, aes(x = group, y = value)) +
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
    name = expression(Trim2^"+" * Iba1^"+" * OLM/Iba1^"+" ~ MG),
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
ggsave(file.path(output_dir, "S3C_OLM_proportion.png"), p_olm, width = 3.0, height = 4.3, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "S3C_OLM_proportion.pdf"), p_olm, width = 3.0, height = 4.3, device = cairo_pdf)
write.csv(olm, file.path(output_dir, "S3C_plot_values.csv"), row.names = FALSE)
cat("S3C values:\n")
print(olm, digits = 17)
cat("\nSession:\n")
print(sessionInfo())
