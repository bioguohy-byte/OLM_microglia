#!/usr/bin/env Rscript
# Full-slice four-colour MERFISH plotting on the original, unrotated
# center coordinates. Source RDS objects and their metadata are read only.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
data_dir <- file.path(project_dir, "data/processed_data")
out_dir <- file.path(project_dir, "results/figure/code_verification_full_slice_four_colour_center_coordinates_1to1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_colours <- c(
  "Other cells" = "grey80",
  "OLG" = "blue",
  "OPC" = "yellow",
  "ST_MG0" = "red"
)
plot_levels <- names(plot_colours)

for (sample_id in c("s2r1", "s2r2", "s2r3")) {
  obj <- readRDS(file.path(data_dir, paste0("merfish_", sample_id, "_v2.rds")))
  stopifnot(all(c("CellType", "center_x", "center_y") %in% colnames(obj@meta.data)))

  # The original full-slice display labels horizontal axis as y and vertical axis
  # as x. No rotation, mirroring, or use of the Seurat spatial reduction occurs.
  df <- obj@meta.data %>%
    transmute(
      center_x = as.numeric(center_x),
      center_y = as.numeric(center_y),
      four_colour_celltype = if_else(
        CellType %in% c("OLG", "OPC", "ST_MG0"), as.character(CellType), "Other cells"
      )
    ) %>%
    mutate(four_colour_celltype = factor(four_colour_celltype, levels = plot_levels))

  p <- ggplot(df, aes(x = center_y, y = center_x, colour = four_colour_celltype)) +
    geom_point(size = 0.25) +
    scale_colour_manual(values = plot_colours, drop = FALSE) +
    coord_fixed(ratio = 1, expand = FALSE) +
    labs(
      title = paste0("Whole-slice raw-coordinate plot: ", sample_id),
      subtitle = "Original center_y (horizontal) and center_x (vertical); 1:1 coordinate scale",
      x = "y (center_y)", y = "x (center_x)", colour = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )

  stem <- file.path(out_dir, paste0(sample_id, "_whole_slice_OLG_OPC_ST_MG0_other_center_coordinates_1to1"))
  ggsave(paste0(stem, ".pdf"), p, width = 8, height = 8, bg = "white")
  ggsave(paste0(stem, ".png"), p, width = 8, height = 8, dpi = 300, bg = "white")
  rm(obj, df, p)
  gc()
}
