# Whole-slice four-colour spatial plots from the retained MERFISH slice objects.
# No object or metadata file is modified.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
data_dir <- file.path(project_dir, "data/processed_data")
out_dir <- file.path(project_dir, "results/figure/code_verification_full_slice_four_colour")
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
  stopifnot("CellType" %in% colnames(obj@meta.data), "spatial" %in% Reductions(obj))

  obj$four_colour_celltype <- ifelse(
    obj$CellType %in% c("OLG", "OPC", "ST_MG0"),
    as.character(obj$CellType),
    "Other cells"
  )
  obj$four_colour_celltype <- factor(obj$four_colour_celltype, levels = plot_levels)

  p <- DimPlot(
    obj,
    reduction = "spatial",
    group.by = "four_colour_celltype",
    cols = plot_colours,
    pt.size = 0.25,
    shuffle = FALSE
  ) +
    ggtitle(paste0("Whole-slice spatial plot: ", sample_id)) +
    theme(plot.title = element_text(hjust = 0.5))

  stem <- file.path(out_dir, paste0(sample_id, "_whole_slice_OLG_OPC_ST_MG0_other_verified"))
  ggsave(paste0(stem, ".pdf"), p, width = 7, height = 7, bg = "white")
  ggsave(paste0(stem, ".png"), p, width = 7, height = 7, dpi = 300, bg = "white")
}
