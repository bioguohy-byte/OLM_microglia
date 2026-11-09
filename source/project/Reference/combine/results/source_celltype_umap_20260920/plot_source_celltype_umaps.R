#!/usr/bin/env Rscript

# Render source-object native cell-type metadata on each source object's stored UMAP.
# This is a read-only visualization: source Seurat objects are never saved or relabeled.

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: plot_source_celltype_umaps.R <manifest.tsv> <output_dir>")
}

manifest <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "png"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# All source-object annotation columns that represent existing cell-type labels.
# Where two such columns exist, both are plotted rather than selecting or reconciling one.
source_label_columns <- list(
  Chia = c("cell.type"),
  Fullard = c("celltype"),
  Grubman = c("cellType", "CellType"),
  Jakel = c("Celltypes", "cell_type"),
  Su = c("CellType")
)

label_table <- function(values, study, column) {
  values <- as.character(values)
  values[is.na(values) | values == ""] <- "<missing>"
  tab <- sort(table(values), decreasing = TRUE)
  data.frame(
    study = study,
    source_metadata_column = column,
    native_celltype = names(tab),
    n_cells = as.integer(tab),
    stringsAsFactors = FALSE
  )
}

make_source_umap <- function(obj, reduction, source_column, study) {
  values <- obj@meta.data[[source_column]]
  values <- as.character(values)
  values[is.na(values) | values == ""] <- "<missing>"
  ordered_levels <- names(sort(table(values), decreasing = TRUE))
  obj@meta.data$.source_celltype_for_plot <- factor(values, levels = ordered_levels)
  n_labels <- length(ordered_levels)
  # Explicit high-chroma colours keep the source labels readable in rasterized PDF pages.
  source_palette <- grDevices::hcl(
    h = seq(15, 375, length.out = n_labels + 1L)[-1L],
    c = 95,
    l = 48,
    fixup = TRUE
  )

  DimPlot(
    obj,
    reduction = reduction,
    group.by = ".source_celltype_for_plot",
    cols = source_palette,
    alpha = 1,
    raster = TRUE,
    raster.dpi = c(512, 512),
    pt.size = 0.30,
    label = FALSE
  ) +
    coord_fixed() +
    ggtitle(paste0(study, ": source `", source_column, "`")) +
    labs(subtitle = paste0(n_labels, " native labels; stored UMAP; labels shown as-is"), colour = "Source cell type") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = if (n_labels > 20) 6.5 else 8),
      legend.key.height = grid::unit(if (n_labels > 20) 0.26 else 0.35, "cm")
    )
}

inventory <- list()
page_plots <- list()
for (i in seq_len(nrow(manifest))) {
  study <- manifest$study[[i]]
  object_path <- manifest$object_path[[i]]
  reduction <- manifest$reduction[[i]]
  message("Rendering ", study)
  obj <- qread(object_path)
  if (!reduction %in% Reductions(obj)) stop(study, ": stored reduction not found: ", reduction)
  columns <- source_label_columns[[study]]
  if (is.null(columns)) stop(study, ": no source annotation columns declared")
  missing_columns <- setdiff(columns, colnames(obj@meta.data))
  if (length(missing_columns)) stop(study, ": source label column(s) absent: ", paste(missing_columns, collapse = ", "))

  panels <- lapply(columns, function(column) {
    inventory[[paste(study, column, sep = "__")]] <<- label_table(obj@meta.data[[column]], study, column)
    make_source_umap(obj, reduction, column, study)
  })
  page <- wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title = paste0(study, ": source-object native CellType metadata"),
      subtitle = "Context only. No project-defined CellType, reviewer suggestion, or source-object metadata was changed."
    ) &
    theme(plot.title = element_text(face = "bold", size = 18), plot.subtitle = element_text(size = 10))
  page_plots[[study]] <- page

  ggsave(
    filename = file.path(output_dir, "png", paste0(study, "_source_native_celltype_UMAP.png")),
    plot = page,
    width = if (length(columns) == 1L) 15 else 28,
    height = 11,
    dpi = 180,
    bg = "white",
    limitsize = FALSE
  )
  rm(obj)
  gc()
}

pdf_path <- file.path(output_dir, "Source_object_native_celltype_UMAPs.pdf")
grDevices::pdf(pdf_path, width = 28, height = 13, onefile = TRUE, useDingbats = FALSE)
for (study in names(page_plots)) print(page_plots[[study]])
grDevices::dev.off()

inventory_out <- do.call(rbind, inventory)
write.csv(inventory_out, file.path(output_dir, "tables", "source_native_celltype_inventory.csv"), row.names = FALSE)
writeLines(c(
  "# Source-object native CellType UMAPs",
  "",
  "Each panel uses a pre-existing cell-type metadata column and the UMAP reduction already stored in that source object.",
  "The values are rendered as-is. No new CellType label is calculated, no reviewer suggestion is used, and no source Seurat object is saved or modified.",
  "Grubman and Jakel each contain two existing annotation columns; both are shown independently rather than being harmonized."
), file.path(output_dir, "README.md"))
