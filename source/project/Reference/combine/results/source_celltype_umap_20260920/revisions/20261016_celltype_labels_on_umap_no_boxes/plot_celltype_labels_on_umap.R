#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(scattermore)
  library(scales)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: plot_celltype_labels_on_umap.R <task_dir>")
out_dir <- normalizePath(args[[1]], mustWork = TRUE)
manifest <- read.delim(file.path(out_dir, "plot_manifest.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
for (d in c("pdf", "png", "supplementary/pdf", "supplementary/png", "tables", "audit")) {
  dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
}
overrides_path <- file.path(out_dir, "label_position_overrides.tsv")
overrides <- if (file.exists(overrides_path)) read.delim(overrides_path, stringsAsFactors = FALSE) else data.frame()
pdf_save <- function(p, path, width, height) {
  ggsave(path, plot = p, device = grDevices::cairo_pdf,
         width = width, height = height, units = "in", bg = "white", limitsize = FALSE)
}
png_save <- function(p, path, width, height) {
  ggsave(path, plot = p, device = ragg::agg_png,
         width = width, height = height, units = "in", dpi = 240, bg = "white", limitsize = FALSE)
}
plots <- list()
inventories <- list()
plot_audit <- list()
label_audit <- list()
dimensions <- list()
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  id <- row$plot_id
  study <- row$study
  message(format(Sys.time(), "%F %T"), " Plotting ", id)
  bundle <- readRDS(file.path(out_dir, "plot_data", row$bundle_file))
  stopifnot(identical(rownames(bundle$umap), bundle$cells),
            identical(rownames(bundle$metadata), bundle$cells),
            row$annotation_column %in% colnames(bundle$metadata),
            all(is.finite(bundle$umap)))
  values <- as.character(bundle$metadata[[row$annotation_column]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  counts <- sort(table(values), decreasing = TRUE)
  levels <- names(counts)
  n_types <- length(counts)
  palette <- setNames(grDevices::hcl(h = seq(15, 375, length.out = n_types + 1L)[-1L],
                                    c = 95, l = 48, fixup = TRUE), levels)
  coords <- bundle$umap
  xr <- range(coords[, 1])
  yr <- range(coords[, 2])
  plot_data <- data.frame(x = coords[, 1], y = coords[, 2],
                          celltype = factor(values, levels = levels), row.names = bundle$cells)
  anchors <- lapply(levels, function(label) {
    idx <- which(values == label)
    med <- apply(coords[idx, , drop = FALSE], 2, stats::median)
    dist <- ((coords[idx, 1] - med[1]) / max(diff(xr), .Machine$double.eps)) ^ 2 +
      ((coords[idx, 2] - med[2]) / max(diff(yr), .Machine$double.eps)) ^ 2
    anchor <- idx[[which.min(dist)]]
    display <- paste(strwrap(label, width = 25), collapse = "\n")
    data.frame(celltype = label, n_cells = unname(counts[[label]]),
      x = coords[anchor, 1], y = coords[anchor, 2], anchor_cell = bundle$cells[[anchor]],
      label_text = paste0(display, "\n(", comma(unname(counts[[label]])), ")"),
      nudge_x = 0, nudge_y = 0, colour = palette[[label]], stringsAsFactors = FALSE)
  })
  anchors <- do.call(rbind, anchors)
  if (nrow(overrides)) {
    for (j in which(overrides$plot_id == id)) {
      k <- match(overrides$celltype[[j]], anchors$celltype)
      if (is.na(k)) stop("Unknown override label for ", id, ": ", overrides$celltype[[j]])
      anchors$nudge_x[[k]] <- overrides$nudge_x[[j]]
      anchors$nudge_y[[k]] <- overrides$nudge_y[[j]]
    }
  }
  width <- 12.5
  height <- 10
  label_size <- if (n_types > 20) 3.2 else 3.7
  p <- ggplot(plot_data, aes(x = x, y = y, colour = celltype)) +
    scattermore::geom_scattermore(pointsize = 2.2, pixels = c(1400, 1400), show.legend = FALSE) +
    scale_colour_manual(values = palette, drop = FALSE) +
    ggrepel::geom_text_repel(data = anchors, aes(x = x, y = y, label = label_text),
      inherit.aes = FALSE, colour = "black",
      size = label_size, lineheight = 0.95, seed = 20261016,
      nudge_x = anchors$nudge_x, nudge_y = anchors$nudge_y,
      box.padding = grid::unit(0.45, "lines"), point.padding = grid::unit(0.15, "lines"),
      min.segment.length = 0, segment.colour = "#454545", segment.size = 0.25,
      force = 2, force_pull = 0.3, max.overlaps = Inf, max.time = 5, max.iter = 100000,
      show.legend = FALSE) +
    coord_fixed(ratio = 1, xlim = xr, ylim = yr, clip = "off") +
    labs(title = paste0(study, " | ", row$title_note),
      subtitle = paste0("Total: ", comma(length(values)), " cells; ", n_types, " cell types; saved UMAP"),
      x = "UMAP_1", y = "UMAP_2") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 19),
      plot.subtitle = element_text(size = 11, colour = "#333333"),
      axis.text = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11),
      legend.position = "none", plot.margin = margin(12, 25, 12, 15))
  folder <- if (row$display_group == "main") "." else "supplementary"
  pdf_name <- file.path(folder, "pdf", paste0(id, "_celltypes_on_UMAP.pdf"))
  png_name <- file.path(folder, "png", paste0(id, "_celltypes_on_UMAP.png"))
  pdf_save(p, file.path(out_dir, pdf_name), width, height)
  png_save(p, file.path(out_dir, png_name), width, height)
  if (interactive()) print(p)
  plots[[id]] <- p
  inventories[[id]] <- data.frame(plot_id = id, study = study, display_group = row$display_group,
    annotation_column = row$annotation_column, celltype = levels,
    n_cells = as.integer(counts), colour = unname(palette), stringsAsFactors = FALSE)
  anchors$plot_id <- id
  anchors$annotation_column <- row$annotation_column
  label_audit[[id]] <- anchors
  dimensions[[id]] <- data.frame(plot_id = id, pdf = pdf_name, png = png_name,
    width_in = width, height_in = height, png_dpi = 240, n_cells = length(values), n_labels = n_types)
  plot_audit[[id]] <- data.frame(plot_id = id, study = study, annotation_column = row$annotation_column,
    n_cells = length(values), n_labels = n_types, x_min = xr[1], x_max = xr[2], y_min = yr[1], y_max = yr[2],
    coordinate_ratio = 1, coordinates_identical_to_bundle = TRUE, labels_identical_to_bundle = TRUE,
    umap_recomputed = FALSE, labels_reassigned = FALSE, legend_shown = FALSE,
    pointsize = 2.2, raster_x = 1400, raster_y = 1400, stringsAsFactors = FALSE)
}
book <- function(group, filename) {
  ids <- manifest$plot_id[manifest$display_group == group]
  if (!length(ids)) return(invisible(NULL))
  grDevices::cairo_pdf(file.path(out_dir, filename), width = 12.5, height = 10, onefile = TRUE)
  for (id in ids) print(plots[[id]])
  grDevices::dev.off()
}
book("main", "Six_study_celltypes_labeled_on_UMAPs.pdf")
book("supplementary", "Supplementary_existing_annotation_columns_on_UMAPs.pdf")
write.csv(do.call(rbind, inventories), file.path(out_dir, "tables/celltype_counts_and_colours.csv"), row.names = FALSE)
write.csv(do.call(rbind, label_audit), file.path(out_dir, "tables/label_anchors_and_display.csv"), row.names = FALSE)
write.csv(do.call(rbind, dimensions), file.path(out_dir, "audit/plot_dimensions.csv"), row.names = FALSE)
write.table(do.call(rbind, plot_audit), file.path(out_dir, "audit/plot_preservation.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
saveRDS(plots, file.path(out_dir, "plot_data/rendered_plot_objects.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "audit/plotting_sessionInfo.txt"))
message("LABEL_PLOTTING_COMPLETE")
