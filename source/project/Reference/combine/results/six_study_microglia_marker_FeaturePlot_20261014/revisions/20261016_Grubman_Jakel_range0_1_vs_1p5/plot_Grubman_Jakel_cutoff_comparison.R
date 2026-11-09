#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(scattermore)
  library(scales)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: plot_Grubman_Jakel_cutoff_comparison.R <output_dir>")
out_dir <- args[[1]]
studies <- c("Grubman", "Jakel")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
upper_limits <- c(range0_1p5 = 1.5, range0_1 = 1)
grey <- "#D3D3D3"
red <- "#E31A1C"
for (d in c("audit", "comparison")) dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
pdf_save <- function(plot, path, width, height) {
  ggsave(path, plot = plot, device = grDevices::cairo_pdf,
         width = width, height = height, units = "in", bg = "white", limitsize = FALSE)
}
png_save <- function(plot, path, width, height, dpi = 240) {
  ggsave(path, plot = plot, device = ragg::agg_png,
         width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
}
scale_audit <- list()
coordinate_audit <- list()
for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Plotting ", study)
  bundle <- readRDS(file.path(out_dir, "plot_data", paste0(study, "_stored_umap_three_RNA_markers.rds")))
  stopifnot(identical(rownames(bundle$umap), bundle$cells),
            identical(colnames(bundle$expression), bundle$cells),
            identical(colnames(bundle$counts), bundle$cells))
  metadata <- data.frame(CellType = as.character(bundle$metadata$CellType), row.names = bundle$cells)
  thin <- CreateSeuratObject(counts = bundle$counts, meta.data = metadata,
                            min.cells = 0, min.features = 0, project = study)
  thin <- SetAssayData(thin, assay = "RNA", layer = "data", new.data = bundle$expression)
  coordinates <- bundle$umap
  colnames(coordinates) <- c("UMAP_1", "UMAP_2")
  thin[["umap"]] <- CreateDimReducObject(embeddings = coordinates, key = "UMAP_", assay = "RNA")
  stopifnot(identical(unname(Embeddings(thin, "umap")), unname(bundle$umap)),
            identical(as.matrix(GetAssayData(thin, assay = "RNA", layer = "data")), as.matrix(bundle$expression)))
  xr <- range(coordinates[, 1])
  yr <- range(coordinates[, 2])
  shared_theme <- theme_classic(base_size = 11) + theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9, colour = "#333333"),
    axis.text = element_text(size = 9, colour = "black"),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.height = grid::unit(0.38, "cm"),
    plot.margin = margin(5, 8, 5, 5)
  )
  comparison_plots <- list()
  for (version in names(upper_limits)) {
    upper <- upper_limits[[version]]
    breaks <- seq(0, upper, by = 0.5)
    current <- list()
    for (d in c("pdf_by_study", "pdf_single_marker", "png_by_study", "png_single_marker")) {
      dir.create(file.path(out_dir, version, d), recursive = TRUE, showWarnings = FALSE)
    }
    for (gene in markers) {
      feature <- bundle$features[[gene]]
      vals <- as.numeric(bundle$expression[feature, ])
      p <- FeaturePlot(thin, features = feature, reduction = "umap", dims = c(1, 2),
        cols = c(grey, red), order = TRUE, pt.size = 1.3,
        min.cutoff = 0, max.cutoff = upper, slot = "data",
        raster = TRUE, raster.dpi = c(1200, 1200), combine = FALSE)[[1]] +
        scale_colour_gradient(low = grey, high = red, limits = c(0, upper), oob = scales::squish,
          breaks = breaks,
          guide = guide_colourbar(title = "Log-normalized\nRNA", barheight = grid::unit(2.6, "cm"),
            barwidth = grid::unit(0.30, "cm"))) +
        coord_fixed(ratio = 1, xlim = xr, ylim = yr) +
        labs(title = paste(study, "|", gene), subtitle = paste0(comma(length(vals)), " cells"),
             x = "UMAP_1", y = "UMAP_2") + shared_theme
      current[[gene]] <- p
      comparison_plots[[paste(version, gene, sep = "__")]] <- p
      pdf_save(p, file.path(out_dir, version, "pdf_single_marker", paste0(study, "_", gene, "_FeaturePlot.pdf")), 4.8, 4.6)
      png_save(p, file.path(out_dir, version, "png_single_marker", paste0(study, "_", gene, "_FeaturePlot.png")), 4.8, 4.6, dpi = 300)
      if (interactive()) print(p)
      scale_audit[[paste(study, version, gene, sep = "__")]] <- data.frame(
        study = study, marker = gene, version = version, lower_limit = 0, upper_limit = upper,
        legend_breaks = paste(breaks, collapse = ";"), colour_low = grey, colour_high = red,
        n_cells = length(vals), n_positive = sum(vals > 0), n_saturated = sum(vals > upper),
        percent_saturated_of_positive = 100 * sum(vals > upper) / max(1, sum(vals > 0)),
        pt_size = 1.3, raster_resolution_x = 1200, raster_resolution_y = 1200, coordinate_ratio = 1)
    }
    triptych <- wrap_plots(current, nrow = 1, guides = "keep") +
      plot_annotation(title = paste0(study, ": microglial marker expression"),
        subtitle = "Stored full-cell UMAP; RNA log1p(count / total RNA count x 10,000); grey to red",
        theme = theme(plot.title = element_text(face = "bold", size = 17), plot.subtitle = element_text(size = 10)))
    pdf_save(triptych, file.path(out_dir, version, "pdf_by_study", paste0(study, "_three_microglia_markers_FeaturePlots.pdf")), 14.4, 5.0)
    png_save(triptych, file.path(out_dir, version, "png_by_study", paste0(study, "_three_microglia_markers_FeaturePlots.png")), 14.4, 5.0)
  }
  comparison <- wrap_plots(comparison_plots, ncol = 3, byrow = TRUE, guides = "keep") +
    plot_annotation(title = paste0(study, ": FeaturePlot colour-range comparison"),
      subtitle = "Top row: 0-1.5; bottom row: 0-1. Same expression values, stored UMAP, cells and plotting parameters.",
      theme = theme(plot.title = element_text(face = "bold", size = 17), plot.subtitle = element_text(size = 10)))
  pdf_save(comparison, file.path(out_dir, "comparison", paste0(study, "_range0_1p5_vs_1_FeaturePlots.pdf")), 14.4, 9.6)
  png_save(comparison, file.path(out_dir, "comparison", paste0(study, "_range0_1p5_vs_1_FeaturePlots.png")), 14.4, 9.6)
  coordinate_audit[[study]] <- data.frame(study = study, n_cells = ncol(thin),
    coordinate_values_identical = TRUE, barcode_order_identical = TRUE,
    expression_values_identical_to_plot_bundle = TRUE, coordinate_ratio = 1, umap_recomputed = FALSE)
  rm(thin, bundle)
  gc(verbose = FALSE)
}
write.table(do.call(rbind, scale_audit), file.path(out_dir, "audit/plot_scale_parameters.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
write.table(do.call(rbind, coordinate_audit), file.path(out_dir, "audit/plot_coordinate_validation.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "audit/plotting_sessionInfo.txt"))
message("CUTOFF_COMPARISON_COMPLETE")
