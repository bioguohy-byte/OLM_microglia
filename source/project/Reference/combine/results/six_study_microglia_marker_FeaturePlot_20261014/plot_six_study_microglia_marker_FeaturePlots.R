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
if (length(args) != 1L) stop("Usage: plot_six_study_microglia_marker_FeaturePlots.R <output_dir>")
out_dir <- args[[1]]
for (d in c("pdf_by_study", "pdf_single_marker", "png_by_study", "png_single_marker", "audit")) {
  dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
}
studies <- c("Su", "Chia", "Fullard", "Jakel", "Grubman", "Yao")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
grey <- "#D3D3D3"
red <- "#E31A1C"
ranges <- read.delim(file.path(out_dir, "audit/marker_expression_ranges.tsv"), check.names = FALSE)
stopifnot(nrow(ranges) == 18L, setequal(unique(ranges$study), studies))
normalization_audit <- read.delim(file.path(out_dir, "audit/source_normalization_commands.tsv"), check.names = FALSE)
normalization_audit$parameters[normalization_audit$command == "<not recorded>"] <-
  "Source NormalizeData command not recorded; plotting marker values derived from RNA/counts"
write.table(normalization_audit, file.path(out_dir, "audit/source_normalization_commands.tsv"),
            sep = "\t", row.names = FALSE, quote = TRUE)
# Each gene has one common colour scale across all six studies. The upper
# bound is the largest study-specific 99.5th percentile among positive cells.
# At most 0.5% of expressing cells in each panel exceed this upper bound.
upper_limits <- setNames(vapply(markers, function(gene) {
  max(ranges$positive_q995[ranges$marker == gene])
}, numeric(1)), markers)
if (any(!is.finite(upper_limits)) || any(upper_limits <= 0)) stop("Invalid colour upper limits")

pdf_save <- function(plot, path, width, height) {
  ggsave(path, plot = plot, device = grDevices::cairo_pdf,
         width = width, height = height, units = "in", bg = "white", limitsize = FALSE)
}
png_save <- function(plot, path, width, height, dpi = 240) {
  ggsave(path, plot = plot, device = ragg::agg_png,
         width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
}
reference_plots <- list()
marker_plots <- list()
page_plots <- list()
scale_audit <- list()
coordinate_audit <- list()

for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Plotting ", study)
  bundle <- readRDS(file.path(out_dir, "plot_data", paste0(study, "_stored_umap_three_RNA_markers.rds")))
  stopifnot(identical(rownames(bundle$umap), bundle$cells),
            identical(colnames(bundle$expression), bundle$cells),
            identical(colnames(bundle$counts), bundle$cells))
  metadata <- data.frame(CellType = as.character(bundle$metadata$CellType), row.names = bundle$cells)
  # A small plotting-only object carries the exact three-row RNA matrix and
  # the saved full-cell embedding. No PCA, UMAP, QC, or clustering is run.
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
  mg <- metadata$CellType %in% bundle$mg_labels
  reference_data <- data.frame(x = coordinates[, 1], y = coordinates[, 2],
    colour = ifelse(mg, red, grey))
  reference_data <- reference_data[order(mg), , drop = FALSE]
  reference <- ggplot(reference_data, aes(x = x, y = y, colour = colour)) +
    scattermore::geom_scattermore(pointsize = 1.3, pixels = c(1200, 1200)) +
    scale_colour_identity() +
    coord_fixed(ratio = 1, xlim = xr, ylim = yr) +
    labs(title = paste0(study, " | existing MG label"),
      subtitle = paste0("CellType: ", paste(bundle$mg_labels, collapse = "; "), "; n = ", comma(sum(mg))),
      x = "UMAP_1", y = "UMAP_2") + shared_theme
  reference_plots[[study]] <- reference
  current <- list()
  for (gene in markers) {
    feature <- bundle$features[[gene]]
    upper <- upper_limits[[gene]]
    vals <- as.numeric(bundle$expression[feature, ])
    p <- FeaturePlot(thin, features = feature, reduction = "umap", dims = c(1, 2),
      cols = c(grey, red), order = TRUE, pt.size = 1.3,
      min.cutoff = 0, max.cutoff = upper, slot = "data",
      raster = TRUE, raster.dpi = c(1200, 1200), combine = FALSE)[[1]] +
      scale_colour_gradient(low = grey, high = red, limits = c(0, upper), oob = scales::squish,
        breaks = pretty(c(0, upper), n = 3),
        guide = guide_colourbar(title = "Log-normalized\nRNA", barheight = grid::unit(2.6, "cm"),
          barwidth = grid::unit(0.30, "cm"))) +
      coord_fixed(ratio = 1, xlim = xr, ylim = yr) +
      labs(title = paste(study, "|", gene), subtitle = paste0(comma(length(vals)), " cells"),
           x = "UMAP_1", y = "UMAP_2") + shared_theme
    current[[gene]] <- p
    marker_plots[[paste(study, gene, sep = "__")]] <- p
    pdf_save(p, file.path(out_dir, "pdf_single_marker", paste0(study, "_", gene, "_FeaturePlot.pdf")), 4.8, 4.6)
    png_save(p, file.path(out_dir, "png_single_marker", paste0(study, "_", gene, "_FeaturePlot.png")), 4.8, 4.6, dpi = 300)
    scale_audit[[paste(study, gene, sep = "__")]] <- data.frame(
      study = study, marker = gene, colour_low = grey, colour_high = red,
      lower_limit = 0, upper_limit = upper, scale_shared_across_studies = TRUE,
      n_cells = length(vals), n_positive = sum(vals > 0), n_saturated = sum(vals > upper),
      percent_saturated_of_positive = 100 * sum(vals > upper) / max(1, sum(vals > 0)),
      pt_size = 1.3, raster_resolution_x = 1200, raster_resolution_y = 1200,
      coordinate_ratio = 1, stringsAsFactors = FALSE
    )
  }
  triptych <- wrap_plots(current, nrow = 1, guides = "keep") +
    plot_annotation(title = paste0(study, ": microglial marker expression"),
      subtitle = "Stored full-cell UMAP; RNA log1p(count / total RNA count x 10,000); grey to red",
      theme = theme(plot.title = element_text(face = "bold", size = 17), plot.subtitle = element_text(size = 10)))
  pdf_save(triptych, file.path(out_dir, "pdf_by_study", paste0(study, "_three_microglia_markers_FeaturePlots.pdf")), 14.4, 5.0)
  png_save(triptych, file.path(out_dir, "png_by_study", paste0(study, "_three_microglia_markers_FeaturePlots.png")), 14.4, 5.0)
  page_plots[[study]] <- wrap_plots(c(list(reference), current), nrow = 1) +
    plot_annotation(title = paste0(study, ": existing microglial annotation and marker evidence"),
      subtitle = paste0("Same ", comma(ncol(thin)), " cells and stored UMAP in all panels; coordinate ratio 1:1"),
      theme = theme(plot.title = element_text(face = "bold", size = 17), plot.subtitle = element_text(size = 10)))
  coordinate_audit[[study]] <- data.frame(study = study, n_cells = ncol(thin),
    coordinate_values_identical = TRUE, barcode_order_identical = TRUE,
    expression_values_identical_to_plot_bundle = TRUE, coordinate_ratio = 1,
    umap_recomputed = FALSE, stringsAsFactors = FALSE)
  rm(thin, bundle)
  gc(verbose = FALSE)
}

main <- wrap_plots(marker_plots, ncol = 3, byrow = TRUE, guides = "keep") +
  plot_annotation(title = "Microglial marker expression across six human datasets",
    subtitle = "P2RY12, CX3CR1 and CSF1R on stored full-cell UMAPs; RNA log-normalized expression; grey to red",
    caption = "Each marker uses one colour scale across studies. Upper bound: maximum study-specific positive-cell 99.5th percentile. Coordinate ratio 1:1.",
    theme = theme(plot.title = element_text(face = "bold", size = 19),
      plot.subtitle = element_text(size = 11), plot.caption = element_text(size = 9, hjust = 0)))
pdf_save(main, file.path(out_dir, "Six_study_P2RY12_CX3CR1_CSF1R_FeaturePlots.pdf"), 14.4, 25.0)
png_save(main, file.path(out_dir, "Six_study_P2RY12_CX3CR1_CSF1R_FeaturePlots.png"), 14.4, 25.0, dpi = 220)
book_path <- file.path(out_dir, "Six_study_existing_MG_annotation_and_marker_audit.pdf")
grDevices::cairo_pdf(book_path, width = 19.2, height = 5.3, onefile = TRUE)
for (study in studies) print(page_plots[[study]])
grDevices::dev.off()
write.table(do.call(rbind, scale_audit), file.path(out_dir, "audit/plot_scale_parameters.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
write.table(do.call(rbind, coordinate_audit), file.path(out_dir, "audit/plot_coordinate_validation.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "audit/plotting_sessionInfo.txt"))
writeLines(c(
  "# Six-study microglial marker FeaturePlots",
  "",
  "The main figure contains 18 panels: six existing full-cell UMAPs x P2RY12/CX3CR1/CSF1R.",
  "All cells are plotted. UMAP coordinates, cell barcodes, and existing CellType labels are preserved.",
  "Su uses the saved 2026-09-18 full-cell Harmony audit checkpoint because its base source object has no stored UMAP.",
  "The other five studies use their original full-cell v2 objects and stored umap reduction.",
  "Only three marker rows are normalized for plotting as log1p(count / full RNA library size x 10000).",
  "Some stored RNA/data layers contain raw counts; the original layers are preserved in the plotting bundles.",
  "For each marker, one colour scale is shared across all six studies; its upper limit is the maximum positive-cell q99.5 across studies.",
  "The six-page audit book adds an existing MG label reference panel, with exactly the same cells, axes, and coordinate ratio.",
  "Jakel's reference highlights only the existing Microglia_Macrophages label. Historical inclusion of Immune_cells/ImOlGs is not automatically endorsed.",
  "Marker expression is review evidence; no final celltype decision or metadata modification is made here.",
  "pdf_by_study/: six individual triptychs. pdf_single_marker/: 18 standalone marker panels.",
  "The corresponding PNGs, exact plotting data, input inventory, expression-by-label audit, colour limits, and R session information are also saved."
), file.path(out_dir, "README.md"))
message("PLOTTING_COMPLETE")
