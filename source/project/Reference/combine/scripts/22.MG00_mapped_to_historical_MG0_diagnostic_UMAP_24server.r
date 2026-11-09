#!/usr/bin/env Rscript

# Project the reanalysis-derived MG00 label and its three reanalysis doublet
# calls onto the *saved historical MG0--MG5 diagnostic UMAP coordinates*.
#
# This script performs no reprocessing and does not alter a Seurat object,
# doublet call, barcode, or label.  It joins records only by Study plus the
# source barcode, audits every match, then writes individual and composite
# figures as PNG and PDF.  The original historical coordinates are retained
# exactly as they occur in the historical per-cell results table.

rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
invisible(gc())

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
combine_root <- file.path(project_root, "combine")

historical_file <- file.path(
  combine_root,
  "results/MG0_MG5_doublet_diagnosis/table/six_study_cell_metadata_with_all_doublet_calls.csv"
)
reanalysis_file <- file.path(
  combine_root,
  "results/full_count_reanalysis/03_current_MG00_doublet_diagnosis/table/six_study_cell_metadata_with_all_three_doublet_calls_MG00.csv"
)
result_root <- file.path(
  combine_root,
  "results/full_count_reanalysis/05_MG00_mapped_to_historical_MG0_diagnostic_UMAPs"
)
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
runtime_dir <- file.path(result_root, "runtime")
for (directory in c(figure_dir, table_dir, runtime_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

if (any(!file.exists(c(historical_file, reanalysis_file)))) {
  stop("Required historical or reanalysis per-cell table is missing.")
}

study_order <- c(
  "Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
  "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"
)
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width = 7.5, height = 4.5) {
  ggsave(
    paste0(stem, ".png"), plot, width = width, height = height,
    dpi = 320, bg = "white", limitsize = FALSE
  )
  ggsave(
    paste0(stem, ".pdf"), plot, width = width, height = height,
    device = grDevices::pdf, bg = "white", limitsize = FALSE
  )
}

historical <- read.csv(historical_file, check.names = FALSE, stringsAsFactors = FALSE)
reanalysis <- read.csv(reanalysis_file, check.names = FALSE, stringsAsFactors = FALSE)
required_historical <- c("Cell", "Study", "CellType1", "UMAP_1", "UMAP_2")
required_reanalysis <- c(
  "Cell", "Study", "original_barcode", "MG00_status",
  "DoubletFinder_call", "scDblFinder_call", "Scrublet_call"
)
if (length(setdiff(required_historical, names(historical)))) {
  stop("Historical table lacks: ", paste(setdiff(required_historical, names(historical)), collapse = ", "))
}
if (length(setdiff(required_reanalysis, names(reanalysis)))) {
  stop("Reanalysis table lacks: ", paste(setdiff(required_reanalysis, names(reanalysis)), collapse = ", "))
}
if (!setequal(study_order, unique(historical$Study))) stop("Historical table does not contain exactly six expected studies.")
if (!setequal(study_order, unique(reanalysis$Study))) stop("Reanalysis table does not contain exactly six expected studies.")

# `historical$Cell` is the source barcode.  The reanalysis records preserve the
# same source barcode in `original_barcode`; Study prevents collisions across
# studies.  The historical table remains the plotting base, so its UMAP never
# changes and cells absent from the reanalysis remain visible in gray.
historical$key <- paste(historical$Study, historical$Cell, sep = "||")
reanalysis$key <- paste(reanalysis$Study, reanalysis$original_barcode, sep = "||")
if (anyDuplicated(historical$key)) stop("Historical Study/barcode keys are duplicated.")
if (anyDuplicated(reanalysis$key)) stop("Reanalysis Study/original_barcode keys are duplicated.")

idx <- match(historical$key, reanalysis$key)
historical$reanalysis_record_available <- !is.na(idx)
historical$MG00_mapped <- historical$reanalysis_record_available &
  reanalysis$MG00_status[idx] == "MG00 (manual cluster 2 transfer)"
historical$reanalysis_DoubletFinder_call <- ifelse(
  historical$reanalysis_record_available, reanalysis$DoubletFinder_call[idx], NA_character_
)
historical$reanalysis_scDblFinder_call <- ifelse(
  historical$reanalysis_record_available, reanalysis$scDblFinder_call[idx], NA_character_
)
historical$reanalysis_Scrublet_call <- ifelse(
  historical$reanalysis_record_available, reanalysis$Scrublet_call[idx], NA_character_
)
historical$historical_MG0 <- historical$CellType1 == "MG0"

binary_overlay <- function(dat, positive, positive_label, title, subtitle, positive_colour) {
  # Cells unavailable in the reanalysis have no caller result.  They are
  # historical background cells for visualization, not a third "NA" class.
  positive <- !is.na(positive) & positive
  dat$overlay <- factor(
    ifelse(positive, positive_label, "Other historical cell"),
    levels = c("Other historical cell", positive_label)
  )
  dat <- dat[order(dat$overlay), , drop = FALSE]
  ggplot(dat, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(colour = overlay), size = 0.48, alpha = 0.88) +
    # Build the named palette dynamically.  Using `positive_label =` here
    # would create a literal palette name rather than the caller-supplied
    # label, causing positive cells to fall back to ggplot's default gray.
    scale_colour_manual(
      values = setNames(
        c("#D9D9D9", positive_colour),
        c("Other historical cell", positive_label)
      )
    ) +
    coord_equal() +
    labs(
      title = title, subtitle = subtitle,
      x = "Historical diagnostic UMAP 1", y = "Historical diagnostic UMAP 2"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      plot.margin = margin(t = 12, r = 14, b = 10, l = 14),
      legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 8)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
}

plot_specs <- list(
  historical_MG0 = list(
    selector = function(dat) dat$historical_MG0,
    label = "Historical MG0", title = "Historical MG0", colour = "#2C7FB8"
  ),
  MG00_mapped = list(
    selector = function(dat) dat$MG00_mapped,
    label = "MG00 mapped from reanalysis", title = "MG00 mapped to historical diagnostic UMAP", colour = "#C51B7D"
  ),
  DoubletFinder_reanalysis_call_mapped = list(
    selector = function(dat) dat$reanalysis_DoubletFinder_call == "Doublet",
    label = "DoubletFinder-called doublet in reanalysis",
    title = "Reanalysis DoubletFinder call mapped to historical diagnostic UMAP", colour = "#D7301F"
  ),
  scDblFinder_reanalysis_call_mapped = list(
    selector = function(dat) dat$reanalysis_scDblFinder_call == "Doublet",
    label = "scDblFinder-called doublet in reanalysis",
    title = "Reanalysis scDblFinder call mapped to historical diagnostic UMAP", colour = "#E6550D"
  ),
  Scrublet_reanalysis_call_mapped = list(
    selector = function(dat) dat$reanalysis_Scrublet_call == "Doublet",
    label = "Scrublet-called doublet in reanalysis",
    title = "Reanalysis Scrublet call mapped to historical diagnostic UMAP", colour = "#756BB1"
  )
)

plot_store <- lapply(names(plot_specs), function(x) list())
names(plot_store) <- names(plot_specs)
for (study in study_order) {
  dat <- historical[historical$Study == study, , drop = FALSE]
  for (plot_name in names(plot_specs)) {
    spec <- plot_specs[[plot_name]]
    plot <- binary_overlay(dat, spec$selector(dat), spec$label, spec$title, study, spec$colour)
    plot_store[[plot_name]][[study]] <- plot +
      labs(title = study, subtitle = NULL) +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5), legend.position = "right")
    save_plot_both(plot, file.path(figure_dir, paste0("UMAP_", safe_name(study), "_", plot_name)))
  }
}

for (plot_name in names(plot_specs)) {
  composite <- wrap_plots(plot_store[[plot_name]][study_order], ncol = 3) +
    plot_annotation(title = paste0("Six-study historical diagnostic UMAPs: ", plot_specs[[plot_name]]$title))
  save_plot_both(
    composite,
    file.path(figure_dir, paste0("Composite_UMAP_", plot_name, "_six_studies")),
    width = 15, height = 9
  )
  print(composite)
}

per_cell <- historical[, c(
  "Cell", "Study", "Sample", "CellType1", "UMAP_1", "UMAP_2",
  "reanalysis_record_available", "historical_MG0", "MG00_mapped",
  "reanalysis_DoubletFinder_call", "reanalysis_scDblFinder_call", "reanalysis_Scrublet_call"
)]
write.csv(per_cell, file.path(table_dir, "historical_coordinate_projection_per_cell.csv"), row.names = FALSE)

summary_by_study <- do.call(rbind, lapply(study_order, function(study) {
  dat <- historical[historical$Study == study, , drop = FALSE]
  data.frame(
    Study = study,
    historical_cells = nrow(dat),
    matched_reanalysis_cells = sum(dat$reanalysis_record_available),
    unmatched_historical_cells = sum(!dat$reanalysis_record_available),
    mapping_rate = mean(dat$reanalysis_record_available),
    historical_MG0_cells = sum(dat$historical_MG0),
    MG00_mapped_cells = sum(dat$MG00_mapped),
    DoubletFinder_reanalysis_doublets_mapped = sum(dat$reanalysis_DoubletFinder_call == "Doublet", na.rm = TRUE),
    scDblFinder_reanalysis_doublets_mapped = sum(dat$reanalysis_scDblFinder_call == "Doublet", na.rm = TRUE),
    Scrublet_reanalysis_doublets_mapped = sum(dat$reanalysis_Scrublet_call == "Doublet", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.csv(summary_by_study, file.path(table_dir, "historical_coordinate_projection_summary_by_study.csv"), row.names = FALSE)

figure_files <- list.files(figure_dir, pattern = "\\.(png|pdf)$", full.names = FALSE)
png_stems <- sub("\\.png$", "", figure_files[grepl("\\.png$", figure_files)])
pdf_stems <- sub("\\.pdf$", "", figure_files[grepl("\\.pdf$", figure_files)])
pair_check <- data.frame(
  figure_stem = sort(unique(c(png_stems, pdf_stems))),
  PNG = sort(unique(c(png_stems, pdf_stems))) %in% png_stems,
  PDF = sort(unique(c(png_stems, pdf_stems))) %in% pdf_stems,
  stringsAsFactors = FALSE
)
write.csv(pair_check, file.path(table_dir, "historical_coordinate_projection_PNG_PDF_pair_check.csv"), row.names = FALSE)

writeLines(c(
  "Purpose: project MG00 and reanalysis doublet calls onto historical MG0--MG5 diagnostic UMAP coordinates.",
  "Coordinates were read exclusively from the historical MG0--MG5 per-cell table; they were not recomputed.",
  "MG00 and all three reanalysis caller labels were matched only by Study plus original source barcode.",
  paste0("Historical cells: ", nrow(historical)),
  paste0("Historical cells with a reanalysis record: ", sum(historical$reanalysis_record_available)),
  paste0("Reanalysis records available for matching: ", nrow(reanalysis)),
  "Gray points include historical cells without a matching reanalysis record, so the historical diagnostic UMAP remains complete.",
  "No Seurat object, raw count matrix, cell barcode, doublet score, doublet call or MG00 label was changed."
), file.path(runtime_dir, "method_manifest.txt"))
writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "R_sessionInfo.txt"))

message(
  "Historical-coordinate projection complete. Matched ",
  sum(historical$reanalysis_record_available), " / ", nrow(historical),
  " historical cells. Output: ", figure_dir
)
