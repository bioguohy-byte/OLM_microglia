#!/usr/bin/env Rscript

# Render separate, binary overlays on the SAME saved diagnostic UMAP coordinates.
# Each study receives one Original-MG0 PNG/PDF and one MG00 PNG/PDF.  The two
# labels are deliberately not combined into overlap categories in one plot.

rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
invisible(gc())

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
combine_root <- file.path(project_root, "combine")
old_table_file <- file.path(
  combine_root,
  "results/MG0_MG5_doublet_diagnosis/table/six_study_cell_metadata_with_all_doublet_calls.csv"
)
MG00_table_file <- file.path(
  combine_root,
  "results/full_count_reanalysis/03_current_MG00_doublet_diagnosis/table/six_study_cell_metadata_with_all_three_doublet_calls_MG00.csv"
)
result_root <- file.path(
  combine_root,
  "results/full_count_reanalysis/04_original_MG0_and_MG00_separate_UMAPs"
)
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
runtime_dir <- file.path(result_root, "runtime")
for (directory in c(figure_dir, table_dir, runtime_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}
if (any(!file.exists(c(old_table_file, MG00_table_file)))) {
  stop("One or both required per-cell tables are missing.")
}

study_order <- c(
  "Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
  "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"
)
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width = 7.5, height = 4.5) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
    dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
    device = grDevices::pdf, bg = "white", limitsize = FALSE)
}

old <- read.csv(old_table_file, check.names = FALSE, stringsAsFactors = FALSE)
MG00 <- read.csv(MG00_table_file, check.names = FALSE, stringsAsFactors = FALSE)
required_old <- c("Cell", "Study", "CellType1")
required_MG00 <- c("Cell", "Study", "original_barcode", "MG00_status", "UMAP_1", "UMAP_2")
if (length(setdiff(required_old, names(old)))) {
  stop("Original table lacks: ", paste(setdiff(required_old, names(old)), collapse = ", "))
}
if (length(setdiff(required_MG00, names(MG00)))) {
  stop("MG00 table lacks: ", paste(setdiff(required_MG00, names(MG00)), collapse = ", "))
}

# The same source barcode must be paired inside each study.  Unpaired cells
# remain gray in the Original-MG0 panel: they are not assigned an old label.
old$key <- paste(old$Study, old$Cell, sep = "||")
MG00$key <- paste(MG00$Study, MG00$original_barcode, sep = "||")
if (anyDuplicated(old$key)) stop("Original table has duplicate Study/barcode keys.")
if (anyDuplicated(MG00$key)) stop("MG00 table has duplicate Study/barcode keys.")
old_index <- match(MG00$key, old$key)
if (!any(!is.na(old_index))) stop("No MG00 cells linked to original labels.")
MG00$original_label_available <- !is.na(old_index)
MG00$original_MG0 <- MG00$original_label_available & old$CellType1[old_index] == "MG0"
MG00$is_MG00 <- MG00$MG00_status == "MG00 (manual cluster 2 transfer)"
MG00$original_MG0_overlay <- factor(
  ifelse(MG00$original_MG0, "Original MG0", "Other cells"),
  levels = c("Other cells", "Original MG0")
)
MG00$MG00_overlay <- factor(
  ifelse(MG00$is_MG00, "MG00", "Other cells"),
  levels = c("Other cells", "MG00")
)
if (!setequal(study_order, unique(as.character(MG00$Study)))) {
  stop("Expected exactly the six selected studies.")
}

umap_base <- function(dat, title, subtitle) {
  ggplot(dat, aes(UMAP_1, UMAP_2)) +
    coord_equal() +
    labs(title = title, subtitle = subtitle, x = "Diagnostic UMAP 1", y = "Diagnostic UMAP 2") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      plot.margin = margin(t = 12, r = 14, b = 10, l = 14),
      legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 8)
    )
}
make_original_plot <- function(dat, study) {
  dat <- dat[order(dat$original_MG0_overlay), , drop = FALSE]
  umap_base(dat, "Original MG0", study) +
    geom_point(aes(colour = original_MG0_overlay), size = 0.48, alpha = 0.88) +
    scale_colour_manual(values = c("Other cells" = "#D9D9D9", "Original MG0" = "#2C7FB8")) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
}
make_MG00_plot <- function(dat, study) {
  dat <- dat[order(dat$MG00_overlay), , drop = FALSE]
  umap_base(dat, "MG00", study) +
    geom_point(aes(colour = MG00_overlay), size = 0.48, alpha = 0.88) +
    scale_colour_manual(values = c("Other cells" = "#D9D9D9", MG00 = "#C51B7D")) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
}

plot_store <- list(original_MG0 = list(), MG00 = list())
for (study in study_order) {
  dat <- MG00[MG00$Study == study, , drop = FALSE]
  original_plot <- make_original_plot(dat, study)
  MG00_plot <- make_MG00_plot(dat, study)
  plot_store$original_MG0[[study]] <- original_plot +
    labs(title = study, subtitle = NULL) + theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5), legend.position = "right")
  plot_store$MG00[[study]] <- MG00_plot +
    labs(title = study, subtitle = NULL) + theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5), legend.position = "right")
  save_plot_both(original_plot, file.path(figure_dir, paste0("UMAP_", safe_name(study), "_original_MG0")))
  save_plot_both(MG00_plot, file.path(figure_dir, paste0("UMAP_", safe_name(study), "_MG00")))
}
for (plot_type in c("original_MG0", "MG00")) {
  composite_title <- if (plot_type == "original_MG0") {
    "Six-study diagnostic UMAPs: original MG0"
  } else {
    "Six-study diagnostic UMAPs: MG00"
  }
  composite <- wrap_plots(plot_store[[plot_type]][study_order], ncol = 3) +
    plot_annotation(title = composite_title)
  save_plot_both(
    composite,
    file.path(figure_dir, paste0("Composite_UMAP_", plot_type, "_six_studies")),
    15, 9
  )
  print(composite)
}

mapping_table <- MG00[, c(
  "Cell", "Study", "original_barcode", "original_label_available", "original_MG0",
  "MG00_status", "is_MG00", "UMAP_1", "UMAP_2"
)]
write.csv(mapping_table, file.path(table_dir, "per_cell_original_MG0_and_MG00_mapping.csv"), row.names = FALSE)
summary_table <- aggregate(
  Cell ~ Study + original_label_available + original_MG0 + is_MG00,
  data = mapping_table,
  FUN = length
)
names(summary_table)[5] <- "n_cells"
write.csv(summary_table, file.path(table_dir, "original_MG0_and_MG00_counts_by_study.csv"), row.names = FALSE)

writeLines(c(
  "Each study has two separate binary overlays on the same MG00 diagnostic UMAP coordinates.",
  "Original-MG0 overlay: blue marks only cells whose linked historical CellType1 label is MG0; all remaining cells, including unavailable old labels, are gray.",
  "MG00 overlay: magenta marks only the renamed manual integrated-cluster-2 transfer; all other cells are gray.",
  paste0("Linked original labels: ", sum(MG00$original_label_available), " / ", nrow(MG00), " MG00 diagnostic cells."),
  "No overlap-category plot is generated and no object, count, QC status, UMAP coordinate, or doublet call is changed."
), file.path(runtime_dir, "separate_overlay_method_manifest.txt"))
writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "separate_overlay_R_sessionInfo.txt"))
message("Separate Original-MG0 and MG00 UMAPs are ready: ", figure_dir)
