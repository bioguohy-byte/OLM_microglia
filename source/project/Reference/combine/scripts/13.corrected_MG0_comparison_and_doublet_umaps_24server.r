#!/usr/bin/env Rscript

# Presentation figure correction: compare the original MG0--MG5 MG0 label with
# the user-selected reanalysis MG0 (integrated cluster 2) on the same diagnostic
# UMAP.  For the same reanalysis UMAP, also render three independently computed
# doublet calls.  This script only reads saved per-cell result tables and writes
# new figures/tables; it does not alter any analytical Seurat object.
#
# Every figure is saved as PNG and PDF.  Four six-study composite views are
# printed to the RStudio Plots history after all individual figures are saved.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
combine_root <- file.path(project_root, "combine")
old_table <- file.path(
  combine_root,
  "results/MG0_MG5_doublet_diagnosis/table/six_study_cell_metadata_with_all_doublet_calls.csv"
)
reanalysis_table <- file.path(
  combine_root,
  "results/full_count_reanalysis/03_current_MG0_doublet_diagnosis/table/six_study_cell_metadata_with_all_three_doublet_calls.csv"
)
result_root <- file.path(
  combine_root,
  "results/full_count_reanalysis/04_presentation_MG0_comparison_and_doublet_UMAPs"
)
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
runtime_dir <- file.path(result_root, "runtime")
for (directory in c(figure_dir, table_dir, runtime_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

study_order <- c(
  "Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
  "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"
)

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width = 7.5, height = 4.5) {
  ggsave(
    filename = paste0(stem, ".png"), plot = plot, width = width, height = height,
    dpi = 320, bg = "white", limitsize = FALSE
  )
  ggsave(
    filename = paste0(stem, ".pdf"), plot = plot, width = width, height = height,
    device = grDevices::pdf, bg = "white", limitsize = FALSE
  )
}

if (!file.exists(old_table)) stop("Original MG0--MG5 per-cell table is missing: ", old_table)
if (!file.exists(reanalysis_table)) stop("Reanalysis per-cell table is missing: ", reanalysis_table)

old <- read.csv(old_table, check.names = FALSE, stringsAsFactors = FALSE)
current <- read.csv(reanalysis_table, check.names = FALSE, stringsAsFactors = FALSE)

required_old <- c("Cell", "Study", "CellType1")
required_current <- c(
  "Cell", "Study", "original_barcode", "MG0_status", "UMAP_1", "UMAP_2",
  "DoubletFinder_call", "scDblFinder_call", "Scrublet_call"
)
if (length(setdiff(required_old, names(old)))) {
  stop("Original table lacks: ", paste(setdiff(required_old, names(old)), collapse = ", "))
}
if (length(setdiff(required_current, names(current)))) {
  stop("Reanalysis table lacks: ", paste(setdiff(required_current, names(current)), collapse = ", "))
}

# `Cell` in the original table is the raw source barcode.  Reanalysis prefixes
# its active `Cell` name, while preserving the source barcode as original_barcode.
# Include Study in the key because raw barcodes can repeat across studies.
old$key <- paste(old$Study, old$Cell, sep = "||")
current$key <- paste(current$Study, current$original_barcode, sep = "||")
if (anyDuplicated(old$key)) stop("The original table has duplicated Study/barcode keys.")
if (anyDuplicated(current$key)) stop("The reanalysis table has duplicated Study/original_barcode keys.")

old_index <- match(current$key, old$key)
if (!any(!is.na(old_index))) stop("No reanalysis cells matched original Study/barcode labels.")

current$old_label_available <- !is.na(old_index)
current$original_MG0 <- current$old_label_available & old$CellType1[old_index] == "MG0"
current$reanalysis_MG0 <- current$MG0_status == "MG0 (manual cluster 2 transfer)"
current$MG0_definition_status <- ifelse(
  !current$old_label_available & current$reanalysis_MG0, "Old label unavailable; reanalysis MG0",
  ifelse(
    !current$old_label_available, "Old label unavailable",
    ifelse(
      current$original_MG0 & current$reanalysis_MG0, "Both definitions",
      ifelse(
        current$original_MG0, "Original MG0 only",
        ifelse(current$reanalysis_MG0, "Reanalysis MG0 only", "Neither definition")
      )
    )
  )
)
current$MG0_definition_status <- factor(
  current$MG0_definition_status,
  levels = c(
    "Old label unavailable", "Neither definition", "Original MG0 only",
    "Reanalysis MG0 only", "Old label unavailable; reanalysis MG0", "Both definitions"
  )
)
current$DoubletFinder_call <- factor(current$DoubletFinder_call, levels = c("Singlet", "Doublet"))
current$scDblFinder_call <- factor(current$scDblFinder_call, levels = c("Singlet", "Doublet"))
current$Scrublet_call <- factor(current$Scrublet_call, levels = c("Singlet", "Doublet"))

if (!setequal(study_order, unique(current$Study))) {
  stop("Reanalysis table does not contain exactly the expected six studies.")
}

umap_base <- function(dat, title, subtitle = NULL) {
  ggplot(dat, aes(UMAP_1, UMAP_2)) +
    coord_equal() +
    labs(title = title, subtitle = subtitle, x = "Diagnostic UMAP 1", y = "Diagnostic UMAP 2") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12, lineheight = 0.95),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      plot.margin = margin(t = 12, r = 14, b = 10, l = 14),
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.position = "bottom"
    )
}

make_definition_plot <- function(dat, study) {
  dat <- dat[order(dat$MG0_definition_status), , drop = FALSE]
  umap_base(dat, "Original and reanalysis MG0", subtitle = study) +
    geom_point(aes(colour = MG0_definition_status), size = 0.48, alpha = 0.88) +
    scale_colour_manual(values = c(
      "Old label unavailable" = "#FEE08B",
      "Neither definition" = "#D9D9D9",
      "Original MG0 only" = "#2C7FB8",
      "Reanalysis MG0 only" = "#C51B7D",
      "Old label unavailable; reanalysis MG0" = "#E7298A",
      "Both definitions" = "#756BB1"
    ), labels = c(
      "Old label unavailable", "Neither MG0", "Original MG0 only",
      "Reanalysis MG0 only", "Old label unavailable + reanalysis MG0", "Both definitions"
    ), drop = FALSE) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.2, alpha = 1)))
}

make_doublet_plot <- function(dat, column, method, study) {
  dat <- dat[order(dat[[column]]), , drop = FALSE]
  umap_base(dat, method, subtitle = study) +
    geom_point(aes(colour = .data[[column]]), size = 0.48, alpha = 0.88) +
    scale_colour_manual(values = c(Singlet = "#BDBDBD", Doublet = "#D7301F")) +
    guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1)))
}

plot_types <- c("MG0_definition_comparison", "DoubletFinder", "scDblFinder", "Scrublet")
plot_store <- setNames(lapply(plot_types, function(x) list()), plot_types)
for (study in study_order) {
  dat <- current[current$Study == study, , drop = FALSE]
  plots <- list(
    MG0_definition_comparison = make_definition_plot(dat, study),
    DoubletFinder = make_doublet_plot(dat, "DoubletFinder_call", "DoubletFinder", study),
    scDblFinder = make_doublet_plot(dat, "scDblFinder_call", "scDblFinder", study),
    Scrublet = make_doublet_plot(dat, "Scrublet_call", "Scrublet", study)
  )
  for (plot_type in names(plots)) {
    plot_store[[plot_type]][[study]] <- plots[[plot_type]] +
      labs(title = study) +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
    save_plot_both(
      plots[[plot_type]],
      file.path(figure_dir, paste0("UMAP_", safe_name(study), "_", plot_type))
    )
  }
}

for (plot_type in plot_types) {
  composite_title <- switch(
    plot_type,
    MG0_definition_comparison = "Six-study diagnostic UMAPs: original MG0 label and reanalysis MG0",
    paste0("Six-study diagnostic UMAPs: ", plot_type)
  )
  composite <- wrap_plots(plot_store[[plot_type]][study_order], ncol = 3) +
    plot_annotation(title = composite_title)
  save_plot_both(
    composite,
    file.path(figure_dir, paste0("Composite_UMAP_", plot_type, "_six_studies")),
    width = 15, height = 9
  )
  print(composite)
}

mapping_table <- current[, c(
  "Cell", "Study", "original_barcode", "original_MG0", "reanalysis_MG0",
  "MG0_definition_status", "DoubletFinder_call", "scDblFinder_call", "Scrublet_call",
  "UMAP_1", "UMAP_2"
)]
write.csv(
  mapping_table,
  file.path(table_dir, "per_cell_original_vs_reanalysis_MG0_mapping.csv"),
  row.names = FALSE
)
summary_table <- aggregate(
  Cell ~ Study + MG0_definition_status,
  data = mapping_table,
  FUN = length
)
names(summary_table)[3] <- "n_cells"
write.csv(
  summary_table,
  file.path(table_dir, "original_vs_reanalysis_MG0_counts_by_study.csv"),
  row.names = FALSE
)

n_matched <- sum(!is.na(old_index))
writeLines(c(
  "Purpose: presentation correction using only previously saved per-cell tables.",
  paste0("Reanalysis QC-passed cells: ", nrow(current)),
  paste0("Cells linked to the original MG0--MG5 table by Study plus original barcode: ", n_matched),
  "Original MG0 refers to CellType1 == MG0 in the original MG0--MG5-label pipeline.",
  "Reanalysis MG0 refers to the user-selected integrated cluster 2 transferred back to the pre-removal diagnostic UMAP.",
  "Doublet plots use the exact calls already calculated by DoubletFinder, scDblFinder and Scrublet.",
  "No Seurat object, barcode, doublet call or label is modified by this script."
), file.path(runtime_dir, "method_manifest.txt"))

message(
  "Corrected figures complete. Matched ", n_matched, " / ", nrow(current),
  " reanalysis cells to the original MG0--MG5 labels. Output: ", figure_dir
)
