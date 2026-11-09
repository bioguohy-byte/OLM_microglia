#!/usr/bin/env Rscript

# Corrected display-only companion to stage 2.
#
# Stage 2's artificial-profile objects, parent-cell pairing and all comparison
# statistics are unchanged.  This script recreates the two human-definition
# UMAPs per study, assigning the approved parent barcodes through a named
# metadata data.frame and then asserting that every selected barcode is shown.
# The figures are saved under new `03_`/`04_` filenames so the original stage-2
# files remain intact as a record of the earlier display run.

local({
  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(qs)
    library(ggplot2)
  })

  project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
  combine_root <- file.path(project_root, "combine")
  per_study_dir <- file.path(combine_root, "results/full_count_reanalysis/01_per_study/object")
  stage2_root <- file.path(combine_root, "results/OL_MG_artificial_reference/02_artificial_profile_construction")
  figure_dir <- file.path(stage2_root, "figure")
  table_dir <- file.path(stage2_root, "table")
  runtime_dir <- file.path(stage2_root, "runtime")
  ol_selection_file <- file.path(combine_root, "results/OL_MG_artificial_reference/01_OL_parent_review/manual_selections/OL_parent_cluster_selection_APPROVED.csv")
  mg_selection_file <- file.path(combine_root, "results/full_count_reanalysis/manual_selections/per_study_Microglia_cluster_selection_TEMPLATE.csv")
  for (dir_path in c(figure_dir, table_dir, runtime_dir)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  save_plot_both <- function(plot, stem, width = 8, height = 6.5) {
    ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 320, bg = "white")
    ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = grDevices::pdf, bg = "white")
  }
  make_checked_binary_plot <- function(object, positive_cells, status_col, positive_label, title, colour) {
    positive_cells <- intersect(as.character(positive_cells), colnames(object))
    if (!length(positive_cells)) stop(title, ": no selected parent cells are present in this object.")
    status <- rep("Other current DoubletFinder singlet", ncol(object))
    names(status) <- colnames(object)
    status[positive_cells] <- positive_label
    object[[status_col]] <- data.frame(
      value = factor(status, levels = c("Other current DoubletFinder singlet", positive_label)),
      row.names = colnames(object)
    )
    assigned <- as.character(object[[status_col]][, 1])
    if (sum(assigned == positive_label, na.rm = TRUE) != length(positive_cells)) {
      stop(title, ": metadata assignment failed; expected ", length(positive_cells), " selected cells but displayed ", sum(assigned == positive_label, na.rm = TRUE), ".")
    }
    # Do not use Seurat::DimPlot here.  On this server its rendering path could
    # flatten the selected category to grey although ggplot_build() contained
    # the correct colours.  A direct ggplot2 scatter plot keeps the exact UMAP
    # coordinates and writes the colours explicitly to the graphics device.
    coordinates <- as.data.frame(Seurat::Embeddings(object, reduction = "umap"))
    if (ncol(coordinates) < 2L) stop(title, ": UMAP has fewer than two dimensions.")
    coordinates$UMAP_1_display <- coordinates[[1L]]
    coordinates$UMAP_2_display <- coordinates[[2L]]
    coordinates$plot_label <- factor(status[rownames(coordinates)], levels = c("Other current DoubletFinder singlet", positive_label))
    if (anyNA(coordinates$plot_label)) stop(title, ": one or more UMAP cells are missing a plotting label.")
    plot <- ggplot(coordinates, aes(x = UMAP_1_display, y = UMAP_2_display)) +
      geom_point(
        data = coordinates[coordinates$plot_label == "Other current DoubletFinder singlet", , drop = FALSE],
        colour = "#D9D9D9", size = 0.18, alpha = 0.80
      ) +
      geom_point(
        data = coordinates[coordinates$plot_label == positive_label, , drop = FALSE],
        colour = colour, size = 0.20, alpha = 1
      ) +
      ggtitle(title) +
      labs(x = "UMAP 1", y = "UMAP 2") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    list(plot = plot, n_positive = length(positive_cells))
  }

  required_files <- c(ol_selection_file, mg_selection_file)
  if (any(!file.exists(required_files))) stop("Missing required input(s): ", paste(required_files[!file.exists(required_files)], collapse = "; "))
  ol_selection <- read.csv(ol_selection_file, check.names = FALSE, stringsAsFactors = FALSE)
  mg_selection <- read.csv(mg_selection_file, check.names = FALSE, stringsAsFactors = FALSE)
  study_order <- c("Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021", "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024")
  target_study <- Sys.getenv("OLMG_STUDY", unset = "")
  if (nzchar(target_study)) {
    if (!target_study %in% study_order) stop("OLMG_STUDY is not one of the six configured studies: ", target_study)
    study_order <- target_study
  }
  audit_rows <- list()
  figure_rows <- list()

  for (study in study_order) {
    study_safe <- safe_name(study)
    message("\n===== Corrected OL/MG definition UMAP: ", study, " =====")
    object_file <- file.path(per_study_dir, paste0(study_safe, "_full_count_processed.qs"))
    if (!file.exists(object_file)) stop(study, ": missing full-count processed object.")
    object <- qs::qread(object_file)
    if (!"umap" %in% names(object@reductions)) stop(study, ": missing UMAP.")
    metadata <- object[[]]
    if (!"analysis_cluster" %in% names(metadata)) stop(study, ": missing analysis_cluster metadata.")
    ol_clusters <- as.character(ol_selection$analysis_cluster[
      ol_selection$Study == study & ol_selection$select_as_mature_OL %in% TRUE
    ])
    mg_clusters <- as.character(mg_selection$analysis_cluster[
      mg_selection$Study == study & mg_selection$select_as_Microglia %in% TRUE
    ])
    ol_cells <- rownames(metadata)[as.character(metadata$analysis_cluster) %in% ol_clusters]
    mg_cells <- rownames(metadata)[as.character(metadata$analysis_cluster) %in% mg_clusters]
    if (length(intersect(ol_cells, mg_cells))) stop(study, ": approved OL and canonical MG parent sets overlap.")

    ol_result <- make_checked_binary_plot(object, ol_cells, "corrected_mature_OL_parent", "Confirmed mature OL", paste0(study, ": reviewer-confirmed mature OL"), "#B2182B")
    mg_result <- make_checked_binary_plot(object, mg_cells, "corrected_canonical_MG_parent", "Confirmed canonical MG", paste0(study, ": prior canonical Microglia definition"), "#2166AC")
    ol_stem <- file.path(figure_dir, paste0(study_safe, "_03_confirmed_mature_OL_UMAP_corrected"))
    mg_stem <- file.path(figure_dir, paste0(study_safe, "_04_confirmed_canonical_MG_UMAP_corrected"))
    save_plot_both(ol_result$plot, ol_stem)
    save_plot_both(mg_result$plot, mg_stem)
    print(ol_result$plot)
    print(mg_result$plot)
    audit_rows[[study]] <- data.frame(
      Study = study,
      approved_OL_clusters = paste(ol_clusters, collapse = ";"),
      displayed_OL_parent_singlets = ol_result$n_positive,
      approved_canonical_MG_clusters = paste(mg_clusters, collapse = ";"),
      displayed_canonical_MG_parent_singlets = mg_result$n_positive,
      stringsAsFactors = FALSE
    )
    figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = study, figure = "confirmed_mature_OL_UMAP_corrected", stem = ol_stem, stringsAsFactors = FALSE)
    figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = study, figure = "confirmed_canonical_MG_UMAP_corrected", stem = mg_stem, stringsAsFactors = FALSE)
    write.csv(audit_rows[[study]], file.path(table_dir, paste0(study_safe, "_corrected_OL_MG_parent_definition_UMAP_audit.csv")), row.names = FALSE)
    rm(object, metadata)
    gc(verbose = FALSE)
  }

  audit_files <- list.files(table_dir, pattern = "_corrected_OL_MG_parent_definition_UMAP_audit\\.csv$", full.names = TRUE)
  audit <- do.call(rbind, lapply(audit_files, read.csv, check.names = FALSE, stringsAsFactors = FALSE))
  corrected_pngs <- list.files(figure_dir, pattern = "_confirmed_(mature_OL|canonical_MG)_UMAP_corrected\\.png$", full.names = TRUE)
  figure_manifest <- data.frame(
    Study = sub("_(03|04)_confirmed_.*$", "", basename(corrected_pngs)),
    figure = sub("^[^_]+.*_(03|04)_", "", sub("\\.png$", "", basename(corrected_pngs))),
    stem = sub("\\.png$", "", corrected_pngs),
    stringsAsFactors = FALSE
  )
  figure_manifest$png_exists <- file.exists(paste0(figure_manifest$stem, ".png"))
  figure_manifest$pdf_exists <- file.exists(paste0(figure_manifest$stem, ".pdf"))
  if (!all(figure_manifest$png_exists & figure_manifest$pdf_exists)) stop("At least one corrected definition UMAP lacks a PNG or PDF counterpart.")
  write.csv(audit, file.path(table_dir, "corrected_OL_MG_parent_definition_UMAP_audit.csv"), row.names = FALSE)
  write.csv(figure_manifest, file.path(table_dir, "corrected_OL_MG_parent_definition_UMAP_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "24b_parent_definition_UMAP_refresh_sessionInfo.txt"))
  message("\nCorrected parent-definition UMAPs complete. Every approved OL and canonical MG barcode was assigned and counted before plotting.")
})
