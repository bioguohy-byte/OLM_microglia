#!/usr/bin/env Rscript

# OL–MG artificial-doublet reference, stage 1
# Review mature oligodendrocyte (OL) parent clusters in the six studies.
#
# This script is intentionally READ-ONLY with respect to existing Seurat objects
# and existing microglia selections.  It creates only review figures, tables,
# and an empty user-editable OL-selection template.  No OL cluster is selected
# automatically: the selection is a biological decision reserved for review.
#
# Input: six post-DoubletFinder singlet per-study Seurat objects, retaining raw
#        integer RNA counts and the existing analysis_cluster / batch metadata.
# Output: results/OL_MG_artificial_reference/01_OL_parent_review/

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(qs)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260909)

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
combine_root <- file.path(project_root, "combine")
per_study_object_dir <- file.path(combine_root, "results/full_count_reanalysis/01_per_study/object")
existing_mg_selection <- file.path(
  combine_root,
  "results/full_count_reanalysis/manual_selections/per_study_Microglia_cluster_selection_TEMPLATE.csv"
)
result_root <- file.path(combine_root, "results/OL_MG_artificial_reference/01_OL_parent_review")
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
runtime_dir <- file.path(result_root, "runtime")
manual_dir <- file.path(result_root, "manual_selections")
for (d in c(figure_dir, table_dir, runtime_dir, manual_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

study_order <- c(
  "Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
  "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"
)

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

save_plot_both <- function(plot, stem, width, height) {
  png_file <- paste0(stem, ".png")
  pdf_file <- paste0(stem, ".pdf")
  ggsave(png_file, plot = plot, width = width, height = height, dpi = 320, bg = "white")
  ggsave(pdf_file, plot = plot, width = width, height = height, device = grDevices::pdf, bg = "white")
  if (!file.exists(png_file) || !file.exists(pdf_file)) {
    stop("PNG/PDF pair was not written for: ", stem)
  }
  invisible(c(png = png_file, pdf = pdf_file))
}

get_assay_layer <- function(object, assay = "RNA", layer = "counts") {
  out <- tryCatch(
    SeuratObject::LayerData(object, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (is.null(out)) {
    out <- tryCatch(
      SeuratObject::GetAssayData(object, assay = assay, layer = layer),
      error = function(e) NULL
    )
  }
  if (is.null(out) && identical(layer, "counts")) {
    out <- tryCatch(
      SeuratObject::GetAssayData(object, assay = assay, slot = "counts"),
      error = function(e) NULL
    )
  }
  if (is.null(out) && identical(layer, "data")) {
    out <- tryCatch(
      SeuratObject::GetAssayData(object, assay = assay, slot = "data"),
      error = function(e) NULL
    )
  }
  if (is.null(out)) stop("Cannot retrieve RNA ", layer, " layer.")
  if (!inherits(out, "sparseMatrix")) out <- Matrix::Matrix(out, sparse = TRUE)
  out
}

cluster_order <- function(x) {
  ux <- unique(as.character(x))
  as.character(ux[order(suppressWarnings(as.numeric(ux)), ux, na.last = TRUE)])
}

assert_integer_counts <- function(x, label) {
  if (length(x@x) && any(abs(x@x - round(x@x)) > 1e-8)) {
    stop(label, " RNA counts are not integer values.")
  }
  invisible(TRUE)
}

# These two panels are independent of the prior MG00 definition.  In particular,
# TRIM2/RNF220/ST18/OLIG2/MBP/MOG are deliberately excluded from the primary
# OL-versus-MG comparison, so that parent selection does not become circular.
ol_markers <- c("PLP1", "MOBP", "MAG", "CNP", "CLDN11", "SOX10", "MYRF", "OPALIN")
mg_markers <- c(
  "P2RY12", "TMEM119", "CX3CR1", "GPR34", "CSF1R", "SPI1",
  "C1QA", "C1QB", "C1QC", "AIF1", "TYROBP", "FCER1G"
)

source_files <- setNames(
  file.path(per_study_object_dir, paste0(vapply(study_order, safe_name, character(1)), "_full_count_processed.qs")),
  study_order
)
missing_inputs <- source_files[!file.exists(source_files)]
if (length(missing_inputs)) stop("Missing processed per-study object(s):\n", paste(missing_inputs, collapse = "\n"))

if (!file.exists(existing_mg_selection)) {
  stop("The existing Microglia selection audit file is missing: ", existing_mg_selection)
}
mg_selection_input <- read.csv(existing_mg_selection, check.names = FALSE, stringsAsFactors = FALSE)
write.csv(mg_selection_input,
          file.path(table_dir, "microglia_selection_input_read_only_copy.csv"),
          row.names = FALSE, na = "")

review_rows <- list()
integrity_rows <- list()
figure_rows <- list()
figure_index <- 0L

for (study in study_order) {
  message("\n===== OL parent review: ", study, " =====")
  source_file <- source_files[[study]]
  seu <- qs::qread(source_file)
  if (!inherits(seu, "Seurat")) stop(study, ": source is not a Seurat object.")
  # Explicit namespace avoids a restored-session conflict with
  # SummarizedExperiment::Assays.
  if (!"RNA" %in% SeuratObject::Assays(seu)) stop(study, ": no RNA assay.")
  SeuratObject::DefaultAssay(seu) <- "RNA"
  meta <- seu[[]]
  meta$Cell <- rownames(meta)
  required_columns <- c("analysis_cluster", "batch_id", "DoubletFinder_call")
  absent_columns <- setdiff(required_columns, colnames(meta))
  if (length(absent_columns)) {
    stop(study, ": required metadata absent: ", paste(absent_columns, collapse = ", "))
  }
  if (!"umap" %in% names(seu@reductions)) stop(study, ": UMAP reduction is absent.")
  counts <- get_assay_layer(seu, layer = "counts")
  normalized <- get_assay_layer(seu, layer = "data")
  assert_integer_counts(counts, study)
  if (!identical(colnames(counts), rownames(meta))) stop(study, ": count columns and metadata rows differ.")

  available_ol <- intersect(ol_markers, rownames(normalized))
  available_mg <- intersect(mg_markers, rownames(normalized))
  if (length(available_ol) < 3L || length(available_mg) < 3L) {
    stop(study, ": too few independent OL/MG markers found after feature-name matching.")
  }

  cluster_chr <- as.character(meta$analysis_cluster)
  cluster_levels <- cluster_order(cluster_chr)
  seu$analysis_cluster <- factor(cluster_chr, levels = cluster_levels)

  # 1. Cluster UMAP; it is the primary map used by the reviewer to nominate OL.
  p_cluster <- DimPlot(
    seu, reduction = "umap", group.by = "analysis_cluster", label = TRUE,
    repel = TRUE, pt.size = 0.18, raster = TRUE
  ) + ggtitle(paste0(study, ": post-DoubletFinder-singlet clusters")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  stem <- file.path(figure_dir, paste0(safe_name(study), "_01_cluster_UMAP"))
  save_plot_both(p_cluster, stem, width = 10, height = 8)
  print(p_cluster)
  figure_index <- figure_index + 1L
  figure_rows[[figure_index]] <- data.frame(Study = study, figure = "cluster_UMAP", stem = stem, stringsAsFactors = FALSE)

  # 2–3. Separate DotPlots retain readable independent OL and MG evidence.
  p_ol_dot <- DotPlot(seu, features = available_ol, group.by = "analysis_cluster") +
    RotatedAxis() + scale_color_gradient(low = "grey90", high = "#B2182B") +
    ggtitle(paste0(study, ": independent mature-OL markers")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  stem <- file.path(figure_dir, paste0(safe_name(study), "_02_OL_independent_marker_DotPlot"))
  save_plot_both(p_ol_dot, stem, width = 10, height = 7)
  print(p_ol_dot)
  figure_index <- figure_index + 1L
  figure_rows[[figure_index]] <- data.frame(Study = study, figure = "OL_marker_DotPlot", stem = stem, stringsAsFactors = FALSE)

  p_mg_dot <- DotPlot(seu, features = available_mg, group.by = "analysis_cluster") +
    RotatedAxis() + scale_color_gradient(low = "grey90", high = "#2166AC") +
    ggtitle(paste0(study, ": independent canonical-MG markers")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  stem <- file.path(figure_dir, paste0(safe_name(study), "_03_MG_independent_marker_DotPlot"))
  save_plot_both(p_mg_dot, stem, width = 12, height = 7)
  print(p_mg_dot)
  figure_index <- figure_index + 1L
  figure_rows[[figure_index]] <- data.frame(Study = study, figure = "MG_marker_DotPlot", stem = stem, stringsAsFactors = FALSE)

  # 4–5. FeaturePlot grids provide spatial confirmation independent of cluster averages.
  p_ol_feature <- FeaturePlot(
    seu, features = available_ol, reduction = "umap", ncol = 4,
    min.cutoff = "q05", max.cutoff = "q95", raster = TRUE
  ) + plot_annotation(title = paste0(study, ": independent mature-OL marker UMAPs")) &
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  stem <- file.path(figure_dir, paste0(safe_name(study), "_04_OL_independent_marker_FeaturePlot"))
  save_plot_both(p_ol_feature, stem, width = 14, height = 8)
  print(p_ol_feature)
  figure_index <- figure_index + 1L
  figure_rows[[figure_index]] <- data.frame(Study = study, figure = "OL_marker_FeaturePlot", stem = stem, stringsAsFactors = FALSE)

  p_mg_feature <- FeaturePlot(
    seu, features = available_mg, reduction = "umap", ncol = 4,
    min.cutoff = "q05", max.cutoff = "q95", raster = TRUE
  ) + plot_annotation(title = paste0(study, ": independent canonical-MG marker UMAPs")) &
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  stem <- file.path(figure_dir, paste0(safe_name(study), "_05_MG_independent_marker_FeaturePlot"))
  save_plot_both(p_mg_feature, stem, width = 14, height = 10)
  print(p_mg_feature)
  figure_index <- figure_index + 1L
  figure_rows[[figure_index]] <- data.frame(Study = study, figure = "MG_marker_FeaturePlot", stem = stem, stringsAsFactors = FALSE)

  # Compact marker summaries are calculated only for the displayed markers.
  # No full matrix is densified, keeping the review reproducible on the server.
  for (cl in cluster_levels) {
    cells <- meta$Cell[cluster_chr == cl]
    ol_matrix <- normalized[available_ol, cells, drop = FALSE]
    mg_matrix <- normalized[available_mg, cells, drop = FALSE]
    ol_score <- mean(Matrix::rowMeans(ol_matrix))
    mg_score <- mean(Matrix::rowMeans(mg_matrix))
    ol_detect <- mean(Matrix::rowMeans(ol_matrix > 0))
    mg_detect <- mean(Matrix::rowMeans(mg_matrix > 0))
    review_rows[[length(review_rows) + 1L]] <- data.frame(
      Study = study,
      analysis_cluster = cl,
      n_cells = length(cells),
      n_batches = length(unique(as.character(meta$batch_id[cluster_chr == cl]))),
      median_nCount_RNA = if ("nCount_RNA" %in% colnames(meta)) median(meta$nCount_RNA[cluster_chr == cl], na.rm = TRUE) else NA_real_,
      median_nFeature_RNA = if ("nFeature_RNA" %in% colnames(meta)) median(meta$nFeature_RNA[cluster_chr == cl], na.rm = TRUE) else NA_real_,
      median_percent.mt = if ("percent.mt" %in% colnames(meta)) median(meta$percent.mt[cluster_chr == cl], na.rm = TRUE) else NA_real_,
      OL_independent_score = ol_score,
      MG_independent_score = mg_score,
      OL_minus_MG = ol_score - mg_score,
      OL_marker_fraction_detected = ol_detect,
      MG_marker_fraction_detected = mg_detect,
      OL_markers_present = paste(available_ol, collapse = ";"),
      MG_markers_present = paste(available_mg, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }

  integrity_rows[[length(integrity_rows) + 1L]] <- data.frame(
    Study = study,
    source_object = source_file,
    n_features = nrow(counts),
    n_cells = ncol(counts),
    count_nnz = length(counts@x),
    min_nonzero_count = if (length(counts@x)) min(counts@x) else NA_real_,
    max_count = if (length(counts@x)) max(counts@x) else NA_real_,
    all_counts_integer = TRUE,
    all_current_DoubletFinder_singlet = all(as.character(meta$DoubletFinder_call) == "Singlet"),
    n_batches = length(unique(as.character(meta$batch_id))),
    n_analysis_clusters = length(cluster_levels),
    stringsAsFactors = FALSE
  )

  rm(seu, meta, counts, normalized)
  gc(verbose = FALSE)
}

review_summary <- do.call(rbind, review_rows)
review_summary <- review_summary[order(match(review_summary$Study, study_order),
                                      suppressWarnings(as.numeric(review_summary$analysis_cluster)),
                                      review_summary$analysis_cluster), , drop = FALSE]
write.csv(review_summary, file.path(table_dir, "OL_parent_cluster_marker_review_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, integrity_rows), file.path(table_dir, "RNA_integer_count_integrity_audit.csv"), row.names = FALSE)

selection_template <- file.path(manual_dir, "OL_parent_cluster_selection_TEMPLATE.csv")
selection_base <- review_summary[, c(
  "Study", "analysis_cluster", "n_cells", "n_batches", "OL_independent_score",
  "MG_independent_score", "OL_minus_MG", "OL_marker_fraction_detected",
  "MG_marker_fraction_detected"
)]
selection_base$select_as_mature_OL <- NA
selection_base$review_rationale <- NA_character_
selection_base$reviewer <- NA_character_
selection_base$review_date <- NA_character_
if (!file.exists(selection_template)) {
  write.csv(selection_base, selection_template, row.names = FALSE, na = "")
  message("Created empty OL parent selection template: ", selection_template)
} else {
  existing_selection <- read.csv(selection_template, check.names = FALSE, stringsAsFactors = FALSE)
  required_selection_cols <- c("Study", "analysis_cluster", "select_as_mature_OL")
  missing_selection_cols <- setdiff(required_selection_cols, colnames(existing_selection))
  if (length(missing_selection_cols)) {
    stop("Existing OL selection table has missing required columns: ", paste(missing_selection_cols, collapse = ", "))
  }
  expected_keys <- paste(selection_base$Study, selection_base$analysis_cluster, sep = "||")
  existing_keys <- paste(existing_selection$Study, existing_selection$analysis_cluster, sep = "||")
  table_matches_current_clusters <- identical(sort(expected_keys), sort(existing_keys))
  empty_flag <- is.na(existing_selection$select_as_mature_OL) |
    trimws(as.character(existing_selection$select_as_mature_OL)) == ""
  empty_rationale <- !"review_rationale" %in% colnames(existing_selection) |
    is.na(existing_selection$review_rationale) |
    trimws(as.character(existing_selection$review_rationale)) == ""
  # A template generated by an interrupted/older script is safe to replace only
  # when it contains no reviewer decisions.  Any nonempty reviewer content is
  # never overwritten; a mismatching table then stops the run for resolution.
  if (!table_matches_current_clusters && all(empty_flag) && all(empty_rationale)) {
    write.csv(existing_selection,
              file.path(manual_dir, "OL_parent_cluster_selection_TEMPLATE_invalid_pre_repair.csv"),
              row.names = FALSE, na = "")
    write.csv(selection_base, selection_template, row.names = FALSE, na = "")
    message("Replaced an empty invalid template; its pre-repair copy was retained for audit.")
  } else if (!table_matches_current_clusters) {
    stop("Existing OL selection table does not match the current study/cluster set and contains reviewer content; it was not overwritten.")
  } else {
    write.csv(existing_selection, file.path(manual_dir, "OL_parent_cluster_selection_CURRENT.csv"), row.names = FALSE, na = "")
    message("Existing OL parent selection template was preserved; no reviewer entries were overwritten.")
  }
}

figure_manifest <- do.call(rbind, figure_rows)
figure_manifest$png_exists <- file.exists(paste0(figure_manifest$stem, ".png"))
figure_manifest$pdf_exists <- file.exists(paste0(figure_manifest$stem, ".pdf"))
if (!all(figure_manifest$png_exists & figure_manifest$pdf_exists)) {
  stop("At least one review figure does not have a complete PNG/PDF pair.")
}
write.csv(figure_manifest, file.path(table_dir, "figure_PNG_PDF_pair_manifest.csv"), row.names = FALSE)

script_name <- "23.OL_parent_cluster_review_24server.r"
manifest <- data.frame(
  key = c(
    "script", "run_timestamp", "seed", "input_objects", "existing_microglia_selection_read_only",
    "OL_primary_markers", "MG_primary_markers", "explicitly_excluded_from_primary_panel",
    "selection_status", "real_cell_objects_modified", "artificial_profiles_created", "n_png_pdf_figure_pairs"
  ),
  value = c(
    script_name,
    format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), "20260909",
    paste(source_files, collapse = ";"), existing_mg_selection,
    paste(ol_markers, collapse = ";"), paste(mg_markers, collapse = ";"),
    "TRIM2;RNF220;ST18;OLIG2;MBP;MOG",
    "Awaiting manual mature-OL cluster selection", "FALSE", "FALSE", nrow(figure_manifest)
  ),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(runtime_dir, "stage1_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "sessionInfo.txt"))

message("\nStage 1 completed.")
message("Review figures: ", figure_dir)
message("Please populate ONLY select_as_mature_OL (TRUE/FALSE) and review_rationale in: ", selection_template)
message("No OL parent clusters have been selected automatically, and no artificial doublets have been created.")
