#!/usr/bin/env Rscript

# Stage 1 of the full count-based reanalysis. Each complete study is rebuilt
# from integer RNA counts; source Seurat processing is never reused. The script
# pauses after generating cluster/marker evidence for USER selection of the
# microglia cluster(s) in each study.

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(harmony)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(qs)
})

set.seed(20260903)
project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
result_root <- file.path(project_root, "combine/results/full_count_reanalysis/01_per_study")
figure_root <- file.path(result_root, "figure")
table_root <- file.path(result_root, "table")
object_root <- file.path(result_root, "object")
runtime_root <- file.path(result_root, "runtime")
selection_root <- file.path(project_root, "combine/results/full_count_reanalysis/manual_selections")
for (d in c(figure_root, table_root, object_root, runtime_root, selection_root)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# This analysis is deliberately started from persisted count objects, not from
# anything in the interactive workspace. Release only large artifacts left by
# prior plotting/diagnostic steps before reading a source dataset. This keeps
# the foreground RStudio session below its practical memory limit; the removed
# objects already have their corresponding saved results on /data4.
workspace_names <- ls(envir = .GlobalEnv, all.names = TRUE)
workspace_bytes <- vapply(workspace_names, function(nm) {
  as.numeric(utils::object.size(get(nm, envir = .GlobalEnv, inherits = FALSE)))
}, numeric(1))
stale_large_names <- names(workspace_bytes)[workspace_bytes >= 100 * 1024^2]
cleanup_audit <- data.frame(object = stale_large_names,
                            bytes_released = unname(workspace_bytes[stale_large_names]),
                            stringsAsFactors = FALSE)
utils::write.csv(cleanup_audit, file.path(runtime_root, "workspace_memory_cleanup.csv"), row.names = FALSE)
if (length(stale_large_names)) {
  rm(list = stale_large_names, envir = .GlobalEnv)
  invisible(gc())
}

# Batch fields were verified in the 00_preflight audit. Yao currently contains
# one recorded technical batch, so Harmony is intentionally skipped there.
input_manifest <- data.frame(
  Study = c("Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
            "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"),
  Source_object = c(
    file.path(project_root, "Su_et_al./data/processed_data/Su_CSC_alldata_v2.qs"),
    file.path(project_root, "Chia_et_al./data/processed_data/Chia.qs"),
    file.path(project_root, "Fullard_et_al./data/processed_data/Fullard.qs"),
    file.path(project_root, "Jakel_et_al./data/processed_data/Jakel.qs"),
    file.path(project_root, "Grubman_et_al./data/processed_data/Grubman.qs"),
    file.path(project_root, "Yao_et_al./data/processed_data/Yao_stroke.qs")
  ),
  Batch_field = c("Sample", "sample", "subject_id", "sample", "patient", "Batch"),
  stringsAsFactors = FALSE
)

microglia_markers <- c("P2RY12", "TMEM119", "CX3CR1", "AIF1", "C1QA", "C1QB", "C1QC", "GPR34",
                       "SALL1", "CSF1R", "SPI1", "TYROBP", "FCER1G", "TREM2", "APOE", "LPL", "CTSD")

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 320,
         bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
         device = grDevices::pdf, bg = "white", limitsize = FALSE)
}

extract_counts <- function(obj) {
  if (!"RNA" %in% Seurat::Assays(obj)) stop("RNA assay is absent.")
  assay <- obj[["RNA"]]
  counts <- tryCatch(assay@counts, error = function(e) NULL)
  if (is.null(counts) || !length(counts)) {
    counts <- tryCatch(Seurat::GetAssayData(obj, assay = "RNA", slot = "counts"), error = function(e) NULL)
  }
  if (is.null(counts) || !length(counts)) stop("RNA counts slot/layer is unavailable.")
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  if (length(counts@x) && any(abs(counts@x - round(counts@x)) > 1e-8)) stop("RNA count matrix is not integer-valued.")
  counts
}

make_thresholds <- function(metadata, mt_available) {
  med_ncount <- median(metadata$nCount_RNA, na.rm = TRUE)
  med_nfeature <- median(metadata$nFeature_RNA, na.rm = TRUE)
  mad_ncount <- mad(metadata$nCount_RNA, constant = 1.4826, na.rm = TRUE)
  mad_nfeature <- mad(metadata$nFeature_RNA, constant = 1.4826, na.rm = TRUE)
  max_mt <- if (!mt_available) Inf else min(25, max(10, ceiling(median(metadata$percent.mt, na.rm = TRUE) + 3 * mad(metadata$percent.mt, constant = 1.4826, na.rm = TRUE))))
  c(min_nCount_RNA = max(200, floor(med_ncount - 3 * mad_ncount)),
    min_nFeature_RNA = max(200, floor(med_nfeature - 3 * mad_nfeature)),
    max_percent_mt = max_mt)
}

run_doubletfinder_partition <- function(counts, metadata, partition_name) {
  if (ncol(counts) < 300L) stop("DoubletFinder partition ", partition_name, " has fewer than 300 QC-passed cells.")
  obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = metadata, project = partition_name)
  obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  n_hvg <- min(3000L, nrow(obj) - 1L)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = n_hvg, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(obj), vars.to.regress = "percent.mt", verbose = FALSE)
  n_pcs <- min(30L, length(Seurat::VariableFeatures(obj)) - 1L, ncol(obj) - 1L)
  obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(obj), npcs = n_pcs, seed.use = 20260903, verbose = FALSE)
  dims_use <- seq_len(min(30L, ncol(Seurat::Embeddings(obj, "pca"))))
  obj <- Seurat::FindNeighbors(obj, dims = dims_use, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = 0.4, random.seed = 20260903, verbose = FALSE)
  obj$DF_cluster <- as.character(obj$seurat_clusters)
  expected_rate <- min(0.20, 0.008 * ncol(obj) / 1000)
  # pN/pK scanning has a very high computational cost. pK is therefore tuned
  # on a deterministic representative subset of at most 3,000 cells, then the
  # selected pK is used for the formal DoubletFinder call on every cell in the
  # technical batch. This preserves full-batch classification while keeping
  # the six-study analysis feasible in a foreground RStudio session.
  tuning_cells <- colnames(obj)
  if (length(tuning_cells) > 3000L) {
    set.seed(20260903)
    tuning_cells <- sample(tuning_cells, 3000L)
  }
  tuning_obj <- subset(obj, cells = tuning_cells)
  sweep <- DoubletFinder::paramSweep(tuning_obj, PCs = dims_use, sct = FALSE)
  sweep_stats <- DoubletFinder::summarizeSweep(sweep, GT = FALSE)
  pk_table <- DoubletFinder::find.pK(sweep_stats)
  pk_values <- suppressWarnings(as.numeric(as.character(pk_table$pK)))
  best_pk <- pk_values[which.max(pk_table$BCmetric)]
  if (!is.finite(best_pk)) stop("No numeric optimal pK for ", partition_name)
  homotypic_fraction <- DoubletFinder::modelHomotypic(obj$DF_cluster)
  n_exp_raw <- max(1L, round(expected_rate * ncol(obj)))
  n_exp_adjusted <- max(1L, round(n_exp_raw * (1 - homotypic_fraction)))
  metadata_columns_before_df <- colnames(obj[[]])
  obj <- DoubletFinder::doubletFinder(obj, PCs = dims_use, pN = 0.25, pK = best_pk,
                                      nExp = n_exp_adjusted, reuse.pANN = NULL, sct = FALSE)
  md <- obj[[]]
  # Some author-provided metadata already contain historical DoubletFinder
  # columns. Select only the two columns created by this analysis call.
  new_df_columns <- setdiff(colnames(md), metadata_columns_before_df)
  class_column <- grep("^DF.classifications_", new_df_columns, value = TRUE)
  score_column <- grep("^pANN_", new_df_columns, value = TRUE)
  if (length(class_column) != 1L || length(score_column) != 1L) {
    stop("Cannot identify newly created DoubletFinder output for ", partition_name,
         "; new metadata columns: ", paste(new_df_columns, collapse = ", "))
  }
  # Seurat metadata columns are vectors in the installed DoubletFinder/Seurat
  # combination (rather than one-column matrices), so index them as vectors.
  calls <- data.frame(Cell = rownames(md), DoubletFinder_score = as.numeric(md[[score_column]]),
    DoubletFinder_call = ifelse(grepl("Doublet", as.character(md[[class_column]])), "Doublet", "Singlet"),
    DF_partition = partition_name, stringsAsFactors = FALSE)
  parameters <- data.frame(DF_partition = partition_name, n_cells = ncol(obj), expected_doublet_rate = expected_rate,
    pN = 0.25, pK = best_pk, nExp_raw = n_exp_raw, nExp_homotypic_adjusted = n_exp_adjusted,
    homotypic_fraction = homotypic_fraction, pK_tuning_n_cells = length(tuning_cells),
    pK_tuned_by_subsample = ncol(obj) > length(tuning_cells), stringsAsFactors = FALSE)
  rm(obj, tuning_obj, sweep, sweep_stats, pk_table)
  gc(verbose = FALSE)
  list(calls = calls, parameters = parameters)
}

plot_qc_metric <- function(qc_cells, thresholds, metric, threshold_column, direction) {
  p <- ggplot(qc_cells, aes(x = .data[[metric]])) +
    geom_histogram(bins = 70, fill = "#9ecae1", colour = "white", linewidth = 0.1) +
    facet_wrap(~batch_id, scales = "free", ncol = 4) +
    labs(title = paste0(metric, " before QC filtering"), subtitle = paste0("Red line: ", direction, " QC threshold selected per technical batch"),
         x = metric, y = "Cells") + theme_bw(base_size = 9) + theme(strip.text = element_text(face = "bold"))
  threshold_data <- thresholds[is.finite(thresholds[[threshold_column]]), c("batch_id", threshold_column), drop = FALSE]
  if (nrow(threshold_data)) p <- p + geom_vline(data = threshold_data, aes(xintercept = .data[[threshold_column]]), colour = "#d7301f", linewidth = 0.4)
  p
}

process_study <- function(study, source_file, batch_field) {
  study_safe <- safe_name(study)
  checkpoint_file <- file.path(object_root, paste0(study_safe, "_full_count_processed.qs"))
  if (file.exists(checkpoint_file)) {
    message("Checkpoint exists; retaining prior completed study: ", study)
    return(list(summary = read.csv(file.path(table_root, paste0(study_safe, "_cluster_summary.csv")), stringsAsFactors = FALSE),
                marker_summary = read.csv(file.path(table_root, paste0(study_safe, "_top_markers.csv")), stringsAsFactors = FALSE)))
  }
  message("\n===== Full count reanalysis: ", study, " =====")
  raw_object <- qs::qread(source_file)
  counts <- extract_counts(raw_object)
  metadata <- raw_object[[]]
  # Keep the immutable integer count matrix and metadata only. The source
  # Seurat object can also carry prior reductions/graphs/normalized layers,
  # which are intentionally not reused and would otherwise inflate peak RAM.
  rm(raw_object)
  invisible(gc())
  if (!batch_field %in% colnames(metadata)) stop(study, ": required batch field is absent: ", batch_field)
  # Retain any author-supplied historical doublet calls as provenance, but
  # rename them so this analysis has an unambiguous pANN/DF output namespace.
  author_df_columns <- grep("^(pANN_|DF.classifications_)", colnames(metadata), value = TRUE)
  if (length(author_df_columns)) {
    colnames(metadata)[match(author_df_columns, colnames(metadata))] <- make.unique(paste0("author_", author_df_columns))
  }
  metadata <- metadata[colnames(counts), , drop = FALSE]
  original_barcode <- colnames(counts)
  new_cell_ids <- make.unique(paste(study_safe, original_barcode, sep = "__"))
  colnames(counts) <- new_cell_ids
  rownames(metadata) <- new_cell_ids
  metadata$original_barcode <- original_barcode
  metadata$Study <- study
  metadata$batch_id <- as.character(metadata[[batch_field]])
  metadata$batch_id[is.na(metadata$batch_id) | metadata$batch_id == ""] <- "missing_batch"
  mt_genes <- grep("^MT-", rownames(counts), value = TRUE)
  ribo_genes <- grep("^RP[SL]", rownames(counts), value = TRUE)
  metadata$nCount_RNA <- Matrix::colSums(counts)
  metadata$nFeature_RNA <- Matrix::colSums(counts > 0)
  metadata$percent.mt <- if (length(mt_genes)) Matrix::colSums(counts[mt_genes, , drop = FALSE]) / pmax(metadata$nCount_RNA, 1) * 100 else 0
  metadata$percent.ribo <- if (length(ribo_genes)) Matrix::colSums(counts[ribo_genes, , drop = FALSE]) / pmax(metadata$nCount_RNA, 1) * 100 else 0

  threshold_rows <- lapply(split(rownames(metadata), metadata$batch_id), function(cells) {
    threshold <- make_thresholds(metadata[cells, , drop = FALSE], mt_available = length(mt_genes) > 0)
    data.frame(batch_id = metadata[cells[[1]], "batch_id"], n_before = length(cells),
      min_nCount_RNA = threshold[["min_nCount_RNA"]], min_nFeature_RNA = threshold[["min_nFeature_RNA"]],
      max_percent_mt = threshold[["max_percent_mt"]], stringsAsFactors = FALSE)
  })
  thresholds <- do.call(rbind, threshold_rows)
  threshold_index <- match(metadata$batch_id, thresholds$batch_id)
  metadata$passes_basic_QC <- metadata$nCount_RNA >= thresholds$min_nCount_RNA[threshold_index] &
    metadata$nFeature_RNA >= thresholds$min_nFeature_RNA[threshold_index] &
    metadata$percent.mt <= thresholds$max_percent_mt[threshold_index]
  thresholds$n_after_basic_QC <- vapply(thresholds$batch_id, function(x) sum(metadata$passes_basic_QC & metadata$batch_id == x), numeric(1))
  thresholds$n_removed_basic_QC <- thresholds$n_before - thresholds$n_after_basic_QC
  if (sum(metadata$passes_basic_QC) < 500L) stop(study, ": fewer than 500 cells passed basic QC.")

  study_figure_dir <- file.path(figure_root, study_safe)
  dir.create(study_figure_dir, recursive = TRUE, showWarnings = FALSE)
  p_qc_ncount <- plot_qc_metric(metadata, thresholds, "nCount_RNA", "min_nCount_RNA", "lower")
  p_qc_nfeature <- plot_qc_metric(metadata, thresholds, "nFeature_RNA", "min_nFeature_RNA", "lower")
  p_qc_mt <- plot_qc_metric(metadata, thresholds, "percent.mt", "max_percent_mt", "upper")
  qc_plot_list <- list(nCount_RNA = p_qc_ncount, nFeature_RNA = p_qc_nfeature, percent_mt = p_qc_mt)
  for (plot_name in names(qc_plot_list)) {
    save_plot_both(qc_plot_list[[plot_name]], file.path(study_figure_dir, paste0("QC_", plot_name)), 15, 10)
    print(qc_plot_list[[plot_name]])
  }

  qc_cells <- rownames(metadata)[metadata$passes_basic_QC]
  qc_counts <- counts[, qc_cells, drop = FALSE]
  qc_metadata <- metadata[qc_cells, , drop = FALSE]
  # DoubletFinder needs enough cells for its pK sweep. Very small technical
  # batches are pooled only for doublet calling; their original batch_id is
  # retained for all subsequent QC and Harmony steps.
  batch_sizes <- table(qc_metadata$batch_id)
  qc_metadata$DF_partition_for_run <- qc_metadata$batch_id
  small_batch_names <- names(batch_sizes)[batch_sizes < 300L]
  if (length(small_batch_names)) {
    qc_metadata$DF_partition_for_run[qc_metadata$batch_id %in% small_batch_names] <- "pooled_small_batches"
    # A collection of small batches can itself still be too small for the
    # DoubletFinder pK sweep. In that case, add the largest valid technical
    # batch to this *calling-only* partition. The original batch_id is never
    # changed, so later QC reporting and Harmony retain the real provenance.
    if (sum(batch_sizes[small_batch_names]) < 300L) {
      reference_batch <- names(batch_sizes)[which.max(batch_sizes)]
      fallback_partition <- paste0("pooled_small_batches_with_", safe_name(reference_batch))
      qc_metadata$DF_partition_for_run[qc_metadata$batch_id %in% c(small_batch_names, reference_batch)] <- fallback_partition
    }
  }
  df_partition_mapping <- unique(qc_metadata[, c("batch_id", "DF_partition_for_run"), drop = FALSE])
  rownames(df_partition_mapping) <- NULL
  partition_cells <- split(qc_cells, qc_metadata$DF_partition_for_run)
  df_results <- lapply(names(partition_cells), function(partition_name) {
    cells <- partition_cells[[partition_name]]
    run_doubletfinder_partition(qc_counts[, cells, drop = FALSE], qc_metadata[cells, , drop = FALSE], partition_name)
  })
  df_calls <- do.call(rbind, lapply(df_results, `[[`, "calls"))
  df_parameters <- do.call(rbind, lapply(df_results, `[[`, "parameters"))
  rownames(df_calls) <- df_calls$Cell
  qc_metadata$DoubletFinder_score <- df_calls[rownames(qc_metadata), "DoubletFinder_score"]
  qc_metadata$DoubletFinder_call <- factor(df_calls[rownames(qc_metadata), "DoubletFinder_call"], levels = c("Singlet", "Doublet"))
  qc_metadata$DF_partition <- df_calls[rownames(qc_metadata), "DF_partition"]

  diagnostic <- Seurat::CreateSeuratObject(counts = qc_counts, meta.data = qc_metadata, project = study)
  diagnostic <- Seurat::NormalizeData(diagnostic, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  diagnostic <- Seurat::FindVariableFeatures(diagnostic, selection.method = "vst", nfeatures = min(3000L, nrow(diagnostic) - 1L), verbose = FALSE)
  diagnostic <- Seurat::ScaleData(diagnostic, features = Seurat::VariableFeatures(diagnostic), vars.to.regress = "percent.mt", verbose = FALSE)
  diagnostic <- Seurat::RunPCA(diagnostic, features = Seurat::VariableFeatures(diagnostic), npcs = min(30L, length(Seurat::VariableFeatures(diagnostic)) - 1L), seed.use = 20260903, verbose = FALSE)
  diagnostic_dims <- seq_len(min(30L, ncol(Seurat::Embeddings(diagnostic, "pca"))))
  diagnostic <- Seurat::FindNeighbors(diagnostic, dims = diagnostic_dims, verbose = FALSE)
  diagnostic <- Seurat::FindClusters(diagnostic, resolution = 0.4, random.seed = 20260903, verbose = FALSE)
  diagnostic$diagnostic_cluster <- diagnostic$seurat_clusters
  diagnostic <- Seurat::RunUMAP(diagnostic, dims = diagnostic_dims, reduction.name = "umap_diagnostic", seed.use = 20260903, verbose = FALSE)
  qs::qsave(diagnostic, file.path(object_root, paste0(study_safe, "_QC_before_doublet_removal.qs")), preset = "high")
  p_diagnostic_cluster <- Seurat::DimPlot(diagnostic, reduction = "umap_diagnostic", group.by = "diagnostic_cluster", label = TRUE, repel = TRUE) + ggtitle(paste0(study, ": QC-passed diagnostic UMAP"))
  p_diagnostic_df <- Seurat::DimPlot(diagnostic, reduction = "umap_diagnostic", group.by = "DoubletFinder_call", cols = c(Singlet = "#bdbdbd", Doublet = "#d7301f")) + ggtitle(paste0(study, ": DoubletFinder calls before removal"))
  save_plot_both(p_diagnostic_cluster, file.path(study_figure_dir, "Diagnostic_UMAP_clusters_pre_doublet"), 8, 6)
  save_plot_both(p_diagnostic_df, file.path(study_figure_dir, "Diagnostic_UMAP_DoubletFinder_pre_doublet"), 8, 6)
  print(p_diagnostic_cluster)
  print(p_diagnostic_df)

  retained_cells <- rownames(qc_metadata)[qc_metadata$DoubletFinder_call == "Singlet"]
  if (length(retained_cells) < 500L) stop(study, ": fewer than 500 cells remain after DoubletFinder removal.")
  final_object <- Seurat::CreateSeuratObject(counts = qc_counts[, retained_cells, drop = FALSE], meta.data = qc_metadata[retained_cells, , drop = FALSE], project = study)
  final_object <- Seurat::NormalizeData(final_object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  final_object <- Seurat::FindVariableFeatures(final_object, selection.method = "vst", nfeatures = min(3000L, nrow(final_object) - 1L), verbose = FALSE)
  final_object <- Seurat::ScaleData(final_object, features = Seurat::VariableFeatures(final_object), vars.to.regress = "percent.mt", verbose = FALSE)
  final_object <- Seurat::RunPCA(final_object, features = Seurat::VariableFeatures(final_object), npcs = min(30L, length(Seurat::VariableFeatures(final_object)) - 1L), seed.use = 20260903, verbose = FALSE)
  final_dims <- seq_len(min(30L, ncol(Seurat::Embeddings(final_object, "pca"))))
  use_reduction <- "pca"
  batch_note <- "Single recorded batch: no within-study batch correction applied."
  if (length(unique(final_object$batch_id)) > 1L) {
    final_object <- harmony::RunHarmony(final_object, group.by.vars = "batch_id", reduction = "pca", dims.use = final_dims,
                                        reduction.save = "harmony", verbose = FALSE)
    use_reduction <- "harmony"
    batch_note <- "Harmony correction by verified technical batch."
  }
  final_object <- Seurat::FindNeighbors(final_object, reduction = use_reduction, dims = final_dims, verbose = FALSE)
  final_object <- Seurat::FindClusters(final_object, resolution = 0.4, random.seed = 20260903, verbose = FALSE)
  final_object <- Seurat::RunUMAP(final_object, reduction = use_reduction, dims = final_dims, seed.use = 20260903, verbose = FALSE)
  final_object$analysis_cluster <- final_object$seurat_clusters
  marker_table <- Seurat::FindAllMarkers(final_object, assay = "RNA", only.pos = TRUE, min.pct = 0.10,
    logfc.threshold = 0.25, max.cells.per.ident = 3000, random.seed = 20260903, verbose = FALSE)
  fold_change_column <- intersect(c("avg_log2FC", "avg_logFC"), colnames(marker_table))[[1]]
  top_markers <- do.call(rbind, lapply(split(marker_table, marker_table$cluster), function(x) head(x[order(x$p_val_adj, -x[[fold_change_column]]), , drop = FALSE], 30)))
  if (is.null(top_markers)) top_markers <- marker_table[0, , drop = FALSE]
  cluster_summary <- as.data.frame(table(Study = rep(study, ncol(final_object)), analysis_cluster = final_object$analysis_cluster), stringsAsFactors = FALSE)
  colnames(cluster_summary)[3] <- "n_cells"
  cluster_summary$batch_correction <- batch_note

  p_final_cluster <- Seurat::DimPlot(final_object, group.by = "analysis_cluster", label = TRUE, repel = TRUE) + ggtitle(paste0(study, ": post-doublet, ", use_reduction, " UMAP"))
  p_final_batch <- Seurat::DimPlot(final_object, group.by = "batch_id") + ggtitle(paste0(study, ": batch distribution after ", use_reduction))
  available_mg_markers <- intersect(microglia_markers, rownames(final_object))
  p_mg_dotplot <- Seurat::DotPlot(final_object, features = available_mg_markers, group.by = "analysis_cluster") + RotatedAxis() + ggtitle(paste0(study, ": canonical microglia marker dot plot"))
  p_mg_feature <- Seurat::FeaturePlot(final_object, features = available_mg_markers, ncol = 4, min.cutoff = "q05", max.cutoff = "q95") + plot_annotation(title = paste0(study, ": canonical microglia marker expression"))
  final_plot_list <- list(Final_UMAP_clusters = p_final_cluster, Final_UMAP_batch = p_final_batch,
    Microglia_marker_dotplot = p_mg_dotplot, Microglia_marker_UMAP = p_mg_feature)
  for (plot_name in names(final_plot_list)) {
    save_plot_both(final_plot_list[[plot_name]], file.path(study_figure_dir, plot_name), 12, 8)
    print(final_plot_list[[plot_name]])
  }

  per_cell <- final_object[[]]
  per_cell$Cell <- rownames(per_cell)
  utils::write.csv(thresholds, file.path(table_root, paste0(study_safe, "_QC_thresholds_and_retention.csv")), row.names = FALSE)
  utils::write.csv(metadata, file.path(table_root, paste0(study_safe, "_all_cells_QC_metrics.csv")), row.names = TRUE)
  utils::write.csv(df_partition_mapping, file.path(table_root, paste0(study_safe, "_DoubletFinder_partition_mapping.csv")), row.names = FALSE)
  utils::write.csv(df_calls, file.path(table_root, paste0(study_safe, "_DoubletFinder_per_cell.csv")), row.names = FALSE)
  utils::write.csv(df_parameters, file.path(table_root, paste0(study_safe, "_DoubletFinder_parameters.csv")), row.names = FALSE)
  utils::write.csv(marker_table, file.path(table_root, paste0(study_safe, "_all_cluster_markers.csv")), row.names = FALSE)
  utils::write.csv(top_markers, file.path(table_root, paste0(study_safe, "_top_markers.csv")), row.names = FALSE)
  utils::write.csv(cluster_summary, file.path(table_root, paste0(study_safe, "_cluster_summary.csv")), row.names = FALSE)
  utils::write.csv(per_cell, file.path(table_root, paste0(study_safe, "_retained_singlet_metadata.csv")), row.names = FALSE)
  qs::qsave(final_object, checkpoint_file, preset = "high")
  rm(counts, metadata, qc_counts, qc_metadata, diagnostic, final_object, df_results, df_calls, df_parameters, marker_table)
  gc(verbose = FALSE)
  list(summary = cluster_summary, marker_summary = top_markers)
}

study_results <- Map(process_study, input_manifest$Study, input_manifest$Source_object, input_manifest$Batch_field)
cluster_summary_all <- do.call(rbind, lapply(study_results, `[[`, "summary"))
marker_summary_all <- do.call(rbind, lapply(study_results, `[[`, "marker_summary"))
utils::write.csv(input_manifest, file.path(table_root, "per_study_input_manifest.csv"), row.names = FALSE)
utils::write.csv(cluster_summary_all, file.path(table_root, "all_studies_cluster_summary.csv"), row.names = FALSE)
utils::write.csv(marker_summary_all, file.path(table_root, "all_studies_top_markers.csv"), row.names = FALSE)

# This is intentionally blank: it is the auditable hand-off point where the
# user records which evidence-supported cluster(s) should become Microglia.
selection_template <- cluster_summary_all[, c("Study", "analysis_cluster", "n_cells"), drop = FALSE]
selection_template$select_as_Microglia <- NA
selection_template$user_note <- NA_character_
utils::write.csv(selection_template, file.path(selection_root, "per_study_Microglia_cluster_selection_TEMPLATE.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(runtime_root, "per_study_R_sessionInfo.txt"))
writeLines(c(
  "Input: verified integer RNA counts from six complete study objects.",
  "QC: nCount_RNA, nFeature_RNA, and percent.mt thresholds selected independently per verified technical batch using median +/- 3 MAD; high UMI/gene outliers are not removed before DoubletFinder.",
  "Doublet removal: DoubletFinder per technical batch; pK is selected using a deterministic at-most-3,000-cell representative subset and the resulting pK is applied to the full technical batch. Score, call, pK, expected rate and homotypic adjustment are retained.",
  "Batch correction: Harmony by the verified technical batch field where multiple levels exist; skipped and explicitly recorded for single-batch Yao.",
  "Markers: FindAllMarkers on retained singlets, capped at 3,000 cells per cluster for tractable, reproducible testing; all cells remain in clustering/UMAP.",
  "Manual gate: this stage does not define Microglia. The user must populate the selection template before stage 2 extraction/integration."
), file.path(runtime_root, "per_study_method_manifest.txt"))
message("Stage 1 complete. Review per-study UMAP/dot plots and populate: ", file.path(selection_root, "per_study_Microglia_cluster_selection_TEMPLATE.csv"))
