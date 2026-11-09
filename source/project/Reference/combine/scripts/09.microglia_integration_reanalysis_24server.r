#!/usr/bin/env Rscript

# Stage 2 of the full count-based reanalysis. This script requires the user to
# complete the per-study Microglia cluster-selection CSV produced by stage 1.
# It never automatically chooses either Microglia or MG0.

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
stage1_root <- file.path(project_root, "combine/results/full_count_reanalysis/01_per_study")
selection_file <- file.path(project_root, "combine/results/full_count_reanalysis/manual_selections/per_study_Microglia_cluster_selection_TEMPLATE.csv")
result_root <- file.path(project_root, "combine/results/full_count_reanalysis/02_microglia_integration")
figure_root <- file.path(result_root, "figure")
table_root <- file.path(result_root, "table")
object_root <- file.path(result_root, "object")
runtime_root <- file.path(result_root, "runtime")
manual_root <- file.path(project_root, "combine/results/full_count_reanalysis/manual_selections")
for (d in c(figure_root, table_root, object_root, runtime_root, manual_root)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Previous figures/objects are persisted on /data4. Release only large
# interactive artifacts so this foreground RStudio job starts from a bounded
# memory footprint.
workspace_names <- ls(envir = .GlobalEnv, all.names = TRUE)
workspace_bytes <- vapply(workspace_names, function(nm) as.numeric(utils::object.size(get(nm, envir = .GlobalEnv, inherits = FALSE))), numeric(1))
stale_large_names <- names(workspace_bytes)[workspace_bytes >= 100 * 1024^2]
utils::write.csv(data.frame(object = stale_large_names, bytes_released = unname(workspace_bytes[stale_large_names])),
                 file.path(runtime_root, "workspace_memory_cleanup.csv"), row.names = FALSE)
if (length(stale_large_names)) { rm(list = stale_large_names, envir = .GlobalEnv); invisible(gc()) }

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = grDevices::pdf, bg = "white", limitsize = FALSE)
}
get_counts <- function(obj) {
  if (!"RNA" %in% Seurat::Assays(obj)) stop("RNA assay is unavailable.")
  # Seurat v5 stores raw counts in an Assay5 layer; older objects use a slot.
  # Try the current layer API first, then retain a legacy fallback for old objects.
  x <- tryCatch(SeuratObject::LayerData(obj[["RNA"]], layer = "counts"), error = function(e) NULL)
  if (is.null(x) || !length(x)) {
    x <- tryCatch(Seurat::GetAssayData(obj, assay = "RNA", layer = "counts"), error = function(e) NULL)
  }
  if (is.null(x) || !length(x)) {
    x <- tryCatch(Seurat::GetAssayData(obj, assay = "RNA", slot = "counts"), error = function(e) NULL)
  }
  if (is.null(x) || !length(x)) stop("RNA integer count layer is unavailable.")
  if (!inherits(x, "dgCMatrix")) x <- as(x, "dgCMatrix")
  if (length(x@x) && any(abs(x@x - round(x@x)) > 1e-8)) stop("Input counts are not integer-valued.")
  x
}
make_thresholds <- function(md, mt_available) {
  median_mad_lower <- function(x) max(200, floor(median(x, na.rm = TRUE) - 3 * mad(x, constant = 1.4826, na.rm = TRUE)))
  max_mt <- if (!mt_available) Inf else min(25, max(10, ceiling(median(md$percent.mt, na.rm = TRUE) + 3 * mad(md$percent.mt, constant = 1.4826, na.rm = TRUE))))
  c(min_nCount_RNA = median_mad_lower(md$nCount_RNA), min_nFeature_RNA = median_mad_lower(md$nFeature_RNA), max_percent_mt = max_mt)
}
scale_by_mt <- function(obj, features = Seurat::VariableFeatures(obj)) {
  regress <- if (stats::sd(obj$percent.mt, na.rm = TRUE) > 0) "percent.mt" else NULL
  Seurat::ScaleData(obj, features = features, vars.to.regress = regress, verbose = FALSE)
}
run_df_partition <- function(counts, md, partition_name) {
  if (ncol(counts) < 300L) stop("DoubletFinder partition ", partition_name, " has fewer than 300 cells.")
  # The globally stored count matrix is a union of study-specific feature
  # axes.  All-zero rows for this calling partition carry no information and
  # are removed before DoubletFinder is fit.
  counts <- counts[Matrix::rowSums(counts) > 0, , drop = FALSE]
  if (nrow(counts) < 500L) stop("DoubletFinder partition ", partition_name, " has fewer than 500 detected features.")
  message("DoubletFinder partition: ", partition_name, " (", ncol(counts), " QC-passed cells)")
  x <- Seurat::CreateSeuratObject(counts = counts, meta.data = md, project = partition_name)
  x <- Seurat::NormalizeData(x, verbose = FALSE)
  x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = min(3000L, nrow(x) - 1L), verbose = FALSE)
  x <- scale_by_mt(x)
  x <- Seurat::RunPCA(x, features = Seurat::VariableFeatures(x), npcs = min(30L, length(Seurat::VariableFeatures(x)) - 1L, ncol(x) - 1L), seed.use = 20260903, verbose = FALSE)
  dims <- seq_len(min(30L, ncol(Seurat::Embeddings(x, "pca"))))
  x <- Seurat::FindNeighbors(x, dims = dims, verbose = FALSE)
  x <- Seurat::FindClusters(x, resolution = 0.4, random.seed = 20260903, verbose = FALSE)
  x$DF_cluster <- as.character(x$seurat_clusters)
  tuning_cells <- colnames(x)
  if (length(tuning_cells) > 3000L) { set.seed(20260903); tuning_cells <- sample(tuning_cells, 3000L) }
  tuning_obj <- subset(x, cells = tuning_cells)
  sweep <- DoubletFinder::paramSweep(tuning_obj, PCs = dims, sct = FALSE)
  stats <- DoubletFinder::summarizeSweep(sweep, GT = FALSE)
  pk_table <- DoubletFinder::find.pK(stats)
  pks <- suppressWarnings(as.numeric(as.character(pk_table$pK)))
  best_pk <- pks[which.max(pk_table$BCmetric)]
  if (!is.finite(best_pk)) stop("No numeric pK found for ", partition_name)
  rate <- min(0.20, 0.008 * ncol(x) / 1000)
  n_raw <- max(1L, round(rate * ncol(x)))
  homo <- DoubletFinder::modelHomotypic(x$DF_cluster)
  n_adjusted <- max(1L, round(n_raw * (1 - homo)))
  metadata_columns_before_df <- colnames(x[[]])
  x <- DoubletFinder::doubletFinder(x, PCs = dims, pN = 0.25, pK = best_pk, nExp = n_adjusted, reuse.pANN = NULL, sct = FALSE)
  out <- x[[]]
  new_df_columns <- setdiff(colnames(out), metadata_columns_before_df)
  score_col <- grep("^pANN_", new_df_columns, value = TRUE)
  class_col <- grep("^DF.classifications_", new_df_columns, value = TRUE)
  if (length(score_col) != 1L || length(class_col) != 1L) stop("Cannot identify newly created DoubletFinder columns for ", partition_name)
  calls <- data.frame(Cell = rownames(out), DoubletFinder_score = as.numeric(out[[score_col]]),
    DoubletFinder_call = ifelse(grepl("Doublet", as.character(out[[class_col]])), "Doublet", "Singlet"),
    DF_partition = partition_name, stringsAsFactors = FALSE)
  params <- data.frame(DF_partition = partition_name, n_cells = ncol(x), expected_doublet_rate = rate,
    pN = 0.25, pK = best_pk, nExp_raw = n_raw, nExp_homotypic_adjusted = n_adjusted,
    homotypic_fraction = homo, pK_tuning_n_cells = length(tuning_cells),
    pK_tuned_by_subsample = ncol(x) > length(tuning_cells), stringsAsFactors = FALSE)
  rm(x, tuning_obj, sweep, stats, pk_table); invisible(gc())
  list(calls = calls, parameters = params)
}

# Select only features actually observed in every study.  This prevents the
# 50k raw-feature union from being interpreted as one common expression space
# during cross-study PCA/Harmony.  The full union remains untouched in RNA
# counts for QC and traceability.
select_shared_integration_features <- function(counts, md, nfeatures = 3000L) {
  study_cells <- split(colnames(counts), as.character(md[colnames(counts), "Study"]))
  study_objects <- lapply(names(study_cells), function(study) {
    cells <- study_cells[[study]]
    study_counts <- counts[, cells, drop = FALSE]
    study_counts <- study_counts[Matrix::rowSums(study_counts) > 0, , drop = FALSE]
    x <- Seurat::CreateSeuratObject(counts = study_counts, meta.data = md[cells, , drop = FALSE], project = study)
    x <- Seurat::NormalizeData(x, verbose = FALSE)
    x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = min(3000L, nrow(x) - 1L), verbose = FALSE)
    x
  })
  names(study_objects) <- names(study_cells)
  shared_features <- Reduce(intersect, lapply(study_objects, rownames))
  if (length(shared_features) < 1000L) stop("Fewer than 1,000 genes are detected in every selected study after QC.")
  candidate_features <- Seurat::SelectIntegrationFeatures(
    object.list = study_objects,
    nfeatures = min(as.integer(nfeatures), length(shared_features))
  )
  features <- intersect(candidate_features, shared_features)
  if (length(features) < 1000L) {
    # Fall back only to the cross-study shared variable features; never use a
    # study-private feature merely to reach a target panel size.
    features <- Reduce(intersect, lapply(study_objects, Seurat::VariableFeatures))
    features <- intersect(features, shared_features)
  }
  if (length(features) < 1000L) stop("Fewer than 1,000 shared variable features are available for integration.")
  summary <- data.frame(
    Study = names(study_objects),
    n_cells = vapply(study_objects, ncol, numeric(1)),
    n_detected_features = vapply(study_objects, nrow, numeric(1)),
    stringsAsFactors = FALSE
  )
  list(features = features, shared_features = shared_features, study_summary = summary)
}

if (!file.exists(selection_file)) stop("Stage-1 Microglia selection file is missing: ", selection_file)
selection <- read.csv(selection_file, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("Study", "analysis_cluster", "select_as_Microglia")
if (!all(required %in% names(selection))) stop("Selection CSV must contain: ", paste(required, collapse = ", "))
selection$selected <- tolower(trimws(as.character(selection$select_as_Microglia))) %in% c("true", "t", "1", "yes", "y")
if (!any(selection$selected)) stop("No Microglia cluster has been selected. Populate select_as_Microglia then rerun.")

selected_rows <- selection[selection$selected, c("Study", "analysis_cluster"), drop = FALSE]
stage1_checkpoints <- file.path(stage1_root, "object", paste0(safe_name(selected_rows$Study), "_full_count_processed.qs"))
if (any(!file.exists(stage1_checkpoints))) stop("Missing stage-1 checkpoint(s): ", paste(stage1_checkpoints[!file.exists(stage1_checkpoints)], collapse = "; "))

selected_objects <- list()
for (study in unique(selected_rows$Study)) {
  object_file <- file.path(stage1_root, "object", paste0(safe_name(study), "_full_count_processed.qs"))
  stage1_obj <- qs::qread(object_file)
  selected_clusters <- as.character(selected_rows$analysis_cluster[selected_rows$Study == study])
  cells <- colnames(stage1_obj)[as.character(stage1_obj$analysis_cluster) %in% selected_clusters]
  if (!length(cells)) stop("No cells matched selected clusters for ", study)
  counts <- get_counts(stage1_obj)[, cells, drop = FALSE]
  md <- stage1_obj[[]][cells, , drop = FALSE]
  author_df_columns <- grep("^(pANN_|DF.classifications_)", colnames(md), value = TRUE)
  if (length(author_df_columns)) {
    colnames(md)[match(author_df_columns, colnames(md))] <- make.unique(paste0("author_", author_df_columns))
  }
  md$stage1_analysis_cluster <- as.character(md$analysis_cluster)
  md$integration_batch <- paste(safe_name(md$Study), md$batch_id, sep = "__")
  # The Seurat merge S3 method unions features across studies, unlike direct matrix cbind.
  # (the six source matrices do not have identical row sets).
  selected_objects[[study]] <- Seurat::CreateSeuratObject(counts = counts, meta.data = md, project = study)
  rm(stage1_obj, counts, md); invisible(gc())
}
selected_input <- Reduce(function(x, y) merge(x, y), selected_objects)
# In Seurat v5, merging Assay5 objects leaves one counts.* layer per input.
# The default JoinLayers call folds each layer type (counts/data/scale.data)
# across studies, yielding a single complete "counts" layer for repeated QC.
selected_input <- SeuratObject::JoinLayers(selected_input, assay = "RNA")
all_counts <- get_counts(selected_input)
all_metadata <- selected_input[[]]
rm(selected_objects, selected_input); invisible(gc())

# Repeat QC on the selected Microglia counts. Thresholds are recomputed per
# original technical batch; percent.ribo remains a reviewed metric only.
mt_genes <- grep("^MT-", rownames(all_counts), value = TRUE)
ribo_genes <- grep("^RP[SL]", rownames(all_counts), value = TRUE)
all_metadata$nCount_RNA <- Matrix::colSums(all_counts)
all_metadata$nFeature_RNA <- Matrix::colSums(all_counts > 0)
all_metadata$percent.mt <- if (length(mt_genes)) Matrix::colSums(all_counts[mt_genes, , drop = FALSE]) / pmax(all_metadata$nCount_RNA, 1) * 100 else 0
all_metadata$percent.ribo <- if (length(ribo_genes)) Matrix::colSums(all_counts[ribo_genes, , drop = FALSE]) / pmax(all_metadata$nCount_RNA, 1) * 100 else 0
thresholds <- do.call(rbind, lapply(split(rownames(all_metadata), all_metadata$integration_batch), function(cells) {
  z <- make_thresholds(all_metadata[cells, , drop = FALSE], length(mt_genes) > 0)
  data.frame(integration_batch = all_metadata[cells[[1]], "integration_batch"], n_before = length(cells),
    min_nCount_RNA = z[["min_nCount_RNA"]], min_nFeature_RNA = z[["min_nFeature_RNA"]], max_percent_mt = z[["max_percent_mt"]])
}))
idx <- match(all_metadata$integration_batch, thresholds$integration_batch)
all_metadata$passes_basic_QC <- all_metadata$nCount_RNA >= thresholds$min_nCount_RNA[idx] &
  all_metadata$nFeature_RNA >= thresholds$min_nFeature_RNA[idx] & all_metadata$percent.mt <= thresholds$max_percent_mt[idx]
thresholds$n_after_basic_QC <- vapply(thresholds$integration_batch, function(x) sum(all_metadata$passes_basic_QC & all_metadata$integration_batch == x), numeric(1))
thresholds$n_removed_basic_QC <- thresholds$n_before - thresholds$n_after_basic_QC
if (sum(all_metadata$passes_basic_QC) < 1000L) stop("Fewer than 1,000 selected Microglia passed repeated QC.")

for (metric in c("nCount_RNA", "nFeature_RNA", "percent.mt")) {
  threshold_col <- c(nCount_RNA = "min_nCount_RNA", nFeature_RNA = "min_nFeature_RNA", percent.mt = "max_percent_mt")[[metric]]
  threshold_data <- thresholds[is.finite(thresholds[[threshold_col]]), c("integration_batch", threshold_col), drop = FALSE]
  p <- ggplot(all_metadata, aes(x = .data[[metric]])) + geom_histogram(bins = 70, fill = "#9ecae1", colour = "white", linewidth = 0.1) +
    facet_wrap(~Study, scales = "free", ncol = 3) + geom_vline(data = threshold_data, aes(xintercept = .data[[threshold_col]]), colour = "#d7301f", linewidth = 0.35) +
    theme_bw(base_size = 9) + labs(title = paste0("Selected Microglia: ", metric, " before repeated QC"), x = metric, y = "Cells")
  save_plot_both(p, file.path(figure_root, paste0("QC_", safe_name(metric))), 14, 8); print(p)
}

qc_cells <- rownames(all_metadata)[all_metadata$passes_basic_QC]
qc_counts <- all_counts[, qc_cells, drop = FALSE]
qc_md <- all_metadata[qc_cells, , drop = FALSE]
sizes <- table(qc_md$integration_batch)
qc_md$DF_partition_for_run <- qc_md$integration_batch
# DoubletFinder is run in technical batches where possible.  Small batches are
# pooled only within their own study, never across studies with different gene
# annotations or biological composition.  If a study's small-batch pool is
# still too small, the whole study is used for the call while the original
# integration_batch remains unchanged for all subsequent analysis.
for (study in unique(as.character(qc_md$Study))) {
  study_cells <- rownames(qc_md)[as.character(qc_md$Study) == study]
  study_batches <- unique(as.character(qc_md[study_cells, "integration_batch"]))
  study_small_batches <- study_batches[sizes[study_batches] < 300L]
  if (!length(study_small_batches)) next
  pooled_name <- paste0("pooled_small_batches__", safe_name(study))
  small_cells <- study_cells[as.character(qc_md[study_cells, "integration_batch"]) %in% study_small_batches]
  qc_md[small_cells, "DF_partition_for_run"] <- pooled_name
  if (length(small_cells) < 300L) {
    qc_md[study_cells, "DF_partition_for_run"] <- paste0("pooled_all_batches__", safe_name(study))
  }
}
df_partition_mapping <- unique(qc_md[, c("integration_batch", "DF_partition_for_run"), drop = FALSE])
rownames(df_partition_mapping) <- NULL
partitions <- split(qc_cells, qc_md$DF_partition_for_run)
df_checkpoint_root <- file.path(runtime_root, "DoubletFinder_partitions_shared_feature_v1")
dir.create(df_checkpoint_root, recursive = TRUE, showWarnings = FALSE)
run_or_load_df_partition <- function(partition_name) {
  stem <- file.path(df_checkpoint_root, safe_name(partition_name))
  calls_file <- paste0(stem, "_calls.csv")
  params_file <- paste0(stem, "_parameters.csv")
  if (file.exists(calls_file) && file.exists(params_file)) {
    calls <- read.csv(calls_file, stringsAsFactors = FALSE, check.names = FALSE)
    params <- read.csv(params_file, stringsAsFactors = FALSE, check.names = FALSE)
    return(list(calls = calls, parameters = params))
  }
  out <- run_df_partition(
    qc_counts[, partitions[[partition_name]], drop = FALSE],
    qc_md[partitions[[partition_name]], , drop = FALSE],
    partition_name
  )
  utils::write.csv(out$calls, calls_file, row.names = FALSE)
  utils::write.csv(out$parameters, params_file, row.names = FALSE)
  out
}
df_out <- lapply(names(partitions), run_or_load_df_partition)
df_calls <- do.call(rbind, lapply(df_out, `[[`, "calls")); rownames(df_calls) <- df_calls$Cell
df_params <- do.call(rbind, lapply(df_out, `[[`, "parameters"))
qc_md$DoubletFinder_score <- df_calls[rownames(qc_md), "DoubletFinder_score"]
qc_md$DoubletFinder_call <- factor(df_calls[rownames(qc_md), "DoubletFinder_call"], levels = c("Singlet", "Doublet"))
qc_md$DF_partition <- df_calls[rownames(qc_md), "DF_partition"]

# A pre-removal UMAP must be based on genes that are comparable across all six
# studies.  It remains an uncorrected diagnostic view; batch correction is
# applied only after DoubletFinder singlets have been retained below.
diagnostic_panel <- select_shared_integration_features(qc_counts, qc_md)
diagnostic <- Seurat::CreateSeuratObject(counts = qc_counts, meta.data = qc_md, project = "integrated_selected_Microglia")
diagnostic <- Seurat::NormalizeData(diagnostic, verbose = FALSE)
Seurat::VariableFeatures(diagnostic) <- diagnostic_panel$features
diagnostic <- scale_by_mt(diagnostic, features = diagnostic_panel$features)
diagnostic <- Seurat::RunPCA(diagnostic, features = diagnostic_panel$features, npcs = min(30L, length(diagnostic_panel$features) - 1L), seed.use = 20260903, verbose = FALSE)
dims <- seq_len(min(30L, ncol(Seurat::Embeddings(diagnostic, "pca"))))
diagnostic <- Seurat::FindNeighbors(diagnostic, dims = dims, verbose = FALSE)
diagnostic <- Seurat::FindClusters(diagnostic, resolution = 0.4, random.seed = 20260903, verbose = FALSE)
diagnostic$diagnostic_cluster <- diagnostic$seurat_clusters
diagnostic <- Seurat::RunUMAP(diagnostic, dims = dims, reduction.name = "umap_diagnostic", seed.use = 20260903, verbose = FALSE)
qs::qsave(diagnostic, file.path(object_root, "selected_Microglia_QC_before_doublet_removal.qs"), preset = "high")
p_diag_cluster <- Seurat::DimPlot(diagnostic, reduction = "umap_diagnostic", group.by = "diagnostic_cluster", label = TRUE, repel = TRUE) + ggtitle("Selected Microglia: QC-passed diagnostic UMAP")
p_diag_df <- Seurat::DimPlot(diagnostic, reduction = "umap_diagnostic", group.by = "DoubletFinder_call", cols = c(Singlet = "#bdbdbd", Doublet = "#d7301f")) + ggtitle("Selected Microglia: DoubletFinder before removal")
diagnostic_plots <- list(
  Diagnostic_UMAP_clusters_pre_doublet = p_diag_cluster,
  Diagnostic_UMAP_DoubletFinder_pre_doublet = p_diag_df
)
for (plot_name in names(diagnostic_plots)) {
  plot_object <- diagnostic_plots[[plot_name]]
  save_plot_both(plot_object, file.path(figure_root, plot_name), 10, 7)
  print(plot_object)
}

singlets <- rownames(qc_md)[qc_md$DoubletFinder_call == "Singlet"]
if (length(singlets) < 1000L) stop("Fewer than 1,000 Microglia singlets remain after DoubletFinder.")
integration_panel <- select_shared_integration_features(qc_counts[, singlets, drop = FALSE], qc_md[singlets, , drop = FALSE])
mg <- Seurat::CreateSeuratObject(counts = qc_counts[, singlets, drop = FALSE], meta.data = qc_md[singlets, , drop = FALSE], project = "integrated_Microglia")
mg <- Seurat::NormalizeData(mg, verbose = FALSE)
Seurat::VariableFeatures(mg) <- integration_panel$features
mg <- scale_by_mt(mg, features = integration_panel$features)
mg <- Seurat::RunPCA(mg, features = integration_panel$features, npcs = min(30L, length(integration_panel$features) - 1L), seed.use = 20260903, verbose = FALSE)
dims <- seq_len(min(30L, ncol(Seurat::Embeddings(mg, "pca"))))
mg <- harmony::RunHarmony(mg, group.by.vars = "integration_batch", reduction = "pca", dims.use = dims, reduction.save = "harmony", verbose = FALSE)
mg <- Seurat::FindNeighbors(mg, reduction = "harmony", dims = dims, verbose = FALSE)
mg <- Seurat::FindClusters(mg, resolution = 0.4, random.seed = 20260903, verbose = FALSE)
mg <- Seurat::RunUMAP(mg, reduction = "harmony", dims = dims, seed.use = 20260903, verbose = FALSE)
mg$integrated_cluster <- mg$seurat_clusters
markers <- Seurat::FindAllMarkers(mg, features = integration_panel$shared_features, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25, max.cells.per.ident = 3000, random.seed = 20260903, verbose = FALSE)
fc <- intersect(c("avg_log2FC", "avg_logFC"), names(markers))[[1]]
top_markers <- do.call(rbind, lapply(split(markers, markers$cluster), function(x) head(x[order(x$p_val_adj, -x[[fc]]), , drop = FALSE], 30)))
state_markers <- intersect(c("P2RY12", "TMEM119", "CX3CR1", "AIF1", "C1QA", "C1QB", "C1QC", "GPR34", "SALL1", "TREM2", "APOE", "LPL", "GPNMB", "SPP1", "LGALS3", "IFITM3", "HLA-DRA", "HLA-DRB1", "CTSD"), rownames(mg))
p_cluster <- Seurat::DimPlot(mg, group.by = "integrated_cluster", label = TRUE, repel = TRUE) + ggtitle("Integrated Microglia: Harmony UMAP")
p_batch <- Seurat::DimPlot(mg, group.by = "Study") + ggtitle("Integrated Microglia: study distribution")
p_dot <- Seurat::DotPlot(mg, features = state_markers, group.by = "integrated_cluster") + RotatedAxis() + ggtitle("Integrated Microglia: canonical and state marker dot plot")
p_feature <- Seurat::FeaturePlot(mg, features = state_markers, ncol = 4, min.cutoff = "q05", max.cutoff = "q95") + plot_annotation(title = "Integrated Microglia: marker expression")
final_plots <- list(
  Final_UMAP_clusters = p_cluster,
  Final_UMAP_study = p_batch,
  Microglia_state_marker_dotplot = p_dot,
  Microglia_state_marker_UMAP = p_feature
)
for (plot_name in names(final_plots)) {
  plot_object <- final_plots[[plot_name]]
  save_plot_both(plot_object, file.path(figure_root, plot_name), 13, 9)
  print(plot_object)
}

cluster_summary <- as.data.frame(table(integrated_cluster = mg$integrated_cluster), stringsAsFactors = FALSE); names(cluster_summary)[2] <- "n_cells"
mg0_template <- cluster_summary; mg0_template$select_as_MG0 <- NA; mg0_template$user_note <- NA_character_
utils::write.csv(thresholds, file.path(table_root, "selected_Microglia_QC_thresholds_and_retention.csv"), row.names = FALSE)
utils::write.csv(all_metadata, file.path(table_root, "selected_Microglia_all_cells_QC_metrics.csv"), row.names = TRUE)
utils::write.csv(df_partition_mapping, file.path(table_root, "selected_Microglia_DoubletFinder_partition_mapping.csv"), row.names = FALSE)
utils::write.csv(df_calls, file.path(table_root, "selected_Microglia_DoubletFinder_per_cell.csv"), row.names = FALSE)
utils::write.csv(df_params, file.path(table_root, "selected_Microglia_DoubletFinder_parameters.csv"), row.names = FALSE)
utils::write.csv(data.frame(Feature = diagnostic_panel$features, stringsAsFactors = FALSE), file.path(table_root, "diagnostic_UMAP_shared_integration_features.csv"), row.names = FALSE)
utils::write.csv(diagnostic_panel$study_summary, file.path(table_root, "diagnostic_UMAP_per_study_detected_feature_summary.csv"), row.names = FALSE)
utils::write.csv(data.frame(Feature = integration_panel$features, stringsAsFactors = FALSE), file.path(table_root, "final_integration_shared_variable_features.csv"), row.names = FALSE)
utils::write.csv(data.frame(Feature = integration_panel$shared_features, stringsAsFactors = FALSE), file.path(table_root, "final_integration_all_shared_detected_features.csv"), row.names = FALSE)
utils::write.csv(integration_panel$study_summary, file.path(table_root, "final_integration_per_study_detected_feature_summary.csv"), row.names = FALSE)
utils::write.csv(markers, file.path(table_root, "integrated_Microglia_all_cluster_markers.csv"), row.names = FALSE)
utils::write.csv(top_markers, file.path(table_root, "integrated_Microglia_top_markers.csv"), row.names = FALSE)
utils::write.csv(cluster_summary, file.path(table_root, "integrated_Microglia_cluster_summary.csv"), row.names = FALSE)
utils::write.csv(mg0_template, file.path(manual_root, "integrated_Microglia_MG0_cluster_selection_TEMPLATE.csv"), row.names = FALSE)
qs::qsave(mg, file.path(object_root, "integrated_Microglia_processed.qs"), preset = "high")
writeLines(capture.output(sessionInfo()), file.path(runtime_root, "integrated_Microglia_R_sessionInfo.txt"))
writeLines(c(
  "Input: stage-1 selected Microglia singlets, retaining their original integer RNA counts and provenance.",
  "Repeated QC: nCount_RNA, nFeature_RNA and percent.mt thresholds are recalculated by original technical batch; percent.ribo is retained for review and is not a hard filter.",
  "Repeated doublet removal: DoubletFinder is run within technical batches (small batches pooled only for the call); pK is tuned on a deterministic at-most-3,000-cell subset and then applied to the full batch. This detects residual doublets but has reduced sensitivity for homotypic Microglia doublets.",
  "Integration feature policy: the 50k RNA counts feature union is retained only as raw storage. PCA/Harmony and marker discovery use genes detected in every study; the variable subset is selected study-wise before the merged object is constructed.",
  "Integration: Harmony corrects the combined Study__batch variable after repeated doublet removal.",
  "Manual gate: review UMAP/dot plot/marker tables, then populate integrated_Microglia_MG0_cluster_selection_TEMPLATE.csv before the MG0-dependent diagnosis stage."
), file.path(runtime_root, "integrated_Microglia_method_manifest.txt"))
message("Stage 2 complete. Review figures and populate: ", file.path(manual_root, "integrated_Microglia_MG0_cluster_selection_TEMPLATE.csv"))
