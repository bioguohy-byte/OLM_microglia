#!/usr/bin/env Rscript

# Corrected historical MG0--MG5 diagnosis, stage 1.
#
# Scope correction: every diagnostic UMAP and every DoubletFinder/scDblFinder
# calculation starts from ALL QC-passed real cells in each individual study.
# Historical MG0--MG5 labels are only joined afterwards by Study + original
# barcode; they never restrict the diagnostic/calling input.

rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
invisible(gc())

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(DoubletFinder)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(BiocParallel)
  library(Matrix)
  library(qs)
})

set.seed(20260913)
project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
combine_root <- file.path(project_root, "combine")
source_object_root <- file.path(combine_root, "results/full_count_reanalysis/01_per_study/object")
history_file <- file.path(combine_root, "data/processed_data/mg_cca_meta_v2.qs")
mg00_file <- file.path(combine_root, "results/full_count_reanalysis/03_current_MG00_doublet_diagnosis/object/current_MG00_QC_before_doublet_removal_all_three_calls.qs")
result_root <- file.path(combine_root, "results/fullcell_historical_MG0_MG5_doublet_diagnosis_v2")
object_dir <- file.path(result_root, "object")
table_dir <- file.path(result_root, "table")
input_dir <- file.path(result_root, "scrublet_input")
runtime_dir <- file.path(result_root, "runtime")
for (directory in c(object_dir, table_dir, input_dir, runtime_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)

study_order <- c("Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
                 "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024")
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
get_counts <- function(object) {
  counts <- tryCatch(SeuratObject::LayerData(object[["RNA"]], layer = "counts"), error = function(e) NULL)
  if (is.null(counts) || !length(counts)) counts <- tryCatch(Seurat::GetAssayData(object, assay = "RNA", layer = "counts"), error = function(e) NULL)
  if (is.null(counts) || !length(counts)) counts <- tryCatch(Seurat::GetAssayData(object, assay = "RNA", slot = "counts"), error = function(e) NULL)
  if (is.null(counts) || !length(counts)) stop("RNA integer count layer is unavailable.")
  counts <- as(counts, "dgCMatrix")
  if (length(counts@x) && any(abs(counts@x - round(counts@x)) > 1e-8)) stop("RNA count layer contains non-integer values.")
  counts
}
write_sparse_bundle <- function(counts, folder, cell_ids) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  Matrix::writeMM(counts, file.path(folder, "counts_genes_by_cells.mtx"))
  write.table(data.frame(cell_id = cell_ids), file.path(folder, "cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(gene_id = rownames(counts)), file.path(folder, "genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
make_df_partition <- function(metadata) {
  if (!"batch_id" %in% names(metadata)) metadata$batch_id <- "single_batch"
  metadata$batch_id <- as.character(metadata$batch_id)
  metadata$batch_id[is.na(metadata$batch_id) | !nzchar(metadata$batch_id)] <- "missing_batch"
  sizes <- table(metadata$batch_id)
  metadata$DF_partition <- metadata$batch_id
  small <- names(sizes)[sizes < 300L]
  if (length(small)) {
    metadata$DF_partition[metadata$batch_id %in% small] <- "pooled_small_batches"
    if (sum(sizes[small]) < 300L) {
      largest <- names(sizes)[which.max(sizes)]
      metadata$DF_partition[metadata$batch_id %in% c(small, largest)] <- paste0("pooled_small_batches_with_", safe_name(largest))
    }
  }
  if (any(table(metadata$DF_partition) < 300L)) stop("A corrected DoubletFinder partition still has <300 QC-passed cells.")
  metadata
}
run_df_partition <- function(counts, metadata, partition_name) {
  object <- Seurat::CreateSeuratObject(counts = counts, meta.data = metadata, project = partition_name)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object, selection.method = "vst", nfeatures = min(3000L, nrow(object) - 1L), verbose = FALSE)
  features <- Seurat::VariableFeatures(object)
  object <- Seurat::ScaleData(object, features = features, vars.to.regress = "percent.mt", verbose = FALSE)
  npcs <- min(30L, length(features) - 1L, ncol(object) - 1L)
  object <- Seurat::RunPCA(object, features = features, npcs = npcs, seed.use = 20260913, verbose = FALSE)
  dims <- seq_len(min(30L, ncol(Seurat::Embeddings(object, "pca"))))
  object <- Seurat::FindNeighbors(object, dims = dims, verbose = FALSE)
  object <- Seurat::FindClusters(object, resolution = 0.4, random.seed = 20260913, verbose = FALSE)
  object$DF_cluster <- as.character(object$seurat_clusters)
  tuning_cells <- colnames(object)
  if (length(tuning_cells) > 3000L) { set.seed(20260913); tuning_cells <- sample(tuning_cells, 3000L) }
  tuning_object <- subset(object, cells = tuning_cells)
  sweep <- DoubletFinder::paramSweep(tuning_object, PCs = dims, sct = FALSE)
  sweep_summary <- DoubletFinder::summarizeSweep(sweep, GT = FALSE)
  pk_table <- DoubletFinder::find.pK(sweep_summary)
  pk_values <- suppressWarnings(as.numeric(as.character(pk_table$pK)))
  pK <- pk_values[which.max(pk_table$BCmetric)]
  if (!is.finite(pK)) stop("DoubletFinder pK optimisation returned no numeric pK for ", partition_name)
  rate <- min(0.20, 0.008 * ncol(object) / 1000)
  n_raw <- max(1L, round(rate * ncol(object)))
  homotypic <- DoubletFinder::modelHomotypic(object$DF_cluster)
  n_adjusted <- max(1L, round(n_raw * (1 - homotypic)))
  existing <- colnames(object[[]])
  object <- DoubletFinder::doubletFinder(object, PCs = dims, pN = 0.25, pK = pK, nExp = n_adjusted, reuse.pANN = NULL, sct = FALSE)
  added <- setdiff(colnames(object[[]]), existing)
  class_column <- grep("^DF.classifications_", added, value = TRUE)
  score_column <- grep("^pANN_", added, value = TRUE)
  if (length(class_column) != 1L || length(score_column) != 1L) stop("Could not identify new DoubletFinder output for ", partition_name)
  calls <- data.frame(Cell = colnames(object), DF_partition = partition_name,
    DoubletFinder_score = as.numeric(object[[score_column]][, 1]),
    DoubletFinder_call = ifelse(grepl("Doublet", as.character(object[[class_column]][, 1])), "Doublet", "Singlet"),
    stringsAsFactors = FALSE)
  parameters <- data.frame(DF_partition = partition_name, n_cells = ncol(object), expected_doublet_rate = rate,
    pN = 0.25, pK = pK, nExp_raw = n_raw, nExp_homotypic_adjusted = n_adjusted,
    homotypic_fraction = homotypic, pK_tuning_n_cells = length(tuning_cells), stringsAsFactors = FALSE)
  rm(object, tuning_object, sweep, sweep_summary, pk_table); invisible(gc())
  list(calls = calls, parameters = parameters)
}
run_sc_partition <- function(counts, metadata, partition_name) {
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(metadata)
  rate <- unique(metadata$expected_doublet_rate)
  if (length(rate) != 1L) stop("Expected one doublet-rate prior per partition.")
  sce <- scDblFinder::scDblFinder(sce, clusters = as.character(metadata$diagnostic_cluster), dbr = rate,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE), verbose = TRUE)
  data.frame(Cell = colnames(sce), DF_partition = partition_name,
    scDblFinder_score = as.numeric(SummarizedExperiment::colData(sce)$scDblFinder.score),
    scDblFinder_call = ifelse(tolower(as.character(SummarizedExperiment::colData(sce)$scDblFinder.class)) == "doublet", "Doublet", "Singlet"),
    scDblFinder_threshold = as.numeric(S4Vectors::metadata(sce)$scDblFinder.threshold), stringsAsFactors = FALSE)
}

if (any(!file.exists(c(history_file, mg00_file)))) stop("Historical MG0--MG5 map or MG00 map is missing.")
historical_object <- qs::qread(history_file)
historical_meta <- historical_object[[]]
if (length(setdiff(c("Study", "CellType1"), names(historical_meta)))) stop("Historical object lacks Study or CellType1.")
historical_map <- data.frame(Study = as.character(historical_meta$Study), original_barcode = rownames(historical_meta),
  historical_MG_state = as.character(historical_meta$CellType1), stringsAsFactors = FALSE)
historical_map <- historical_map[historical_map$historical_MG_state %in% paste0("MG", 0:5), , drop = FALSE]
historical_map$key <- paste(historical_map$Study, historical_map$original_barcode, sep = "||")
if (anyDuplicated(historical_map$key)) stop("Historical Study/barcode map is not unique.")
rm(historical_object, historical_meta)
mg00_object <- qs::qread(mg00_file)
mg00_meta <- mg00_object[[]]
if (length(setdiff(c("Study", "original_barcode", "MG00_status"), names(mg00_meta)))) stop("MG00 diagnostic map is incomplete.")
mg00_map <- data.frame(Study = as.character(mg00_meta$Study), original_barcode = as.character(mg00_meta$original_barcode),
  MG00 = grepl("MG00", as.character(mg00_meta$MG00_status)), stringsAsFactors = FALSE)
mg00_map$key <- paste(mg00_map$Study, mg00_map$original_barcode, sep = "||")
if (anyDuplicated(mg00_map$key)) stop("MG00 Study/barcode map is not unique.")
rm(mg00_object, mg00_meta); invisible(gc())

all_metadata <- list(); all_df_parameters <- list(); scrublet_manifest <- list(); mapping_summary <- list()
for (study in study_order) {
  safe <- safe_name(study)
  source_file <- file.path(source_object_root, paste0(safe, "_QC_before_doublet_removal.qs"))
  if (!file.exists(source_file)) stop("Missing all-QC-passed source object: ", source_file)
  message("\n===== corrected full-cell diagnosis: ", study, " =====")
  source <- qs::qread(source_file)
  counts <- get_counts(source)
  meta <- source[[]]
  if (!identical(colnames(counts), rownames(meta))) stop(study, ": count columns and metadata rows disagree.")
  if ("passes_basic_QC" %in% names(meta) && any(!meta$passes_basic_QC)) stop(study, ": source object contains non-QC-passed cells.")
  if (!"original_barcode" %in% names(meta)) stop(study, ": source object lacks original_barcode.")
  meta$Study <- study
  # Preserve the old reanalysis call strictly as provenance; all unprefixed
  # call columns written below are recalculated on this full-cell universe.
  for (field in intersect(c("DoubletFinder_score", "DoubletFinder_call", "DF_partition"), names(meta))) names(meta)[names(meta) == field] <- paste0("prior_reanalysis_", field)
  meta <- make_df_partition(meta)
  partition_cells <- split(rownames(meta), meta$DF_partition)
  partition_rates <- vapply(partition_cells, function(cells) min(0.20, 0.008 * length(cells) / 1000), numeric(1))
  meta$expected_doublet_rate <- unname(partition_rates[meta$DF_partition])
  key <- paste(meta$Study, meta$original_barcode, sep = "||")
  history_index <- match(key, historical_map$key)
  mg00_index <- match(key, mg00_map$key)
  meta$historical_MG_state <- factor(historical_map$historical_MG_state[history_index], levels = paste0("MG", 0:5))
  meta$historical_label_available <- !is.na(history_index)
  meta$historical_MG0 <- !is.na(meta$historical_MG_state) & as.character(meta$historical_MG_state) == "MG0"
  meta$MG00 <- !is.na(mg00_index) & mg00_map$MG00[mg00_index]
  # The diagnostic UMAP is fit from every QC-passed real cell in the study.
  diagnostic <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta, project = study)
  diagnostic <- Seurat::NormalizeData(diagnostic, verbose = FALSE)
  diagnostic <- Seurat::FindVariableFeatures(diagnostic, selection.method = "vst", nfeatures = min(3000L, nrow(diagnostic) - 1L), verbose = FALSE)
  features <- Seurat::VariableFeatures(diagnostic)
  diagnostic <- Seurat::ScaleData(diagnostic, features = features, vars.to.regress = "percent.mt", verbose = FALSE)
  npcs <- min(30L, length(features) - 1L, ncol(diagnostic) - 1L)
  diagnostic <- Seurat::RunPCA(diagnostic, features = features, npcs = npcs, seed.use = 20260913, verbose = FALSE)
  dims <- seq_len(min(30L, ncol(Seurat::Embeddings(diagnostic, "pca"))))
  diagnostic <- Seurat::FindNeighbors(diagnostic, dims = dims, verbose = FALSE)
  diagnostic <- Seurat::FindClusters(diagnostic, resolution = 0.4, random.seed = 20260913, verbose = FALSE)
  diagnostic$diagnostic_cluster <- as.character(diagnostic$seurat_clusters)
  diagnostic <- Seurat::RunUMAP(diagnostic, dims = dims, reduction.name = "umap_fullcell_diagnostic", seed.use = 20260913, verbose = FALSE)
  meta <- diagnostic[[]]
  df_results <- lapply(names(partition_cells), function(partition) run_df_partition(counts[, partition_cells[[partition]], drop = FALSE], meta[partition_cells[[partition]], , drop = FALSE], partition))
  df_calls <- do.call(rbind, lapply(df_results, `[[`, "calls"))
  rownames(df_calls) <- df_calls$Cell
  diagnostic$DoubletFinder_score <- df_calls[colnames(diagnostic), "DoubletFinder_score"]
  diagnostic$DoubletFinder_call <- factor(df_calls[colnames(diagnostic), "DoubletFinder_call"], levels = c("Singlet", "Doublet"))
  diagnostic$DF_partition <- df_calls[colnames(diagnostic), "DF_partition"]
  meta <- diagnostic[[]]
  sc_calls <- do.call(rbind, lapply(names(partition_cells), function(partition) run_sc_partition(counts[, partition_cells[[partition]], drop = FALSE], meta[partition_cells[[partition]], , drop = FALSE], partition)))
  rownames(sc_calls) <- sc_calls$Cell
  diagnostic$scDblFinder_score <- sc_calls[colnames(diagnostic), "scDblFinder_score"]
  diagnostic$scDblFinder_call <- factor(sc_calls[colnames(diagnostic), "scDblFinder_call"], levels = c("Singlet", "Doublet"))
  diagnostic$scDblFinder_threshold <- sc_calls[colnames(diagnostic), "scDblFinder_threshold"]
  qs::qsave(diagnostic, file.path(object_dir, paste0(safe, "_all_QC_passed_fullcell_partial_calls.qs")), preset = "high")
  # Scrublet receives exactly the same complete QC-passed partitioned cells.
  for (partition in names(partition_cells)) {
    cells <- partition_cells[[partition]]
    part_counts <- counts[, cells, drop = FALSE]
    part_counts <- part_counts[Matrix::rowSums(part_counts) > 0, , drop = FALSE]
    part_safe <- paste0(safe, "__", safe_name(partition))
    write_sparse_bundle(part_counts, file.path(input_dir, part_safe), cells)
    scrublet_manifest[[length(scrublet_manifest) + 1L]] <- data.frame(Study = study, safe_study = safe, DF_partition = partition,
      safe_partition = part_safe, n_cells = length(cells), n_genes = nrow(part_counts), expected_doublet_rate = partition_rates[[partition]], stringsAsFactors = FALSE)
  }
  current <- diagnostic[[]]
  current$Cell <- rownames(current)
  umap <- Seurat::Embeddings(diagnostic, "umap_fullcell_diagnostic")
  current$UMAP_1 <- umap[current$Cell, 1]
  current$UMAP_2 <- umap[current$Cell, 2]
  all_metadata[[study]] <- current
  df_parameter_table <- do.call(rbind, lapply(df_results, `[[`, "parameters")); df_parameter_table$Study <- study
  all_df_parameters[[study]] <- df_parameter_table
  mapping_summary[[study]] <- data.frame(Study = study, n_all_QC_passed = nrow(current),
    n_historical_MG0_MG5_mapped = sum(!is.na(current$historical_MG_state)), n_historical_MG0_mapped = sum(current$historical_MG0),
    n_MG00_mapped = sum(current$MG00), stringsAsFactors = FALSE)
  rm(source, counts, meta, diagnostic, df_results, df_calls, sc_calls, current, umap); invisible(gc())
}

per_cell <- do.call(rbind, all_metadata)
utils::write.csv(per_cell, file.path(table_dir, "per_cell_full_QC_passed_partial_calls_with_historical_MG0_MG5_mapping.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, all_df_parameters), file.path(table_dir, "DoubletFinder_parameters_full_QC_passed_cells.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, scrublet_manifest), file.path(input_dir, "scrublet_manifest.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, mapping_summary), file.path(table_dir, "historical_MG0_MG5_and_MG00_mapping_coverage.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "stage1_R_sessionInfo.txt"))
writeLines(c(
  "Scope correction: ALL QC-passed real cells within each study are the only input to the diagnostic UMAP, DoubletFinder, scDblFinder and Scrublet.",
  "Historical MG0--MG5 labels are mapped by Study + original_barcode after calls. They never filter or train a caller.",
  "DoubletFinder and scDblFinder run in within-study technical partitions; small batches are pooled only for calling when necessary.",
  "No cross-study cells are combined for calling or diagnostic UMAPs. No caller output removes or relabels any cell in this branch."
), file.path(runtime_dir, "scope_and_method_manifest.txt"))
message("Corrected full-cell stage 1 complete. Run 27.fullcell_historical_MG0_MG5_doublet_stage2_scrublet_24server.py, then stage 3.")
