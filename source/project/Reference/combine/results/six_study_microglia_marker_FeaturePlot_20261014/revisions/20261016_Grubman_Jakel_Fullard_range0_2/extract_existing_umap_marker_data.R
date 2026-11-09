#!/usr/bin/env Rscript

# Read the six existing full-cell objects and extract only stored coordinates,
# three RNA markers, and existing annotation fields for reproducible plotting.
suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(Matrix)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: extract_existing_umap_marker_data.R <output_dir>")
out_dir <- args[[1]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plot_data"), showWarnings = FALSE)
dir.create(file.path(out_dir, "audit"), showWarnings = FALSE)
project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
studies <- c("Su", "Chia", "Fullard", "Jakel", "Grubman", "Yao")
source_names <- c("Su_CSC_alldata_v2.qs", "Chia_v2.qs", "Fullard_v2.qs",
                  "Jakel_v2.qs", "Grubman_v2.qs", "Yao_v2.qs")
mg_labels <- list(Su = "MG", Chia = "Microglia", Fullard = "Mic",
                  Jakel = "Microglia_Macrophages", Grubman = "Microglia", Yao = "Micro")
inventory <- list()
expression_audit <- list()
label_audit <- list()
normalization_audit <- list()

get_single_layer <- function(object, layer) {
  assay <- object[["RNA"]]
  available <- Layers(assay, search = NA)
  if (!layer %in% available) {
    stop("Exact RNA layer '", layer, "' not present. Layers: ", paste(available, collapse = ", "))
  }
  GetAssayData(object, assay = "RNA", layer = layer)
}

for (i in seq_along(studies)) {
  study <- studies[[i]]
  source_path <- file.path(project_root, paste0(study, "_et_al."), "data/processed_data", source_names[[i]])
  if (study == "Su") {
    source_path <- file.path(project_root, "combine/results/celltype_audit_y20260918/checkpoints/Su_audit_harmony_clustered.qs")
  }
  message(format(Sys.time(), "%F %T"), " Reading ", study, ": ", source_path)
  obj <- qread(source_path, nthreads = 1)
  if (!"umap" %in% Reductions(obj)) stop(study, ": existing UMAP absent")
  if (!"RNA" %in% Assays(obj)) stop(study, ": RNA assay absent")
  if (!"CellType" %in% colnames(obj[[]])) stop(study, ": CellType field absent")
  cells <- colnames(obj)
  coordinates <- Embeddings(obj, "umap")
  if (ncol(coordinates) != 2L || !identical(rownames(coordinates), cells)) {
    stop(study, ": UMAP is not a barcode-aligned two-dimensional embedding")
  }
  if (any(!is.finite(coordinates))) stop(study, ": non-finite UMAP coordinates")
  rna_data <- get_single_layer(obj, "data")
  rna_counts <- get_single_layer(obj, "counts")
  if (!identical(colnames(rna_data), cells) || !identical(colnames(rna_counts), cells)) {
    stop(study, ": RNA columns do not match cell order")
  }
  features <- vapply(markers, function(gene) {
    hits <- rownames(rna_data)[toupper(rownames(rna_data)) == gene]
    if (length(hits) != 1L) stop(study, ": marker not uniquely resolved: ", gene)
    hits[[1]]
  }, character(1))
  if (!all(features %in% rownames(rna_counts))) stop(study, ": count marker absent")
  stored_expression <- rna_data[features, cells, drop = FALSE]
  counts <- rna_counts[features, cells, drop = FALSE]
  if (any(!is.finite(stored_expression@x)) || any(stored_expression@x < 0)) {
    stop(study, ": RNA/data contains negative or non-finite marker values")
  }
  library_sizes <- Matrix::colSums(rna_counts)
  if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) stop(study, ": non-positive RNA library size")
  # Derive the same LogNormalize scale for all six studies. Only three marker
  # rows are transformed; the full RNA count sum supplies the denominator.
  # Source UMAPs, source expression layers, QC, and cell labels are untouched.
  expression <- counts %*% Diagonal(x = 10000 / library_sizes)
  expression@x <- log1p(expression@x)
  dimnames(expression) <- dimnames(counts)
  stored_equals_counts <- isTRUE(all.equal(as.matrix(stored_expression), as.matrix(counts), tolerance = 1e-10, check.attributes = FALSE))
  stored_vs_log_maxabsdiff <- max(abs(as.matrix(stored_expression) - as.matrix(expression)))
  metadata_columns <- intersect(c("CellType", "cell.type", "celltype", "cellType", "Celltypes", "cell_type", "seurat_clusters"), colnames(obj[[]]))
  metadata <- obj[[]][cells, metadata_columns, drop = FALSE]
  normalized_commands <- names(obj@commands)[grepl("NormalizeData", names(obj@commands))]
  normalization_audit[[study]] <- data.frame(
    study = study,
    command = if (length(normalized_commands)) normalized_commands else "<not recorded>",
    parameters = if (length(normalized_commands)) vapply(normalized_commands, function(nm) {
      paste(capture.output(str(obj@commands[[nm]]@params)), collapse = " ")
    }, character(1)) else "Source NormalizeData command not recorded; plotting marker values derived from RNA/counts",
    stringsAsFactors = FALSE
  )
  for (gene in markers) {
    vals <- as.numeric(expression[features[[gene]], ])
    count_vals <- as.numeric(counts[features[[gene]], ])
    positive <- vals[vals > 0]
    expression_audit[[paste(study, gene, sep = "__")]] <- data.frame(
      study = study, marker = gene, actual_feature = features[[gene]], n_cells = length(vals),
      n_positive_data = sum(vals > 0), percent_positive_data = 100 * mean(vals > 0),
      n_positive_counts = sum(count_vals > 0), data_min = min(vals), data_max = max(vals),
      positive_q95 = if (length(positive)) unname(quantile(positive, 0.95)) else 0,
      positive_q99 = if (length(positive)) unname(quantile(positive, 0.99)) else 0,
      positive_q995 = if (length(positive)) unname(quantile(positive, 0.995)) else 0,
      positive_mean = if (length(positive)) mean(positive) else 0,
      stringsAsFactors = FALSE
    )
  }
  for (column in setdiff(metadata_columns, "seurat_clusters")) {
    values <- as.character(metadata[[column]])
    values[is.na(values) | !nzchar(values)] <- "<missing>"
    for (label in sort(unique(values))) {
      idx <- which(values == label)
      for (gene in markers) {
        vals <- as.numeric(expression[features[[gene]], idx])
        label_audit[[paste(study, column, label, gene, sep = "__")]] <- data.frame(
          study = study, annotation_column = column, celltype = label, marker = gene,
          n_cells = length(idx), percent_expressed = 100 * mean(vals > 0),
          mean_log_normalized_expression = mean(vals),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  source_info <- file.info(source_path)
  inventory[[study]] <- data.frame(
    study = study, source_path = source_path, source_size_bytes = source_info$size,
    source_mtime = as.character(source_info$mtime), n_cells = ncol(obj), n_rna_features = nrow(rna_data),
    stored_reduction = "umap", annotation_column = "CellType",
    highlighted_existing_label = mg_labels[[study]],
    n_highlighted_existing_label = sum(as.character(metadata$CellType) %in% mg_labels[[study]]),
    x_min = min(coordinates[, 1]), x_max = max(coordinates[, 1]),
    y_min = min(coordinates[, 2]), y_max = max(coordinates[, 2]),
    barcode_order_verified = TRUE, umap_recomputed = FALSE, source_normalization_modified = FALSE,
    plot_normalization = "log1p(marker RNA count / full-cell RNA total * 10000)",
    source_three_marker_data_equals_counts = stored_equals_counts,
    source_three_marker_data_vs_plot_log_maxabsdiff = stored_vs_log_maxabsdiff,
    stringsAsFactors = FALSE
  )
  bundle <- list(study = study, source_path = source_path, cells = cells,
    umap = coordinates, expression = expression, counts = counts,
    stored_expression = stored_expression, library_sizes = library_sizes,
    features = features, metadata = metadata, mg_labels = mg_labels[[study]])
  bundle_path <- file.path(out_dir, "plot_data", paste0(study, "_stored_umap_three_RNA_markers.rds"))
  saveRDS(bundle, bundle_path, compress = "gzip")
  check <- readRDS(bundle_path)
  stopifnot(identical(check$umap, coordinates), identical(check$expression, expression), identical(check$cells, cells))
  message(study, ": ", ncol(obj), " full cells; all 3 markers present; stored UMAP preserved exactly")
  rm(obj, rna_data, rna_counts, expression, counts, stored_expression, library_sizes, metadata, bundle, check)
  gc(verbose = FALSE)
  write.table(do.call(rbind, inventory), file.path(out_dir, "audit/source_object_inventory.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
  write.table(do.call(rbind, expression_audit), file.path(out_dir, "audit/marker_expression_ranges.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
  write.table(do.call(rbind, label_audit), file.path(out_dir, "audit/marker_expression_by_existing_celltype.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
  write.table(do.call(rbind, normalization_audit), file.path(out_dir, "audit/source_normalization_commands.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "audit/extraction_sessionInfo.txt"))
message("EXTRACTION_COMPLETE")
