#!/usr/bin/env Rscript

# Artificial-doublet visualisation, part 1 of 3.
#
# This script intentionally does NOT run a new doublet call.  It exports the
# three algorithms' method-matched artificial profiles from the exact current
# MG00 diagnostic cell universe (before any diagnostic doublet is removed).
# Existing DoubletFinder/scDblFinder/Scrublet calls remain the calls displayed
# later for the real cells.  All artificial profiles are newly and
# deterministically constructed solely for this visualisation branch.

rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
invisible(gc())

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(qs)
  library(scDblFinder)
})

set.seed(20260908)

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
combine_root <- file.path(project_root, "combine")
input_file <- file.path(
  combine_root,
  "results/full_count_reanalysis/03_current_MG00_doublet_diagnosis/object/current_MG00_QC_before_doublet_removal_all_three_calls.qs"
)
historical_table_file <- file.path(
  combine_root,
  "results/MG0_MG5_doublet_diagnosis/table/six_study_cell_metadata_with_all_doublet_calls.csv"
)
result_root <- file.path(combine_root, "results/artificial_doublet_visualization")
input_dir <- file.path(result_root, "input")
table_dir <- file.path(result_root, "table")
runtime_dir <- file.path(result_root, "runtime")
for (directory in c(input_dir, table_dir, runtime_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

study_order <- c(
  "Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
  "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"
)
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
read_tsv1 <- function(path) read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
get_counts <- function(object) {
  counts <- tryCatch(SeuratObject::LayerData(object[["RNA"]], layer = "counts"), error = function(e) NULL)
  if (is.null(counts) || !length(counts)) {
    counts <- tryCatch(Seurat::GetAssayData(object, assay = "RNA", layer = "counts"), error = function(e) NULL)
  }
  if (is.null(counts) || !length(counts)) {
    counts <- tryCatch(Seurat::GetAssayData(object, assay = "RNA", slot = "counts"), error = function(e) NULL)
  }
  if (is.null(counts) || !length(counts)) stop("The input object has no RNA count layer.")
  counts <- as(counts, "dgCMatrix")
  if (length(counts@x) && any(abs(counts@x - round(counts@x)) > 1e-8)) {
    stop("The current MG00 diagnostic RNA layer is not integer-valued.")
  }
  counts
}
write_sparse_bundle <- function(counts, folder, matrix_file, cell_file, gene_file, cell_ids) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  counts <- as(counts, "dgCMatrix")
  Matrix::writeMM(counts, file.path(folder, matrix_file))
  write.table(data.frame(cell_id = cell_ids), file.path(folder, cell_file),
    sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(gene_id = rownames(counts)), file.path(folder, gene_file),
    sep = "\t", quote = FALSE, row.names = FALSE)
}

if (any(!file.exists(c(input_file, historical_table_file)))) {
  stop("The current MG00 object or historical label table is missing.")
}

message("Loading the current pre-removal MG00 diagnostic object ...")
current <- qs::qread(input_file)
counts_all <- get_counts(current)
meta <- current[[]]
meta$Cell <- rownames(meta)
required <- c(
  "Study", "original_barcode", "DF_partition", "diagnostic_cluster",
  "MG00_status", "DoubletFinder_call", "scDblFinder_call", "Scrublet_call"
)
if (length(setdiff(required, colnames(meta)))) {
  stop("Current MG00 metadata lacks: ", paste(setdiff(required, colnames(meta)), collapse = ", "))
}
if (!setequal(study_order, unique(as.character(meta$Study)))) {
  stop("The current object does not contain exactly the expected six studies.")
}
if (!identical(colnames(counts_all), meta$Cell)) stop("Count columns and metadata rows are misaligned.")

# Link original MG0 only by study plus the original barcode.  A missing old
# label is retained as unavailable, never recoded as a non-MG0 cell.
historical <- read.csv(historical_table_file, stringsAsFactors = FALSE, check.names = FALSE)
if (length(setdiff(c("Cell", "Study", "CellType1"), colnames(historical)))) {
  stop("The historical metadata table does not contain Cell, Study and CellType1.")
}
historical$key <- paste(historical$Study, historical$Cell, sep = "||")
meta$key <- paste(meta$Study, meta$original_barcode, sep = "||")
if (anyDuplicated(historical$key) || anyDuplicated(meta$key)) stop("A study/barcode mapping key is duplicated.")
old_index <- match(meta$key, historical$key)
meta$historical_label_available <- !is.na(old_index)
meta$historical_MG0 <- meta$historical_label_available & historical$CellType1[old_index] == "MG0"

# pN=0.25 is the exact setting used in the formal current DoubletFinder calls.
# The expression below is the number used by DoubletFinder::doubletFinder().
df_pN <- 0.25
partition_cells <- split(meta$Cell, as.character(meta$DF_partition))
partition_names <- sort(names(partition_cells))
if (!length(partition_names) || any(lengths(partition_cells) < 300L)) {
  stop("Each saved DF_partition must contain at least 300 QC-passed cells.")
}

manifest_rows <- list()
df_parent_rows <- list()
sc_origin_rows <- list()
for (i in seq_along(partition_names)) {
  partition_name <- partition_names[[i]]
  cells <- partition_cells[[partition_name]]
  partition_meta <- meta[cells, , drop = FALSE]
  partition_studies <- unique(as.character(partition_meta$Study))
  if (length(partition_studies) != 1L) {
    stop("DF_partition ", partition_name, " crosses studies; do not create artificial profiles across studies.")
  }
  study <- partition_studies[[1]]
  part_safe <- safe_name(partition_name)
  part_dir <- file.path(input_dir, part_safe)
  part_counts <- counts_all[, cells, drop = FALSE]
  detected_rows <- Matrix::rowSums(part_counts) > 0
  if (sum(detected_rows) < 500L) stop(partition_name, " has fewer than 500 detected features.")
  part_counts_detected <- part_counts[detected_rows, , drop = FALSE]

  # DoubletFinder: its public implementation samples with replacement and
  # constructs the arithmetic mean of the two raw-count columns.
  set.seed(2026090800 + i)
  n_df <- round(ncol(part_counts) / (1 - df_pN) - ncol(part_counts))
  parent_1 <- sample(cells, n_df, replace = TRUE)
  parent_2 <- sample(cells, n_df, replace = TRUE)
  df_counts <- (part_counts[, parent_1, drop = FALSE] + part_counts[, parent_2, drop = FALSE]) / 2
  df_cell_ids <- paste0("DFart__", part_safe, "__", seq_len(n_df))
  colnames(df_counts) <- df_cell_ids
  write_sparse_bundle(df_counts, part_dir,
    "DoubletFinder_artificial_counts_genes_by_cells.mtx",
    "DoubletFinder_artificial_cells.tsv", "DoubletFinder_artificial_genes.tsv", df_cell_ids)
  df_parents <- data.frame(
    artificial_cell_id = df_cell_ids, Study = study, DF_partition = partition_name,
    parent_1 = parent_1, parent_2 = parent_2,
    construction = "DoubletFinder pN=0.25: (raw parent-1 + raw parent-2)/2",
    stringsAsFactors = FALSE
  )
  write.csv(df_parents, file.path(part_dir, "DoubletFinder_artificial_parent_pairs.csv"), row.names = FALSE)
  df_parent_rows[[partition_name]] <- df_parents

  # scDblFinder: use its public helper with the same defaults implied by the
  # current call (propRandom=0).  It returns generated counts and origin
  # clusters.  The helper does not expose individual parent barcodes for all
  # of its random/metacell/triplet profiles, so no parent ID is invented.
  set.seed(2026091800 + i)
  clusters <- as.factor(as.character(partition_meta$diagnostic_cluster))
  n_sc <- min(25000L, max(1500L, ceiling(ncol(part_counts_detected) * 0.8), 10L * length(unique(clusters))^2))
  sc_output <- scDblFinder::getArtificialDoublets(
    part_counts_detected, n = n_sc, clusters = clusters, propRandom = 0
  )
  sc_counts <- as(sc_output$counts, "dgCMatrix")
  sc_cell_ids <- paste0("scDblart__", part_safe, "__", seq_len(ncol(sc_counts)))
  colnames(sc_counts) <- sc_cell_ids
  write_sparse_bundle(sc_counts, part_dir,
    "scDblFinder_artificial_counts_genes_by_cells.mtx",
    "scDblFinder_artificial_cells.tsv", "scDblFinder_artificial_genes.tsv", sc_cell_ids)
  sc_origins <- data.frame(
    artificial_cell_id = sc_cell_ids, Study = study, DF_partition = partition_name,
    scDblFinder_origin_cluster = as.character(sc_output$origins),
    construction = "scDblFinder::getArtificialDoublets (official helper; origin cluster retained)",
    stringsAsFactors = FALSE
  )
  write.csv(sc_origins, file.path(part_dir, "scDblFinder_artificial_origin_clusters.csv"), row.names = FALSE)
  sc_origin_rows[[partition_name]] <- sc_origins

  # Scrublet uses the integer count matrix.  Part 2 invokes its
  # simulate_doublets() method directly and records the exposed parent index.
  scrublet_dir <- file.path(part_dir, "scrublet_input")
  write_sparse_bundle(part_counts_detected, scrublet_dir,
    "counts_genes_by_cells.mtx", "cells.tsv", "genes.tsv", cells)

  manifest_rows[[partition_name]] <- data.frame(
    Study = study, DF_partition = partition_name, safe_partition = part_safe,
    n_real_cells = ncol(part_counts), n_detected_features = nrow(part_counts_detected),
    n_DoubletFinder_artificial = ncol(df_counts),
    n_scDblFinder_artificial = ncol(sc_counts),
    Scrublet_sim_doublet_ratio = 1 / 3,
    stringsAsFactors = FALSE
  )
  message("Prepared ", study, " / ", partition_name, ": ", ncol(part_counts),
    " real, ", ncol(df_counts), " DoubletFinder artificial, ", ncol(sc_counts), " scDblFinder artificial.")
  rm(part_counts, part_counts_detected, df_counts, sc_counts, sc_output)
  invisible(gc())
}

manifest <- do.call(rbind, manifest_rows)
manifest <- manifest[order(match(manifest$Study, study_order), manifest$DF_partition), , drop = FALSE]
write.csv(manifest, file.path(input_dir, "partition_manifest.csv"), row.names = FALSE)
write.csv(meta[, c(
  "Cell", "Study", "original_barcode", "DF_partition", "diagnostic_cluster",
  "historical_label_available", "historical_MG0", "MG00_status",
  "DoubletFinder_call", "scDblFinder_call", "Scrublet_call"
)], file.path(table_dir, "real_cell_metadata_and_existing_calls.csv"), row.names = FALSE)
write.csv(do.call(rbind, df_parent_rows), file.path(table_dir, "DoubletFinder_artificial_parent_pairs.csv"), row.names = FALSE)
write.csv(do.call(rbind, sc_origin_rows), file.path(table_dir, "scDblFinder_artificial_origin_clusters.csv"), row.names = FALSE)

mapping_summary <- aggregate(
  Cell ~ Study + historical_label_available + historical_MG0,
  data = meta[, c("Cell", "Study", "historical_label_available", "historical_MG0")], FUN = length
)
names(mapping_summary)[4] <- "n_real_cells"
write.csv(mapping_summary, file.path(table_dir, "historical_MG0_mapping_availability.csv"), row.names = FALSE)
writeLines(c(
  "Input: current MG00 QC-passed, pre-removal diagnostic object (16,931 real cells).",
  "No doublet call is recalculated in this branch. Real-cell calls remain the completed current DoubletFinder/scDblFinder/Scrublet calls.",
  "DoubletFinder artificial profiles: newly generated with pN=0.25, matching the project call parameter; parent IDs are recorded.",
  "scDblFinder artificial profiles: newly generated by scDblFinder::getArtificialDoublets using propRandom=0; official helper output contains origin clusters but does not expose universal individual parent IDs.",
  "Scrublet input: sparse integer counts are exported per saved technical DF_partition; the following visible Python stage calls Scrublet::simulate_doublets with ratio 1/3 and records Scrublet's exposed parent indices.",
  "Artificial profiles are never used to alter QC, labels, or doublet calls for real cells."
), file.path(runtime_dir, "stage1_method_manifest.txt"))
writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "stage1_R_sessionInfo.txt"))
message("Artificial-doublet stage 1 complete. Run stage 2 in the visible RStudio console/terminal, then stage 3.")
