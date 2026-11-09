#!/usr/bin/env Rscript

# Corrected artificial-doublet validation, stage 1.
# Parent universe: every QC-passed real cell in each study, split only by the
# recorded within-study technical calling partition. Historical MG0/MG00 are
# overlays, never parent-pool filters.

rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
invisible(gc())
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix); library(qs); library(scDblFinder) })
set.seed(20260913)

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
combine_root <- file.path(project_root, "combine")
diagnosis_root <- file.path(combine_root, "results/fullcell_historical_MG0_MG5_doublet_diagnosis_v2")
result_root <- file.path(combine_root, "results/fullcell_artificial_doublet_validation_v2")
input_dir <- file.path(result_root, "input")
table_dir <- file.path(result_root, "table")
runtime_dir <- file.path(result_root, "runtime")
for (directory in c(input_dir, table_dir, runtime_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
study_order <- c("Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021", "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024")
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
get_counts <- function(object) {
  counts <- tryCatch(SeuratObject::LayerData(object[["RNA"]], layer = "counts"), error = function(e) NULL)
  if (is.null(counts) || !length(counts)) counts <- tryCatch(Seurat::GetAssayData(object, assay = "RNA", layer = "counts"), error = function(e) NULL)
  if (is.null(counts) || !length(counts)) counts <- tryCatch(Seurat::GetAssayData(object, assay = "RNA", slot = "counts"), error = function(e) NULL)
  if (is.null(counts) || !length(counts)) stop("RNA integer counts are unavailable.")
  counts <- as(counts, "dgCMatrix")
  if (length(counts@x) && any(abs(counts@x - round(counts@x)) > 1e-8)) stop("Count matrix is not integer-valued.")
  counts
}
write_bundle <- function(counts, folder, prefix, cell_ids) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  Matrix::writeMM(as(counts, "dgCMatrix"), file.path(folder, paste0(prefix, "_counts_genes_by_cells.mtx")))
  write.table(data.frame(cell_id = cell_ids), file.path(folder, paste0(prefix, "_cells.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(gene_id = rownames(counts)), file.path(folder, paste0(prefix, "_genes.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

all_manifest <- list(); all_df_parents <- list(); all_sc_origins <- list(); all_real_metadata <- list()
for (study_i in seq_along(study_order)) {
  study <- study_order[[study_i]]; safe <- safe_name(study)
  source_file <- file.path(diagnosis_root, "object", paste0(safe, "_all_QC_passed_fullcell_all_three_calls.qs"))
  if (!file.exists(source_file)) stop("Missing corrected all-QC-passed full-cell object: ", source_file)
  message("\n===== corrected full-cell artificial profiles: ", study, " =====")
  object <- qs::qread(source_file); counts <- get_counts(object); meta <- object[[]]; meta$Cell <- rownames(meta)
  required <- c("Study", "DF_partition", "diagnostic_cluster", "historical_MG0", "MG00", "DoubletFinder_call", "scDblFinder_call", "Scrublet_call")
  if (length(setdiff(required, names(meta)))) stop(study, ": corrected source lacks ", paste(setdiff(required, names(meta)), collapse = ", "))
  if (!identical(colnames(counts), meta$Cell) || any(as.character(meta$Study) != study)) stop(study, ": input counts/metadata are not the complete study universe.")
  partitions <- split(meta$Cell, as.character(meta$DF_partition))
  if (any(lengths(partitions) < 300L)) stop(study, ": a technical calling partition has <300 all-QC-passed cells.")
  for (partition_i in seq_along(partitions)) {
    partition <- names(partitions)[[partition_i]]; cells <- partitions[[partition_i]]
    part_dir <- file.path(input_dir, paste0(safe, "__", safe_name(partition)))
    part_counts <- counts[, cells, drop = FALSE]
    detected <- Matrix::rowSums(part_counts) > 0
    part_counts_detected <- part_counts[detected, , drop = FALSE]
    # Same sampling scale as DoubletFinder pN=0.25. The artificial profiles
    # are generated from every all-QC-passed cell in the partition.
    n_df <- round(ncol(part_counts) / (1 - 0.25) - ncol(part_counts))
    set.seed(2026091300 + study_i * 100 + partition_i)
    parent_1 <- sample(cells, n_df, replace = TRUE); parent_2 <- sample(cells, n_df, replace = TRUE)
    df_counts <- (part_counts[, parent_1, drop = FALSE] + part_counts[, parent_2, drop = FALSE]) / 2
    df_ids <- paste0("DFart__", safe, "__", safe_name(partition), "__", seq_len(n_df)); colnames(df_counts) <- df_ids
    write_bundle(df_counts, part_dir, "DoubletFinder_artificial", df_ids)
    df_parents <- data.frame(artificial_cell_id = df_ids, Study = study, DF_partition = partition, parent_1 = parent_1, parent_2 = parent_2,
      parent_universe = "all QC-passed real cells in this study/technical partition", construction = "DoubletFinder-like pN=0.25: (raw parent-1 + raw parent-2)/2", stringsAsFactors = FALSE)
    write.csv(df_parents, file.path(part_dir, "DoubletFinder_artificial_parent_pairs.csv"), row.names = FALSE)
    all_df_parents[[length(all_df_parents) + 1L]] <- df_parents
    # Official scDblFinder helper on the same all-cell partition and original
    # global diagnostic clusters. Individual parent IDs are not exposed by it.
    set.seed(2026091800 + study_i * 100 + partition_i)
    clusters <- as.factor(as.character(meta[cells, "diagnostic_cluster"]))
    n_sc <- min(12000L, max(1500L, ceiling(ncol(part_counts_detected) * .50), 10L * length(unique(clusters))^2))
    sc_output <- scDblFinder::getArtificialDoublets(part_counts_detected, n = n_sc, clusters = clusters, propRandom = 0)
    sc_counts <- as(sc_output$counts, "dgCMatrix")
    sc_ids <- paste0("scDblart__", safe, "__", safe_name(partition), "__", seq_len(ncol(sc_counts))); colnames(sc_counts) <- sc_ids
    write_bundle(sc_counts, part_dir, "scDblFinder_artificial", sc_ids)
    sc_origins <- data.frame(artificial_cell_id = sc_ids, Study = study, DF_partition = partition,
      scDblFinder_origin_cluster = as.character(sc_output$origins), parent_universe = "all QC-passed real cells in this study/technical partition",
      construction = "scDblFinder::getArtificialDoublets (official helper; origin cluster retained)", stringsAsFactors = FALSE)
    write.csv(sc_origins, file.path(part_dir, "scDblFinder_artificial_origin_clusters.csv"), row.names = FALSE)
    all_sc_origins[[length(all_sc_origins) + 1L]] <- sc_origins
    # Scrublet receives the identical all-cell raw-count matrix. Python stage 2
    # calls its public simulate_doublets() method and records its parent index.
    scrublet_dir <- file.path(part_dir, "scrublet_input")
    write_bundle(part_counts_detected, scrublet_dir, "", cells)
    old_matrix <- file.path(scrublet_dir, "_counts_genes_by_cells.mtx")
    if (file.exists(old_matrix)) file.rename(old_matrix, file.path(scrublet_dir, "counts_genes_by_cells.mtx"))
    for (suffix in c("cells.tsv", "genes.tsv")) {
      old <- file.path(scrublet_dir, paste0("_", suffix)); if (file.exists(old)) file.rename(old, file.path(scrublet_dir, suffix))
    }
    all_manifest[[length(all_manifest) + 1L]] <- data.frame(Study = study, DF_partition = partition, safe_partition = basename(part_dir),
      n_all_QC_passed_real_cells = ncol(part_counts), n_detected_features = nrow(part_counts_detected),
      n_DoubletFinder_artificial = ncol(df_counts), n_scDblFinder_artificial = ncol(sc_counts), Scrublet_sim_doublet_ratio = 1/3,
      stringsAsFactors = FALSE)
    message(study, " / ", partition, ": ", ncol(part_counts), " complete QC-passed real parents")
    rm(part_counts, part_counts_detected, df_counts, sc_counts, sc_output); invisible(gc())
  }
  all_real_metadata[[study]] <- meta[, c("Cell", "Study", "DF_partition", "diagnostic_cluster", "historical_MG0", "MG00", "DoubletFinder_call", "scDblFinder_call", "Scrublet_call"), drop = FALSE]
  rm(object, counts, meta); invisible(gc())
}
write.csv(do.call(rbind, all_manifest), file.path(input_dir, "partition_manifest.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_df_parents), file.path(table_dir, "DoubletFinder_artificial_parent_pairs_all_QC_passed_cells.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_sc_origins), file.path(table_dir, "scDblFinder_artificial_origin_clusters_all_QC_passed_cells.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_real_metadata), file.path(table_dir, "all_QC_passed_real_cell_metadata_with_three_calls.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "stage1_R_sessionInfo.txt"))
writeLines(c(
  "Scope correction: artificial-profile parent pools are every QC-passed real cell in each study and technical calling partition.",
  "Historical MG0 and MG00 columns are retained only for later plot overlays. They do not restrict artificial parents or modify any call.",
  "DoubletFinder artificial: pN=0.25 mean profile. scDblFinder artificial: official getArtificialDoublets helper. Scrublet artificial: official simulate_doublets in stage 2.",
  "No cross-study or cross-partition artificial pairing is allowed."
), file.path(runtime_dir, "stage1_scope_manifest.txt"))
message("Corrected full-cell artificial stage 1 complete. Run 28.fullcell_artificial_doublet_validation_stage2_scrublet_24server.py, then stage 3.")
