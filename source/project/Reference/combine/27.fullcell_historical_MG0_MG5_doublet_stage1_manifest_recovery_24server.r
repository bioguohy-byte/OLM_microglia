#!/usr/bin/env Rscript

# Recover only the stage-1 manifests after the six complete all-QC-passed
# per-study objects have been written.  This deliberately loads one study at a
# time so the audit does not build a giant combined metadata data frame.

rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
invisible(gc())
suppressPackageStartupMessages({ library(qs); library(Seurat); library(SeuratObject); library(Matrix) })

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
result_root <- file.path(project_root, "combine/results/fullcell_historical_MG0_MG5_doublet_diagnosis_v2")
object_dir <- file.path(result_root, "object")
input_dir <- file.path(result_root, "scrublet_input")
table_dir <- file.path(result_root, "table")
runtime_dir <- file.path(result_root, "runtime")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(runtime_dir, recursive = TRUE, showWarnings = FALSE)
study_order <- c("Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
                 "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024")
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

manifest_rows <- list(); coverage_rows <- list(); i <- 0L
for (study in study_order) {
  safe <- safe_name(study)
  file <- file.path(object_dir, paste0(safe, "_all_QC_passed_fullcell_partial_calls.qs"))
  if (!file.exists(file)) stop("Missing completed full-cell partial object: ", file)
  object <- qs::qread(file)
  meta <- object[[]]
  if (!"DF_partition" %in% names(meta) || !"expected_doublet_rate" %in% names(meta)) stop(study, ": partition metadata missing")
  if (!"historical_MG_state" %in% names(meta)) stop(study, ": historical label mapping missing")
  partitions <- split(rownames(meta), as.character(meta$DF_partition))
  for (partition in names(partitions)) {
    cells <- partitions[[partition]]
    part_safe <- paste0(safe, "__", safe_name(partition))
    required_files <- file.path(input_dir, part_safe, c("counts_genes_by_cells.mtx", "cells.tsv", "genes.tsv"))
    if (any(!file.exists(required_files))) stop(study, " / ", partition, ": full-cell Scrublet bundle is incomplete")
    # The matrix file was created by the completed stage-1 call.  Counting its
    # accompanying feature file is sufficient here and avoids duplicating the
    # several-GB RNA matrix in memory during manifest recovery.
    n_genes <- length(readLines(file.path(input_dir, part_safe, "genes.tsv"), warn = FALSE))
    i <- i + 1L
    manifest_rows[[i]] <- data.frame(Study = study, safe_study = safe, DF_partition = partition,
      safe_partition = part_safe, n_cells = length(cells), n_genes = n_genes,
      expected_doublet_rate = unique(meta[cells, "expected_doublet_rate"]), stringsAsFactors = FALSE)
  }
  coverage_rows[[length(coverage_rows) + 1L]] <- data.frame(Study = study,
    n_all_QC_passed = nrow(meta), n_historical_MG0_MG5_mapped = sum(!is.na(meta$historical_MG_state)),
    n_historical_MG0_mapped = sum(as.logical(meta$historical_MG0)), n_MG00_mapped = sum(as.logical(meta$MG00)),
    stringsAsFactors = FALSE)
  rm(object, meta); invisible(gc())
}
manifest <- do.call(rbind, manifest_rows)
coverage <- do.call(rbind, coverage_rows)
utils::write.csv(manifest, file.path(input_dir, "scrublet_manifest.csv"), row.names = FALSE)
utils::write.csv(coverage, file.path(table_dir, "historical_MG0_MG5_and_MG00_mapping_coverage.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "stage1_manifest_recovery_R_sessionInfo.txt"))
writeLines(c(
  "Recovery writes manifests only after six all-QC-passed full-cell partial objects have completed.",
  "It does not recalculate, filter, remove or relabel any cell.",
  "Every Scrublet bundle is verified against the completed all-QC-passed study object."
), file.path(runtime_dir, "stage1_manifest_recovery_scope.txt"))
message("Recovered stage-1 full-cell manifest: ", nrow(manifest), " technical partitions across six studies.")
