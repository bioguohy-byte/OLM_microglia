#!/usr/bin/env Rscript

# One-time terminology migration for the six stage-2 OL–MG reference objects.
# This changes metadata labels only: canonical_MG_parent -> MG_parent,
# historical_MG0 -> MG0, and parent type text canonical_MG -> MG. Counts,
# barcodes, reductions, pairings, QC status and doublet calls remain identical.

local({
  suppressPackageStartupMessages({
    library(qs)
    library(SeuratObject)
    library(Matrix)
  })

  project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
  object_dir <- file.path(project_root, "combine/results/OL_MG_artificial_reference/02_artificial_profile_construction/object")
  runtime_dir <- file.path(project_root, "combine/results/OL_MG_artificial_reference/02_artificial_profile_construction/runtime")
  dir.create(runtime_dir, recursive = TRUE, showWarnings = FALSE)
  expected_studies <- c(
    "Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
    "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"
  )
  safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  object_files <- file.path(object_dir, paste0(vapply(expected_studies, safe_name, character(1)), "_real_cells_plus_artificial_profiles.qs"))
  if (any(!file.exists(object_files))) stop("Missing stage-2 objects: ", paste(object_files[!file.exists(object_files)], collapse = "; "))

  audit <- lapply(seq_along(object_files), function(i) {
    file <- object_files[[i]]
    object <- qs::qread(file)
    meta <- object@meta.data
    cells_before <- rownames(meta)
    count_dim_before <- dim(SeuratObject::GetAssayData(object, assay = "RNA", layer = "counts"))
    has_legacy_field <- "canonical_MG_parent" %in% names(meta)
    if (!"MG_parent" %in% names(meta)) {
      if (!has_legacy_field) stop(expected_studies[[i]], ": neither canonical_MG_parent nor MG_parent exists.")
      meta$MG_parent <- as.logical(meta$canonical_MG_parent)
    }
    meta$canonical_MG_parent <- NULL
    had_legacy_MG0_field <- "historical_MG0" %in% names(meta)
    if (!"MG0" %in% names(meta)) {
      if (!had_legacy_MG0_field) stop(expected_studies[[i]], ": neither historical_MG0 nor MG0 exists.")
      meta$MG0 <- as.logical(meta$historical_MG0)
    }
    meta$historical_MG0 <- NULL
    for (field in intersect(c("parent_1_type", "parent_2_type"), names(meta))) {
      meta[[field]] <- gsub("canonical_MG", "MG", as.character(meta[[field]]), fixed = TRUE)
    }
    object@meta.data <- meta[cells_before, , drop = FALSE]
    qs::qsave(object, file, preset = "high")

    check <- qs::qread(file)
    check_meta <- check@meta.data
    count_dim_after <- dim(SeuratObject::GetAssayData(check, assay = "RNA", layer = "counts"))
    if (!identical(cells_before, rownames(check_meta))) stop(expected_studies[[i]], ": barcode order changed during metadata migration.")
    if (!identical(count_dim_before, count_dim_after)) stop(expected_studies[[i]], ": count dimensions changed during metadata migration.")
    if ("canonical_MG_parent" %in% names(check_meta) || !"MG_parent" %in% names(check_meta) || "historical_MG0" %in% names(check_meta) || !"MG0" %in% names(check_meta)) stop(expected_studies[[i]], ": terminology migration did not validate.")
    type_values <- unlist(check_meta[, intersect(c("parent_1_type", "parent_2_type"), names(check_meta)), drop = FALSE], use.names = FALSE)
    if (any(grepl("canonical_MG", type_values, fixed = TRUE))) stop(expected_studies[[i]], ": legacy parent type text remains.")
    data.frame(
      Study = expected_studies[[i]], object_file = file, had_legacy_MG_parent_field = has_legacy_field,
      had_legacy_MG0_field = had_legacy_MG0_field, MG_parent_present = TRUE, MG0_present = TRUE, barcode_order_preserved = TRUE,
      count_dimensions_preserved = TRUE, stringsAsFactors = FALSE
    )
  })
  audit <- do.call(rbind, audit)
  write.csv(audit, file.path(runtime_dir, "MG_terminology_in_place_migration_audit.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "MG_terminology_in_place_migration_sessionInfo.txt"))
  message("Terminology migration complete: six stage-2 OL–MG .qs objects now use MG_parent and MG0 only.")
})
