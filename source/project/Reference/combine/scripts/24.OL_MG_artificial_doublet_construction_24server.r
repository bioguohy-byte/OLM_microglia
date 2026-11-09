#!/usr/bin/env Rscript

# OL–MG artificial-doublet reference, stage 2
#
# 1) Reads reviewer-confirmed mature OL clusters and the previously
#    reviewer-confirmed MG clusters.
# 2) Draws only within-study, within-technical-batch pairs from current
#    DoubletFinder singlets.  It never changes an existing real-cell object.
# 3) Creates raw-count-sum and parent-depth-matched versions of OL–MG,
#    OL–OL and MG–MG artificial profiles.
# 4) Reuses the PCA that was fitted earlier on QC-passed REAL cells, fits a
#    fresh projection-capable UMAP only on those real-cell PCA coordinates,
#    then projects artificial profiles into that fixed reference.  Artificial
#    profiles therefore cannot alter the real-cell manifold, QC, labels, or
#    doublet calls.
#
# The saved per-study .qs object contains QC-passed real cells and the new
# artificial profiles with raw integer counts, audit metadata, fixed-reference
# PCA/UMAP coordinates, and no modified real-cell annotations.

local({
  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
    library(qs)
    library(ggplot2)
    library(uwot)
  })

  set.seed(20260909)
  project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
  combine_root <- file.path(project_root, "combine")
  per_study_dir <- file.path(combine_root, "results/full_count_reanalysis/01_per_study/object")
  manual_root <- file.path(combine_root, "results/OL_MG_artificial_reference/01_OL_parent_review/manual_selections")
  ol_selection_file <- file.path(manual_root, "OL_parent_cluster_selection_APPROVED.csv")
  mg_selection_file <- file.path(combine_root, "results/full_count_reanalysis/manual_selections/per_study_Microglia_cluster_selection_TEMPLATE.csv")
  historical_file <- file.path(combine_root, "results/MG0_MG5_doublet_diagnosis/table/six_study_cell_metadata_with_all_doublet_calls.csv")
  mg00_file <- file.path(combine_root, "results/full_count_reanalysis/03_current_MG00_doublet_diagnosis/table/six_study_cell_metadata_with_all_three_doublet_calls_MG00.csv")
  result_root <- file.path(combine_root, "results/OL_MG_artificial_reference/02_artificial_profile_construction")
  figure_dir <- file.path(result_root, "figure")
  table_dir <- file.path(result_root, "table")
  object_dir <- file.path(result_root, "object")
  checkpoint_dir <- file.path(result_root, "checkpoint")
  runtime_dir <- file.path(result_root, "runtime")
  for (d in c(figure_dir, table_dir, object_dir, checkpoint_dir, runtime_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  study_order <- c(
    "Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
    "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"
  )
  # Every eligible technical batch contributes at least one pair.  Capping at
  # 3 pair-sets/batch yields at most 18 profiles/batch (3 lineage pair types ×
  # 2 count constructions), retaining batch representation while preventing a
  # few large batches from dominating the artificial reference.
  max_pairs_per_batch <- 3L
  seed_base <- 20260909L
  safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  is_true <- function(x) !is.na(x) & trimws(tolower(as.character(x))) %in% c("true", "t", "1", "yes", "y")

  save_plot_both <- function(plot, stem, width, height) {
    ggsave(paste0(stem, ".png"), plot = plot, width = width, height = height, dpi = 320, bg = "white", limitsize = FALSE)
    ggsave(paste0(stem, ".pdf"), plot = plot, width = width, height = height, device = grDevices::pdf, bg = "white", limitsize = FALSE)
    if (!file.exists(paste0(stem, ".png")) || !file.exists(paste0(stem, ".pdf"))) stop("Incomplete PNG/PDF output: ", stem)
  }

  get_layer <- function(object, layer) {
    out <- tryCatch(SeuratObject::LayerData(object, assay = "RNA", layer = layer), error = function(e) NULL)
    if (is.null(out)) out <- tryCatch(SeuratObject::GetAssayData(object, assay = "RNA", layer = layer), error = function(e) NULL)
    if (is.null(out)) out <- tryCatch(SeuratObject::GetAssayData(object, assay = "RNA", slot = layer), error = function(e) NULL)
    if (is.null(out)) stop("Cannot read RNA ", layer, " layer.")
    if (!inherits(out, "sparseMatrix")) out <- Matrix::Matrix(out, sparse = TRUE)
    out
  }

  assert_integer_counts <- function(x, label) {
    if (length(x@x) && any(abs(x@x - round(x@x)) > 1e-8)) stop(label, " contains non-integer count values.")
  }

  make_binary_plot <- function(object, status_col, positive_label, title, colour) {
    object[[status_col]] <- factor(
      ifelse(as.logical(object[[status_col]][, 1]), positive_label, "Other current singlet"),
      levels = c("Other current singlet", positive_label)
    )
    # Direct ggplot2 rendering preserves the explicit positive colour on this
    # server (the Seurat::DimPlot rendering path previously flattened it).
    coordinates <- as.data.frame(Seurat::Embeddings(object, reduction = "umap"))
    if (ncol(coordinates) < 2L) stop(title, ": UMAP has fewer than two dimensions.")
    coordinates$UMAP_1_display <- coordinates[[1L]]
    coordinates$UMAP_2_display <- coordinates[[2L]]
    coordinates$plot_label <- object[[status_col]][rownames(coordinates), 1]
    if (anyNA(coordinates$plot_label)) stop(title, ": missing plotting label.")
    ggplot2::ggplot(coordinates, ggplot2::aes(x = UMAP_1_display, y = UMAP_2_display)) +
      ggplot2::geom_point(
        data = coordinates[coordinates$plot_label == "Other current singlet", , drop = FALSE],
        colour = "#D9D9D9", size = 0.15, alpha = 0.80
      ) +
      ggplot2::geom_point(
        data = coordinates[coordinates$plot_label == positive_label, , drop = FALSE],
        colour = colour, size = 0.17, alpha = 1
      ) +
      ggtitle(title) +
      labs(x = "UMAP 1", y = "UMAP 2") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  }

  downsample_sparse_profile <- function(profile, target_umi, seed) {
    profile <- as(profile, "dgCMatrix")
    total <- sum(profile@x)
    target_umi <- as.integer(max(1, min(round(target_umi), total)))
    if (target_umi >= total) return(profile)
    set.seed(seed)
    profile@x <- as.numeric(stats::rmultinom(1L, size = target_umi, prob = profile@x))
    Matrix::drop0(profile)
  }

  pair_same_type <- function(cells, n_pairs) {
    if (length(cells) < 2L || n_pairs < 1L) return(NULL)
    parent_1 <- sample(cells, n_pairs, replace = FALSE)
    parent_2 <- sample(cells, n_pairs, replace = FALSE)
    identical_parent <- which(parent_1 == parent_2)
    if (length(identical_parent)) {
      for (i in identical_parent) parent_2[i] <- sample(setdiff(cells, parent_1[i]), 1L)
    }
    data.frame(parent_1 = parent_1, parent_2 = parent_2, stringsAsFactors = FALSE)
  }

  append_profile_pair <- function(profile_list, metadata_list, audit_list, counts, study, batch, pair_set, parent_1, parent_2, parent_1_type, parent_2_type, pair_index, seed_offset) {
    raw <- as(counts[, parent_1, drop = FALSE] + counts[, parent_2, drop = FALSE], "dgCMatrix")
    raw_total <- as.numeric(sum(raw@x))
    parent_1_umi <- as.numeric(sum(counts[, parent_1, drop = FALSE]@x))
    parent_2_umi <- as.numeric(sum(counts[, parent_2, drop = FALSE]@x))
    depth_target <- as.integer(round(mean(c(parent_1_umi, parent_2_umi))))
    for (construction in c("raw_count_sum", "parent_depth_matched")) {
      profile <- if (construction == "raw_count_sum") raw else downsample_sparse_profile(raw, depth_target, seed = seed_offset + pair_index)
      profile_id <- paste0("ART__", safe_name(study), "__", safe_name(batch), "__", pair_set, "__", construction, "__", sprintf("%05d", pair_index))
      colnames(profile) <- profile_id
      profile_list[[length(profile_list) + 1L]] <- profile
      metadata_list[[length(metadata_list) + 1L]] <- data.frame(
        Cell = profile_id, Study = study, batch_id = batch, cell_origin = "artificial_profile",
        profile_class = pair_set, construction = construction, parent_1 = parent_1, parent_2 = parent_2,
        parent_1_type = parent_1_type, parent_2_type = parent_2_type,
        parent_1_umi = parent_1_umi, parent_2_umi = parent_2_umi,
        raw_profile_umi = raw_total, depth_matched_target_umi = depth_target,
        stringsAsFactors = FALSE
      )
      audit_list[[length(audit_list) + 1L]] <- data.frame(
        Study = study, batch_id = batch, profile_id = profile_id, profile_class = pair_set,
        construction = construction, parent_1 = parent_1, parent_2 = parent_2,
        parent_1_type = parent_1_type, parent_2_type = parent_2_type,
        parent_1_umi = parent_1_umi, parent_2_umi = parent_2_umi,
        raw_profile_umi = raw_total, depth_matched_target_umi = depth_target,
        random_seed = if (construction == "parent_depth_matched") seed_offset + pair_index else NA_integer_,
        same_study = TRUE, same_batch = TRUE, parent_1_current_DoubletFinder_singlet = TRUE,
        parent_2_current_DoubletFinder_singlet = TRUE, stringsAsFactors = FALSE
      )
    }
    list(profiles = profile_list, metadata = metadata_list, audit = audit_list)
  }

  project_artificial_pca <- function(art_counts, reference, features, dims, block_size = 400L) {
    # Align all three matrices by explicit feature names.  This matters because
    # a Seurat reduction can retain a feature-loading order that differs from
    # the assay's row order even when the selected HVG set is the same.
    all_loadings <- as.matrix(Seurat::Loadings(reference, reduction = "pca"))
    features <- intersect(features, rownames(all_loadings))
    features <- intersect(features, rownames(art_counts))
    if (length(features) < 2L) stop("Too few shared variable features for artificial-profile PCA projection.")
    reference_data <- get_layer(reference, "data")[features, , drop = FALSE]
    center <- Matrix::rowMeans(reference_data)
    variance <- Matrix::rowMeans(reference_data ^ 2) - center ^ 2
    scale_value <- sqrt(pmax(variance, 1e-8))
    loadings <- all_loadings[features, dims, drop = FALSE]
    pcs <- matrix(NA_real_, nrow = ncol(art_counts), ncol = length(dims),
                  dimnames = list(colnames(art_counts), colnames(loadings)))
    starts <- seq.int(1L, ncol(art_counts), by = block_size)
    for (start in starts) {
      ids <- start:min(ncol(art_counts), start + block_size - 1L)
      block <- as.matrix(art_counts[rownames(loadings), ids, drop = FALSE])
      libsize <- colSums(block)
      block <- log1p(t(t(block) / pmax(libsize, 1)) * 10000)
      block <- sweep(block, 1L, center, "-")
      block <- sweep(block, 1L, scale_value, "/")
      block[!is.finite(block)] <- 0
      # pmin/pmax would drop the matrix dimension on this R build; clip
      # in-place so the features × profiles geometry is preserved for %*%.
      block[block < -10] <- -10
      block[block > 10] <- 10
      if (ncol(t(block)) != nrow(loadings)) stop("Artificial PCA projection feature mismatch after explicit alignment: ", ncol(t(block)), " versus ", nrow(loadings))
      pcs[ids, ] <- t(block) %*% loadings
    }
    pcs
  }

  construct_profile_bundle <- function(final_counts, final_meta, ol_cells, mg_cells, study, study_i, max_pairs_per_batch, seed_base, pairing_plan_file = NULL) {
    parent_cells <- unique(c(ol_cells, mg_cells))
    parent_batch <- as.character(final_meta[parent_cells, "batch_id"])
    if (any(is.na(parent_batch) | parent_batch == "")) stop(study, ": parent cells with missing batch_id cannot be paired.")
    batches <- sort(unique(parent_batch))
    batch_audit <- do.call(rbind, lapply(seq_along(batches), function(batch_i) {
      batch <- batches[[batch_i]]
      ol_batch <- intersect(ol_cells, final_meta$Cell[as.character(final_meta$batch_id) == batch])
      mg_batch <- intersect(mg_cells, final_meta$Cell[as.character(final_meta$batch_id) == batch])
      n_pairs <- min(length(ol_batch), length(mg_batch), max_pairs_per_batch)
      data.frame(
        Study = study, batch_index = batch_i, batch_id = batch, mature_OL_singlets_available = length(ol_batch),
        MG_singlets_available = length(mg_batch), n_pair_sets_constructed = n_pairs,
        expected_artificial_profiles = 6L * n_pairs, max_pairs_per_batch = max_pairs_per_batch,
        stringsAsFactors = FALSE
      )
    }))
    if (!is.null(pairing_plan_file)) write.csv(batch_audit, pairing_plan_file, row.names = FALSE)
    message(study, ": pairing plan has ", sum(batch_audit$n_pair_sets_constructed > 0L),
            " eligible batches, ", sum(batch_audit$n_pair_sets_constructed), " pair-sets, and ",
            sum(batch_audit$expected_artificial_profiles), " artificial profiles.")
    profiles <- list(); artificial_meta <- list(); profile_audit <- list(); pair_counter <- 0L
    for (plan_i in which(batch_audit$n_pair_sets_constructed > 0L)) {
      batch <- batch_audit$batch_id[[plan_i]]
      batch_i <- batch_audit$batch_index[[plan_i]]
      n_pairs <- batch_audit$n_pair_sets_constructed[[plan_i]]
      ol_batch <- intersect(ol_cells, final_meta$Cell[as.character(final_meta$batch_id) == batch])
      mg_batch <- intersect(mg_cells, final_meta$Cell[as.character(final_meta$batch_id) == batch])
      set.seed(seed_base + study_i * 1000L + batch_i)
      ol_mg <- data.frame(parent_1 = sample(ol_batch, n_pairs, replace = FALSE), parent_2 = sample(mg_batch, n_pairs, replace = FALSE), stringsAsFactors = FALSE)
      ol_ol <- pair_same_type(ol_batch, n_pairs)
      mg_mg <- pair_same_type(mg_batch, n_pairs)
      pair_sets <- list(OL_MG = list(pairs = ol_mg, type_1 = "mature_OL", type_2 = "MG"),
                        OL_OL = list(pairs = ol_ol, type_1 = "mature_OL", type_2 = "mature_OL"),
                        MG_MG = list(pairs = mg_mg, type_1 = "MG", type_2 = "MG"))
      for (pair_set in names(pair_sets)) {
        spec <- pair_sets[[pair_set]]
        if (is.null(spec$pairs) || !nrow(spec$pairs)) next
        for (i in seq_len(nrow(spec$pairs))) {
          pair_counter <- pair_counter + 1L
          appended <- append_profile_pair(
            profiles, artificial_meta, profile_audit, final_counts, study, batch, pair_set,
            spec$pairs$parent_1[[i]], spec$pairs$parent_2[[i]], spec$type_1, spec$type_2,
            pair_counter, seed_offset = seed_base + study_i * 100000L + batch_i * 1000L
          )
          profiles <- appended$profiles; artificial_meta <- appended$metadata; profile_audit <- appended$audit
        }
      }
    }
    if (!length(profiles)) stop(study, ": no same-batch OL–MG parent pairs were available.")
    artificial_counts <- as(do.call(cbind, profiles), "dgCMatrix")
    assert_integer_counts(artificial_counts, paste0(study, " artificial profiles"))
    artificial_meta <- do.call(rbind, artificial_meta)
    rownames(artificial_meta) <- artificial_meta$Cell
    profile_audit <- do.call(rbind, profile_audit)
    batch_audit <- do.call(rbind, batch_audit)
    if (!all(profile_audit$same_study & profile_audit$same_batch & profile_audit$parent_1_current_DoubletFinder_singlet & profile_audit$parent_2_current_DoubletFinder_singlet)) stop(study, ": parent audit failed.")
    list(artificial_counts = artificial_counts, artificial_meta = artificial_meta, profile_audit = profile_audit, batch_audit = batch_audit)
  }

  normalize_batch_audit <- function(x) {
    # qs 0.27.3 can deserialize a data.frame nested in a list as a transposed
    # character matrix.  Recover the named columns before resuming a checkpoint.
    if (is.matrix(x)) {
      field_names <- rownames(x)
      if (is.null(field_names)) stop("Checkpoint batch audit matrix is missing field names.")
      x <- as.data.frame(t(x), stringsAsFactors = FALSE)
      names(x) <- field_names
    }
    if (!is.data.frame(x)) stop("Checkpoint batch audit has an unsupported type: ", paste(class(x), collapse = "/"))
    # Checkpoints generated before the terminology revision used the same
    # values under canonical_MG_singlets_available.  Migrate that one
    # metadata heading on read; pair identities and counts are untouched.
    if ("canonical_MG_singlets_available" %in% names(x) && !"MG_singlets_available" %in% names(x)) {
      names(x)[names(x) == "canonical_MG_singlets_available"] <- "MG_singlets_available"
    }
    numeric_fields <- intersect(c("batch_index", "mature_OL_singlets_available", "MG_singlets_available", "n_pair_sets_constructed", "expected_artificial_profiles", "max_pairs_per_batch"), names(x))
    for (field in numeric_fields) x[[field]] <- as.integer(x[[field]])
    x
  }

  required_files <- c(ol_selection_file, mg_selection_file, historical_file, mg00_file)
  if (any(!file.exists(required_files))) stop("Missing required input(s): ", paste(required_files[!file.exists(required_files)], collapse = "; "))
  if (!requireNamespace("uwot", quietly = TRUE)) stop("Package uwot is required for fixed-reference UMAP projection.")

  ol_selection <- read.csv(ol_selection_file, check.names = FALSE, stringsAsFactors = FALSE)
  mg_selection <- read.csv(mg_selection_file, check.names = FALSE, stringsAsFactors = FALSE)
  historical <- read.csv(historical_file, check.names = FALSE, stringsAsFactors = FALSE)
  mg00 <- read.csv(mg00_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("Study", "analysis_cluster", "select_as_mature_OL") %in% names(ol_selection))) stop("Approved OL selection file is malformed.")
  if (!all(c("Study", "analysis_cluster", "select_as_Microglia") %in% names(mg_selection))) stop("Existing MG selection file is malformed.")
  if (!all(c("Study", "Cell", "CellType1") %in% names(historical))) stop("Original MG0 mapping table is malformed.")
  if (!all(c("Study", "original_barcode", "MG00_status", "DoubletFinder_call", "scDblFinder_call", "Scrublet_call") %in% names(mg00))) stop("MG00 mapping table is malformed.")
  ol_selection <- ol_selection[is_true(ol_selection$select_as_mature_OL), , drop = FALSE]
  mg_selection <- mg_selection[is_true(mg_selection$select_as_Microglia), , drop = FALSE]
  if (!setequal(study_order, unique(ol_selection$Study))) stop("Approved OL selection must include every study.")
  if (!setequal(study_order, unique(mg_selection$Study))) stop("Existing Microglia selection must include every study.")
  write.csv(ol_selection, file.path(table_dir, "reviewer_approved_mature_OL_clusters.csv"), row.names = FALSE)
  write.csv(mg_selection, file.path(table_dir, "read_only_MG_clusters.csv"), row.names = FALSE)

  historical$key <- paste(historical$Study, historical$Cell, sep = "||")
  historical_mg0_keys <- historical$key[historical$CellType1 == "MG0"]
  mg00$key <- paste(mg00$Study, mg00$original_barcode, sep = "||")
  if (anyDuplicated(mg00$key)) stop("MG00 table has duplicated Study/original barcode mappings.")

  all_audit <- list()
  all_mapping <- list()
  all_study_summary <- list()
  all_figure_manifest <- list()

  for (study_i in seq_along(study_order)) {
    study <- study_order[[study_i]]
    study_safe <- safe_name(study)
    message("\n===== Artificial OL–MG construction: ", study, " =====")
    completed_object_file <- file.path(object_dir, paste0(study_safe, "_real_cells_plus_artificial_profiles.qs"))
    completed_audit_file <- file.path(table_dir, paste0(study_safe, "_artificial_profile_parent_pair_audit.csv"))
    completed_batch_file <- file.path(table_dir, paste0(study_safe, "_same_batch_pairing_availability.csv"))
    completed_mapping_file <- file.path(table_dir, paste0(study_safe, "_barcode_mapping_audit.csv"))
    completed_ol_stem <- file.path(figure_dir, paste0(study_safe, "_01_confirmed_mature_OL_UMAP"))
    completed_mg_stem <- file.path(figure_dir, paste0(study_safe, "_02_confirmed_MG_UMAP"))
    if (all(file.exists(c(completed_object_file, completed_audit_file, completed_batch_file, completed_mapping_file,
                          paste0(completed_ol_stem, ".png"), paste0(completed_ol_stem, ".pdf"),
                          paste0(completed_mg_stem, ".png"), paste0(completed_mg_stem, ".pdf"))))) {
      # A completed per-study object plus its three audits is an atomic resume
      # boundary.  Reconstruct the aggregate summaries from those audits rather
      # than reloading and re-embedding the already completed study.
      profile_audit_done <- read.csv(completed_audit_file, check.names = FALSE, stringsAsFactors = FALSE)
      batch_audit_done <- read.csv(completed_batch_file, check.names = FALSE, stringsAsFactors = FALSE)
      mapping_audit_done <- read.csv(completed_mapping_file, check.names = FALSE, stringsAsFactors = FALSE)
      all_audit[[study]] <- profile_audit_done
      all_mapping[[study]] <- mapping_audit_done
      all_study_summary[[study]] <- data.frame(
        Study = study, QC_passed_real_cells = nrow(mapping_audit_done),
        current_DoubletFinder_singlet_parent_cells = sum(mapping_audit_done$source_final_singlet %in% TRUE),
        mature_OL_singlet_parents = sum(batch_audit_done$mature_OL_singlets_available),
        MG_singlet_parents = sum(batch_audit_done$MG_singlets_available),
        artificial_profiles = nrow(profile_audit_done),
        same_batch_pair_sets = sum(batch_audit_done$n_pair_sets_constructed),
        MG0_mapped = sum(mapping_audit_done$MG0 %in% TRUE),
        MG00_mapped = sum(mapping_audit_done$MG00 %in% TRUE), stringsAsFactors = FALSE
      )
      all_figure_manifest[[length(all_figure_manifest) + 1L]] <- data.frame(Study = study, figure = "confirmed_mature_OL_UMAP", stem = completed_ol_stem, stringsAsFactors = FALSE)
      all_figure_manifest[[length(all_figure_manifest) + 1L]] <- data.frame(Study = study, figure = "confirmed_MG_UMAP", stem = completed_mg_stem, stringsAsFactors = FALSE)
      message(study, ": completed object and audits detected; skipping recomputation.")
      next
    }
    final_file <- file.path(per_study_dir, paste0(study_safe, "_full_count_processed.qs"))
    qc_file <- file.path(per_study_dir, paste0(study_safe, "_QC_before_doublet_removal.qs"))
    if (!file.exists(final_file) || !file.exists(qc_file)) stop(study, ": missing final-singlet or QC-passed input object.")

    final_object <- qs::qread(final_file)
    qc_object <- qs::qread(qc_file)
    if (!all(c("RNA") %in% SeuratObject::Assays(final_object)) || !all(c("RNA") %in% SeuratObject::Assays(qc_object))) stop(study, ": RNA assay missing.")
    final_meta <- final_object[[]]
    qc_meta <- qc_object[[]]
    final_meta$Cell <- rownames(final_meta)
    qc_meta$Cell <- rownames(qc_meta)
    required_final <- c("analysis_cluster", "batch_id", "DoubletFinder_call", "original_barcode")
    required_qc <- c("batch_id", "DoubletFinder_call", "original_barcode")
    if (length(setdiff(required_final, names(final_meta)))) stop(study, ": final object missing required metadata.")
    if (length(setdiff(required_qc, names(qc_meta)))) stop(study, ": QC object missing required metadata.")
    if (!all(as.character(final_meta$DoubletFinder_call) == "Singlet")) stop(study, ": final parent object contains non-singlets.")
    final_counts <- get_layer(final_object, "counts")
    qc_counts <- get_layer(qc_object, "counts")
    assert_integer_counts(final_counts, paste0(study, " final parents"))
    assert_integer_counts(qc_counts, paste0(study, " QC-passed real cells"))
    if (!setequal(rownames(final_counts), rownames(qc_counts))) stop(study, ": parent and QC expression feature sets differ.")
    # The two saved objects normally inherit the same feature order.  Avoid a
    # full sparse-matrix duplicate unless an explicit realignment is required.
    if (!identical(rownames(final_counts), rownames(qc_counts))) {
      final_counts <- final_counts[rownames(qc_counts), , drop = FALSE]
    }

    ol_clusters <- as.character(ol_selection$analysis_cluster[
      ol_selection$Study == study & ol_selection$select_as_mature_OL %in% TRUE
    ])
    mg_clusters <- as.character(mg_selection$analysis_cluster[
      mg_selection$Study == study & mg_selection$select_as_Microglia %in% TRUE
    ])
    if (!length(ol_clusters) || !length(mg_clusters)) stop(study, ": approved OL or MG cluster selection is empty.")
    final_meta$mature_OL_parent <- as.character(final_meta$analysis_cluster) %in% ol_clusters
    final_meta$MG_parent <- as.character(final_meta$analysis_cluster) %in% mg_clusters
    if (any(final_meta$mature_OL_parent & final_meta$MG_parent)) stop(study, ": reviewer-defined OL and MG parent sets overlap.")
    if (!"umap" %in% names(final_object@reductions)) stop(study, ": final object has no UMAP for definition display.")
    final_object$mature_OL_parent <- final_meta$mature_OL_parent
    final_object$MG_parent <- final_meta$MG_parent

    # These two figures document the human OL and MG definitions on
    # the exact per-study final-singlet UMAP used for prior cluster review.
    p_ol <- make_binary_plot(final_object, "mature_OL_parent", "Confirmed mature OL", paste0(study, ": reviewer-confirmed mature OL"), "#B2182B")
    p_mg <- make_binary_plot(final_object, "MG_parent", "Confirmed MG", paste0(study, ": prior MG definition"), "#2166AC")
    ol_stem <- file.path(figure_dir, paste0(study_safe, "_01_confirmed_mature_OL_UMAP"))
    mg_stem <- file.path(figure_dir, paste0(study_safe, "_02_confirmed_MG_UMAP"))
    save_plot_both(p_ol, ol_stem, 8, 6.5)
    save_plot_both(p_mg, mg_stem, 8, 6.5)
    print(p_ol)
    print(p_mg)
    all_figure_manifest[[length(all_figure_manifest) + 1L]] <- data.frame(Study = study, figure = "confirmed_mature_OL_UMAP", stem = ol_stem, stringsAsFactors = FALSE)
    all_figure_manifest[[length(all_figure_manifest) + 1L]] <- data.frame(Study = study, figure = "confirmed_MG_UMAP", stem = mg_stem, stringsAsFactors = FALSE)

    ol_cells <- final_meta$Cell[final_meta$mature_OL_parent]
    mg_cells <- final_meta$Cell[final_meta$MG_parent]
    if (!length(ol_cells) || !length(mg_cells)) stop(study, ": zero confirmed OL or MG singlet parent cells.")
    checkpoint_file <- file.path(checkpoint_dir, paste0(study_safe, "_artificial_profile_bundle.qs"))
    if (file.exists(checkpoint_file)) {
      profile_bundle <- qs::qread(checkpoint_file)
      profile_bundle$batch_audit <- normalize_batch_audit(profile_bundle$batch_audit)
      for (bundle_table in c("artificial_meta", "profile_audit")) {
        if (!bundle_table %in% names(profile_bundle)) next
        for (type_field in intersect(c("parent_1_type", "parent_2_type"), names(profile_bundle[[bundle_table]]))) {
          profile_bundle[[bundle_table]][[type_field]] <- gsub("canonical_MG", "MG", as.character(profile_bundle[[bundle_table]][[type_field]]), fixed = TRUE)
        }
      }
      message(study, ": resumed audited artificial profiles from checkpoint.")
    } else {
      profile_bundle <- construct_profile_bundle(
        final_counts, final_meta, ol_cells, mg_cells, study, study_i, max_pairs_per_batch, seed_base,
        pairing_plan_file = file.path(table_dir, paste0(study_safe, "_same_batch_pairing_plan_preconstruction.csv"))
      )
      # The profile bundle is a durable checkpoint: a later interruption never
      # requires re-sampling parent pairs, and all random choices remain auditable.
      qs::qsave(profile_bundle, checkpoint_file, preset = "balanced")
      # qs can mutate a nested data.frame during serialization on the installed
      # version, so normalize it here as well as on later checkpoint reads.
      profile_bundle$batch_audit <- normalize_batch_audit(profile_bundle$batch_audit)
    }
    artificial_counts <- profile_bundle$artificial_counts
    artificial_meta <- profile_bundle$artificial_meta
    profile_audit <- profile_bundle$profile_audit
    batch_audit <- profile_bundle$batch_audit

    # Map original-pipeline MG0 and current MG00 strictly by Study + preserved source barcode.
    qc_key <- paste(study, as.character(qc_meta$original_barcode), sep = "||")
    qc_meta$MG0 <- qc_key %in% historical_mg0_keys
    mg00_idx <- match(qc_key, mg00$key)
    qc_meta$MG00 <- !is.na(mg00_idx) & mg00$MG00_status[mg00_idx] == "MG00 (manual cluster 2 transfer)"
    qc_meta$scDblFinder_call <- ifelse(!is.na(mg00_idx), as.character(mg00$scDblFinder_call[mg00_idx]), NA_character_)
    qc_meta$Scrublet_call <- ifelse(!is.na(mg00_idx), as.character(mg00$Scrublet_call[mg00_idx]), NA_character_)
    qc_meta$external_reanalysis_DoubletFinder_call <- ifelse(!is.na(mg00_idx), as.character(mg00$DoubletFinder_call[mg00_idx]), NA_character_)
    final_key <- paste(study, as.character(final_meta$original_barcode), sep = "||")
    qc_meta$source_final_singlet <- qc_key %in% final_key
    qc_meta$mature_OL_parent <- qc_meta$source_final_singlet & qc_key %in% final_key[final_meta$mature_OL_parent]
    qc_meta$MG_parent <- qc_meta$source_final_singlet & qc_key %in% final_key[final_meta$MG_parent]

    # The QC-before-removal objects already retain their real-cell PCA.  Reuse
    # it rather than rebuilding a dense scaled matrix.  This is still a strict
    # real-cell-only reference: a fresh UMAP transform model is fitted solely
    # to those pre-existing QC-passed real-cell PCA coordinates.
    if (!"pca" %in% names(qc_object@reductions)) stop(study, ": QC-passed object has no PCA reduction to use as a real-cell reference.")
    pca_loadings <- Seurat::Loadings(qc_object, reduction = "pca")
    features <- rownames(pca_loadings)
    if (length(features) < 2L) stop(study, ": QC-passed PCA has too few loading features.")
    reference_pca_all <- Seurat::Embeddings(qc_object, reduction = "pca")
    dims <- seq_len(min(30L, ncol(reference_pca_all)))
    reference_pca <- reference_pca_all[, dims, drop = FALSE]
    set.seed(seed_base + study_i)
    umap_fit <- uwot::umap(
      as.matrix(reference_pca), n_neighbors = 30L, n_components = 2L,
      metric = "euclidean", n_threads = 4L, ret_model = TRUE, verbose = FALSE
    )
    # On the installed uwot version, ret_model=TRUE returns a model-list whose
    # embedding is one element of that list (rather than a nested $model).
    # Pass that verified model-list directly to umap_transform().
    if (is.null(umap_fit$embedding)) stop(study, ": real-cell UMAP embedding was not retained; artificial profiles cannot be projected safely.")
    umap_model <- umap_fit
    reference_umap <- umap_fit$embedding
    rownames(reference_umap) <- rownames(reference_pca)
    colnames(reference_umap) <- c("UMAP_1", "UMAP_2")
    artificial_pca <- project_artificial_pca(artificial_counts, qc_object, features, dims)
    artificial_umap <- uwot::umap_transform(artificial_pca, umap_model)
    rownames(artificial_umap) <- rownames(artificial_pca)
    colnames(artificial_umap) <- colnames(reference_umap)

    real_meta_out <- data.frame(
      Cell = colnames(qc_counts), Study = study, batch_id = as.character(qc_meta$batch_id),
      original_barcode = as.character(qc_meta$original_barcode), cell_origin = "real_QC_passed",
      profile_class = "real", construction = "real", parent_1 = NA_character_, parent_2 = NA_character_,
      parent_1_type = NA_character_, parent_2_type = NA_character_, parent_1_umi = NA_real_, parent_2_umi = NA_real_,
      raw_profile_umi = as.numeric(Matrix::colSums(qc_counts)), depth_matched_target_umi = NA_real_,
      DoubletFinder_call = as.character(qc_meta$DoubletFinder_call), MG0 = qc_meta$MG0,
      MG00 = qc_meta$MG00, mature_OL_parent = qc_meta$mature_OL_parent, MG_parent = qc_meta$MG_parent,
      source_final_singlet = qc_meta$source_final_singlet,
      scDblFinder_call = qc_meta$scDblFinder_call, Scrublet_call = qc_meta$Scrublet_call,
      external_reanalysis_DoubletFinder_call = qc_meta$external_reanalysis_DoubletFinder_call,
      stringsAsFactors = FALSE, row.names = colnames(qc_counts)
    )
    artificial_meta$original_barcode <- NA_character_
    artificial_meta$DoubletFinder_call <- NA_character_
    artificial_meta$MG0 <- NA
    artificial_meta$MG00 <- NA
    artificial_meta$mature_OL_parent <- NA
    artificial_meta$MG_parent <- NA
    artificial_meta$source_final_singlet <- NA
    artificial_meta$scDblFinder_call <- NA_character_
    artificial_meta$Scrublet_call <- NA_character_
    artificial_meta$external_reanalysis_DoubletFinder_call <- NA_character_
    artificial_meta <- artificial_meta[, colnames(real_meta_out), drop = FALSE]
    combined_counts <- as(cbind(qc_counts, artificial_counts), "dgCMatrix")
    combined_meta <- rbind(real_meta_out, artificial_meta)
    if (!identical(colnames(combined_counts), rownames(combined_meta))) stop(study, ": combined count columns and metadata rows do not match.")
    combined_object <- Seurat::CreateSeuratObject(counts = combined_counts, meta.data = combined_meta, project = study)
    combined_object[["pca_real_reference_projection"]] <- SeuratObject::CreateDimReducObject(
      embeddings = rbind(reference_pca, artificial_pca),
      key = "realrefPC_", assay = "RNA"
    )
    combined_object[["umap_real_reference_projection"]] <- SeuratObject::CreateDimReducObject(
      embeddings = rbind(reference_umap, artificial_umap),
      key = "realrefUMAP_", assay = "RNA"
    )
    combined_object$artificial_profile <- combined_object$cell_origin == "artificial_profile"
    combined_object$nCount_RNA <- as.numeric(Matrix::colSums(combined_counts))
    combined_object$nFeature_RNA <- as.numeric(Matrix::colSums(combined_counts > 0))
    qs::qsave(combined_object, file.path(object_dir, paste0(study_safe, "_real_cells_plus_artificial_profiles.qs")), preset = "high")

    mapping_audit <- data.frame(
      Study = study, Cell = rownames(real_meta_out), original_barcode = real_meta_out$original_barcode,
      source_final_singlet = real_meta_out$source_final_singlet, MG0 = real_meta_out$MG0,
      MG00 = real_meta_out$MG00, current_DoubletFinder_call = real_meta_out$DoubletFinder_call,
      stringsAsFactors = FALSE
    )
    write.csv(profile_audit, file.path(table_dir, paste0(study_safe, "_artificial_profile_parent_pair_audit.csv")), row.names = FALSE)
    write.csv(batch_audit, file.path(table_dir, paste0(study_safe, "_same_batch_pairing_availability.csv")), row.names = FALSE)
    write.csv(mapping_audit, file.path(table_dir, paste0(study_safe, "_barcode_mapping_audit.csv")), row.names = FALSE)
    all_audit[[study]] <- profile_audit
    all_mapping[[study]] <- mapping_audit
    all_study_summary[[study]] <- data.frame(
      Study = study, QC_passed_real_cells = ncol(qc_counts), current_DoubletFinder_singlet_parent_cells = ncol(final_counts),
      mature_OL_singlet_parents = length(ol_cells), MG_singlet_parents = length(mg_cells),
      artificial_profiles = ncol(artificial_counts), same_batch_pair_sets = sum(batch_audit$n_pair_sets_constructed),
      MG0_mapped = sum(real_meta_out$MG0), MG00_mapped = sum(real_meta_out$MG00),
      stringsAsFactors = FALSE
    )
    rm(final_object, qc_object, final_counts, qc_counts, final_meta, qc_meta, combined_object,
       combined_counts, artificial_counts, artificial_meta, profile_audit, batch_audit, profile_bundle,
       pca_loadings, reference_pca_all, reference_pca, reference_umap, artificial_pca, artificial_umap,
       umap_model, umap_fit)
    gc(verbose = FALSE)
  }

  figure_manifest <- do.call(rbind, all_figure_manifest)
  figure_manifest$png_exists <- file.exists(paste0(figure_manifest$stem, ".png"))
  figure_manifest$pdf_exists <- file.exists(paste0(figure_manifest$stem, ".pdf"))
  if (!all(figure_manifest$png_exists & figure_manifest$pdf_exists)) stop("At least one OL/MG definition UMAP lacks PNG or PDF.")
  write.csv(figure_manifest, file.path(table_dir, "OL_MG_definition_UMAP_PNG_PDF_manifest.csv"), row.names = FALSE)
  write.csv(do.call(rbind, all_audit), file.path(table_dir, "all_studies_artificial_profile_parent_pair_audit.csv"), row.names = FALSE)
  write.csv(do.call(rbind, all_mapping), file.path(table_dir, "all_studies_barcode_mapping_audit.csv"), row.names = FALSE)
  write.csv(do.call(rbind, all_study_summary), file.path(table_dir, "artificial_profile_construction_summary.csv"), row.names = FALSE)
  manifest <- data.frame(
    key = c("script", "run_timestamp", "seed", "max_pairs_per_batch", "OL_parent_selection", "MG_selection", "artificial_profile_classes", "construction_versions", "real_reference_fit", "artificial_effect_on_real_cells", "output_objects"),
    value = c("24.OL_MG_artificial_doublet_construction_24server.r", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), seed_base,
              max_pairs_per_batch, ol_selection_file, mg_selection_file, "OL_MG;OL_OL;MG_MG", "raw_count_sum;parent_depth_matched",
              "Reused existing QC-passed real-cell PCA; UMAP fitted only on its coordinates before artificial projection", "none", "one .qs per study plus one resumable artificial-profile bundle per study"), stringsAsFactors = FALSE
  )
  write.csv(manifest, file.path(runtime_dir, "stage2_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "sessionInfo.txt"))
  message("\nStage 2 complete. Artificial references are in: ", object_dir)
  message("Every parent pair is audited and restricted to the same study, technical batch, and current DoubletFinder singlet parent pool.")
})
