#!/usr/bin/env Rscript

# Read-only audit of feature axes / identifiers for the original MG integration
# and the current count-based reconstruction.  This script does not overwrite
# any Seurat object or analysis result; it writes only audit CSV files.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(qs)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
audit_root <- file.path(project_root, "combine/results/full_count_reanalysis/00_preflight/feature_identifier_audit")
dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)

# Fetch every counts layer without assuming the Seurat v4/v5 assay class.
# A v5 object can retain one counts.* layer per source after a merge; the
# reported feature set is their union, whereas the matrix for a single layer
# is reported separately below.
extract_counts_layers <- function(object, assay = "RNA") {
  if (!assay %in% Seurat::Assays(object)) stop("RNA assay is absent.")
  assay_object <- object[[assay]]
  if (inherits(assay_object, "Assay5")) {
    available_layers <- SeuratObject::Layers(assay_object)
    count_layers <- grep("^counts($|\\.)", available_layers, value = TRUE)
    if (!length(count_layers)) stop("Assay5 has no counts layer.")
    matrices <- lapply(count_layers, function(layer_name) {
      SeuratObject::LayerData(assay_object, layer = layer_name)
    })
    names(matrices) <- count_layers
  } else {
    count_matrix <- tryCatch(
      SeuratObject::LayerData(assay_object, layer = "counts"),
      error = function(e) NULL
    )
    if (is.null(count_matrix) || !length(count_matrix)) {
      count_matrix <- Seurat::GetAssayData(object, assay = assay, layer = "counts")
    }
    matrices <- list(counts = count_matrix)
  }
  matrices
}

feature_style <- function(features) {
  is_ensembl <- grepl("^ENSG[0-9]+(?:\\.[0-9]+)?$", features)
  is_ensembl_versioned <- grepl("^ENSG[0-9]+\\.[0-9]+$", features)
  is_symbol_like <- grepl("^[A-Za-z][A-Za-z0-9-]*$", features)
  is_upper_symbol <- grepl("^[A-Z][A-Z0-9-]*$", features)
  is_numeric <- grepl("^[0-9]+$", features)
  is_makeunique_suffix <- grepl("\\.[0-9]+$", features) & !is_ensembl_versioned
  data.frame(
    n_features = length(features),
    n_unique_features = length(unique(features)),
    n_duplicate_feature_names = sum(duplicated(features)),
    n_ensembl = sum(is_ensembl),
    n_ensembl_versioned = sum(is_ensembl_versioned),
    n_symbol_like = sum(is_symbol_like),
    n_uppercase_symbol_like = sum(is_upper_symbol),
    n_numeric_only = sum(is_numeric),
    n_terminal_dot_number_non_ensembl = sum(is_makeunique_suffix),
    n_mitochondrial_symbol = sum(grepl("^MT-", features)),
    n_ribosomal_symbol = sum(grepl("^RP[SL]", features)),
    stringsAsFactors = FALSE
  )
}

# A final ".number" is not sufficient evidence of a generated suffix: many
# legitimate lncRNA symbols use that syntax (for example AC005550.3).  Mark a
# feature only when its immediate prefix is also present in *the same* feature
# set.  This is a conservative, diagnostic-only flag; no name is changed and
# no rows are summed in this audit.
probable_makeunique_suffix <- function(features) {
  is_ensembl_versioned <- grepl("^ENSG[0-9]+\\.[0-9]+$", features)
  has_terminal_number <- grepl("\\.[0-9]+$", features) & !is_ensembl_versioned
  prefix <- sub("\\.[0-9]+$", "", features)
  has_terminal_number & prefix %in% features
}

object_manifest <- data.frame(
  Object_group = c(
    rep("Current_reanalysis_stage1", 6),
    rep("Original_microglia_input", 6),
    rep("Original_CCA_integration", 2)
  ),
  Label = c(
    "Su", "Chia", "Fullard", "Jakel", "Grubman", "Yao",
    "Su", "Chia", "Fullard", "Jakel", "Grubman", "Yao",
    "mg_before_inte", "mg_cca_meta_v2"
  ),
  File = c(
    file.path(project_root, "combine/results/full_count_reanalysis/01_per_study/object", c(
      "Su_et_al_2022_full_count_processed.qs",
      "Chia_et_al_2025_full_count_processed.qs",
      "Fullard_et_al_2021_full_count_processed.qs",
      "Jakel_et_al_2019_full_count_processed.qs",
      "Grubman_et_al_2019_full_count_processed.qs",
      "Yao_et_al_2024_full_count_processed.qs"
    )),
    file.path(project_root, c(
      "Su_et_al./data/processed_data/Su_mg.qs",
      "Chia_et_al./data/processed_data/Chia_mg.qs",
      "Fullard_et_al./data/processed_data/Fullard_mg.qs",
      "Jakel_et_al./data/processed_data/Jakel_mg.qs",
      "Grubman_et_al./data/processed_data/Grubman_mg.qs",
      "Yao_et_al./data/processed_data/Yao_mg.qs",
      "combine/data/processed_data/mg_before_inte.qs",
      "combine/data/processed_data/mg_cca_meta_v2.qs"
    ))
  ),
  stringsAsFactors = FALSE
)

if (any(!file.exists(object_manifest$File))) {
  missing_files <- object_manifest$File[!file.exists(object_manifest$File)]
  stop("Missing audit inputs:\n", paste(missing_files, collapse = "\n"))
}

object_rows <- list()
layer_rows <- list()
feature_sets <- list()
for (i in seq_len(nrow(object_manifest))) {
  entry <- object_manifest[i, ]
  message(sprintf("[%d/%d] Reading %s: %s", i, nrow(object_manifest), entry$Object_group, entry$Label))
  object <- qs::qread(entry$File)
  count_layers <- extract_counts_layers(object)

  layer_features <- lapply(count_layers, rownames)
  layer_cells <- lapply(count_layers, colnames)
  union_features <- sort(unique(unlist(layer_features, use.names = FALSE)))
  union_cells <- unique(unlist(layer_cells, use.names = FALSE))
  object_key <- paste(entry$Object_group, entry$Label, sep = "::")
  feature_sets[[object_key]] <- union_features

  layer_rows[[length(layer_rows) + 1L]] <- do.call(rbind, lapply(seq_along(count_layers), function(j) {
    count_matrix <- count_layers[[j]]
    data.frame(
      Object_group = entry$Object_group,
      Label = entry$Label,
      Count_layer = names(count_layers)[j],
      n_features = nrow(count_matrix),
      n_cells = ncol(count_matrix),
      n_nonzero = length(count_matrix@x),
      stringsAsFactors = FALSE
    )
  }))

  styles <- feature_style(union_features)
  probable_suffix_count <- sum(probable_makeunique_suffix(union_features))
  object_rows[[length(object_rows) + 1L]] <- cbind(
    data.frame(
      Object_group = entry$Object_group,
      Label = entry$Label,
      File = entry$File,
      Seurat_class = paste(class(object), collapse = ";"),
      Assays = paste(Seurat::Assays(object), collapse = ";"),
      Default_assay = Seurat::DefaultAssay(object),
      n_object_cells = ncol(object),
      n_metadata_rows = nrow(object@meta.data),
      n_counts_layers = length(count_layers),
      n_count_features_union = length(union_features),
      n_count_cells_union = length(union_cells),
      n_probable_makeunique_suffix = probable_suffix_count,
      stringsAsFactors = FALSE
    ),
    styles
  )

  rm(count_layers, layer_features, layer_cells, union_features, union_cells, object)
  invisible(gc())
}

object_summary <- do.call(rbind, object_rows)
layer_summary <- do.call(rbind, layer_rows)

# Compare six original MG input objects with six corresponding objects produced
# by the current stage-1 procedure.  A selected-cell subset does not alter its
# RNA feature axis, so these feature-set unions exactly describe the matrix that
# would be formed after the user-approved MG selection.
study_order <- c("Su", "Chia", "Fullard", "Jakel", "Grubman", "Yao")
group_overlap_rows <- list()
for (group_name in c("Current_reanalysis_stage1", "Original_microglia_input")) {
  group_keys <- paste(group_name, study_order, sep = "::")
  sets <- feature_sets[group_keys]
  if (any(vapply(sets, is.null, logical(1)))) stop("Feature set missing for ", group_name)
  union_set <- sort(unique(unlist(sets, use.names = FALSE)))
  intersection_set <- Reduce(intersect, sets)
  n_studies_present <- vapply(union_set, function(gene) sum(vapply(sets, function(s) gene %in% s, logical(1))), integer(1))
  group_overlap_rows[[length(group_overlap_rows) + 1L]] <- data.frame(
    Object_group = group_name,
    n_studies = length(sets),
    feature_union = length(union_set),
    feature_intersection_all_six = length(intersection_set),
    stringsAsFactors = FALSE
  )
  study_presence <- as.data.frame(table(n_studies_present), stringsAsFactors = FALSE)
  names(study_presence) <- c("n_studies_containing_feature", "n_features")
  study_presence$Object_group <- group_name
  utils::write.csv(study_presence,
    file.path(audit_root, paste0(group_name, "_feature_presence_distribution.csv")),
    row.names = FALSE
  )

  intersection_matrix <- outer(seq_along(sets), seq_along(sets), Vectorize(function(a, b) length(intersect(sets[[a]], sets[[b]]))))
  union_matrix <- outer(seq_along(sets), seq_along(sets), Vectorize(function(a, b) length(union(sets[[a]], sets[[b]]))))
  dimnames(intersection_matrix) <- list(study_order, study_order)
  dimnames(union_matrix) <- list(study_order, study_order)
  jaccard_matrix <- intersection_matrix / union_matrix
  utils::write.csv(intersection_matrix,
    file.path(audit_root, paste0(group_name, "_pairwise_feature_intersection.csv"))
  )
  utils::write.csv(jaccard_matrix,
    file.path(audit_root, paste0(group_name, "_pairwise_feature_jaccard.csv"))
  )
}
group_overlap_summary <- do.call(rbind, group_overlap_rows)

# Directly contrast the current and original six-study unions and their latest
# original CCA object.  This is the requested check against the 50,670 row
# current merged count matrix.
current_union <- sort(unique(unlist(feature_sets[paste("Current_reanalysis_stage1", study_order, sep = "::")], use.names = FALSE)))
original_input_union <- sort(unique(unlist(feature_sets[paste("Original_microglia_input", study_order, sep = "::")], use.names = FALSE)))
original_final_features <- feature_sets[["Original_CCA_integration::mg_cca_meta_v2"]]
original_pre_integration_features <- feature_sets[["Original_CCA_integration::mg_before_inte"]]
comparison_summary <- data.frame(
  Item = c(
    "Current_reanalysis_six_study_feature_union",
    "Original_six_MG_input_feature_union",
    "Original_mg_before_inte_RNA_counts",
    "Original_mg_cca_meta_v2_RNA_counts",
    "Current_vs_original_input_union_intersection",
    "Current_vs_original_final_RNA_intersection",
    "Features_present_in_original_input_but_absent_from_current"
  ),
  n_features = c(
    length(current_union),
    length(original_input_union),
    length(original_pre_integration_features),
    length(original_final_features),
    length(intersect(current_union, original_input_union)),
    length(intersect(current_union, original_final_features)),
    length(setdiff(original_input_union, current_union))
  ),
  stringsAsFactors = FALSE
)

per_study_current_vs_original <- do.call(rbind, lapply(study_order, function(study) {
  current_features <- feature_sets[[paste("Current_reanalysis_stage1", study, sep = "::")]]
  original_features <- feature_sets[[paste("Original_microglia_input", study, sep = "::")]]
  data.frame(
    Study = study,
    current_n_features = length(current_features),
    original_n_features = length(original_features),
    exact_intersection = length(intersect(current_features, original_features)),
    current_only = length(setdiff(current_features, original_features)),
    original_only = length(setdiff(original_features, current_features)),
    current_probable_makeunique_suffix = sum(probable_makeunique_suffix(current_features)),
    original_probable_makeunique_suffix = sum(probable_makeunique_suffix(original_features)),
    stringsAsFactors = FALSE
  )
}))
missing_from_current <- data.frame(
  Feature = sort(setdiff(original_input_union, current_union)),
  stringsAsFactors = FALSE
)
probable_makeunique_candidates <- do.call(rbind, lapply(names(feature_sets), function(object_key) {
  features <- feature_sets[[object_key]]
  candidate <- probable_makeunique_suffix(features)
  if (!any(candidate)) {
    return(data.frame(
      Object = character(0),
      Feature = character(0),
      Immediate_prefix_present_in_same_object = character(0),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    Object = object_key,
    Feature = features[candidate],
    Immediate_prefix_present_in_same_object = sub("\\.[0-9]+$", "", features[candidate]),
    stringsAsFactors = FALSE
  )
}))

utils::write.csv(object_summary, file.path(audit_root, "object_RNA_counts_summary.csv"), row.names = FALSE)
utils::write.csv(layer_summary, file.path(audit_root, "RNA_counts_layer_summary.csv"), row.names = FALSE)
utils::write.csv(group_overlap_summary, file.path(audit_root, "six_study_feature_overlap_summary.csv"), row.names = FALSE)
utils::write.csv(comparison_summary, file.path(audit_root, "current_vs_original_feature_comparison.csv"), row.names = FALSE)
utils::write.csv(per_study_current_vs_original, file.path(audit_root, "per_study_current_vs_original_feature_comparison.csv"), row.names = FALSE)
utils::write.csv(missing_from_current, file.path(audit_root, "features_present_in_original_MG_input_but_absent_from_current.csv"), row.names = FALSE)
utils::write.csv(probable_makeunique_candidates, file.path(audit_root, "probable_makeunique_suffix_candidates.csv"), row.names = FALSE)

message("\n=== Feature-axis audit complete ===")
print(comparison_summary, row.names = FALSE)
message("Results written to: ", audit_root)
