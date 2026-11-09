#!/usr/bin/env Rscript

# Export per-cell QC fields from the exact six full-cell objects used for the
# source-specific CellType QC violin plots.  The source objects are read only.
# For concordance with those plots, historical MG cells are relabelled only in
# the exported QC_CellType column: MG0 -> OLM; MG1-MG5 are retained.

suppressPackageStartupMessages(library(qs))

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
output_dir <- file.path(
  project_root,
  "combine/results/qc_per_cell_nCount_nFeature_by_source_celltype_20260920"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
per_study_dir <- file.path(output_dir, "per_study_csv")
dir.create(per_study_dir, recursive = TRUE, showWarnings = FALSE)

study_manifest <- data.frame(
  Study = c("Su", "Chia", "Fullard", "Jakel", "Grubman", "Yao"),
  historical_Study = c(
    "Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
    "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"
  ),
  source_path = c(
    file.path(project_root, "Su_et_al./data/processed_data/Su_CSC_alldata_v2.qs"),
    file.path(project_root, "Chia_et_al./data/processed_data/Chia_v2.qs"),
    file.path(project_root, "Fullard_et_al./data/processed_data/Fullard_v2.qs"),
    file.path(project_root, "Jakel_et_al./data/processed_data/Jakel_v2.qs"),
    file.path(project_root, "Grubman_et_al./data/processed_data/Grubman_v2.qs"),
    file.path(project_root, "Yao_et_al./data/processed_data/Yao_v2.qs")
  ),
  stringsAsFactors = FALSE
)

final_mg_path <- file.path(project_root, "combine/data/processed_data/mg_cca_meta_v2.qs")
final_mg_metadata <- qread(final_mg_path)[[]]
required_reference_columns <- c("Study", "CellType1")
if (!all(required_reference_columns %in% colnames(final_mg_metadata))) {
  stop("The final microglia reference lacks Study and/or CellType1 metadata.")
}

required_source_columns <- c("CellType", "nCount_RNA", "nFeature_RNA")
expected_header <- c(
  "Study", "historical_Study", "cell_id", "CellType_original", "QC_CellType",
  "nCount_RNA", "nFeature_RNA"
)

total_cells <- 0L
for (i in seq_len(nrow(study_manifest))) {
  study <- study_manifest$Study[[i]]
  historical_study <- study_manifest$historical_Study[[i]]
  source_object <- qread(study_manifest$source_path[[i]])
  source_metadata <- source_object[[]]

  if (!all(required_source_columns %in% colnames(source_metadata))) {
    stop(study, " is missing required QC metadata: ",
         paste(setdiff(required_source_columns, colnames(source_metadata)), collapse = ", "))
  }

  source_cells <- rownames(source_metadata)
  historical_cells <- rownames(final_mg_metadata)[final_mg_metadata$Study == historical_study]
  missing_historical_cells <- setdiff(historical_cells, source_cells)
  if (length(missing_historical_cells) > 0L) {
    stop(study, " lacks historical MG barcodes (first missing barcode: ",
         missing_historical_cells[[1]], ").")
  }

  qc_celltype <- as.character(source_metadata$CellType)
  qc_celltype[is.na(qc_celltype) | qc_celltype == ""] <- "Unannotated"
  historical_subtype <- as.character(
    final_mg_metadata[historical_cells, "CellType1", drop = TRUE]
  )
  qc_celltype[match(historical_cells, source_cells)] <- ifelse(
    historical_subtype == "MG0", "OLM", historical_subtype
  )

  export_table <- data.frame(
    Study = study,
    historical_Study = historical_study,
    cell_id = source_cells,
    CellType_original = as.character(source_metadata$CellType),
    QC_CellType = qc_celltype,
    nCount_RNA = as.numeric(source_metadata$nCount_RNA),
    nFeature_RNA = as.numeric(source_metadata$nFeature_RNA),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  stopifnot(identical(colnames(export_table), expected_header))

  output_path <- file.path(
    per_study_dir, paste0(study, "_per_cell_nCount_nFeature.csv")
  )
  temporary_path <- paste0(output_path, ".tmp")
  if (file.exists(temporary_path)) unlink(temporary_path)
  write.table(
    export_table,
    file = temporary_path,
    sep = ",",
    quote = TRUE,
    row.names = FALSE,
    col.names = TRUE,
    append = FALSE,
    na = ""
  )
  if (!file.rename(temporary_path, output_path)) {
    stop("Could not finalize the QC CSV at ", output_path)
  }
  total_cells <- total_cells + nrow(export_table)
  message(study, ": wrote ", nrow(export_table), " cells to ", output_path)
  rm(source_object, source_metadata, export_table)
  gc(verbose = FALSE)
}

message("Finished ", total_cells, " cells. Output directory: ", per_study_dir)
