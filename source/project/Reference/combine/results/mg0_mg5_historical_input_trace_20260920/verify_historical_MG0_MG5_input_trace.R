#!/usr/bin/env Rscript

# Read-only audit of the original six-study MG0-MG5 input chain.
# Checks exact cell-barcode equality at: source CellType selection -> *_mg.qs -> mg_cca_meta_v2.qs.

suppressPackageStartupMessages(library(qs))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: verify_historical_MG0_MG5_input_trace.R <report.txt>")
report_path <- args[[1]]
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
final_path <- file.path(root, "combine/data/processed_data/mg_cca_meta_v2.qs")

spec <- data.frame(
  study = c("Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021", "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"),
  source_path = c(
    file.path(root, "Su_et_al./data/processed_data/Su_CSC_alldata_v2.qs"),
    file.path(root, "Chia_et_al./data/processed_data/Chia_v2.qs"),
    file.path(root, "Fullard_et_al./data/processed_data/Fullard_v2.qs"),
    file.path(root, "Jakel_et_al./data/processed_data/Jakel_v2.qs"),
    file.path(root, "Grubman_et_al./data/processed_data/Grubman_v2.qs"),
    file.path(root, "Yao_et_al./data/processed_data/Yao_v2.qs")
  ),
  mg_path = c(
    file.path(root, "Su_et_al./data/processed_data/Su_mg.qs"),
    file.path(root, "Chia_et_al./data/processed_data/Chia_mg.qs"),
    file.path(root, "Fullard_et_al./data/processed_data/Fullard_mg.qs"),
    file.path(root, "Jakel_et_al./data/processed_data/Jakel_mg.qs"),
    file.path(root, "Grubman_et_al./data/processed_data/Grubman_mg.qs"),
    file.path(root, "Yao_et_al./data/processed_data/Yao_mg.qs")
  ),
  selected_labels = I(list(
    c("MG"), c("Microglia"), c("Mic"),
    c("Immune_cells", "ImOlGs", "Microglia_Macrophages"),
    c("Microglia"), c("Micro")
  )),
  stringsAsFactors = FALSE
)

final_obj <- qread(final_path)
if (!all(c("Study", "CellType", "CellType1") %in% colnames(final_obj@meta.data))) {
  stop("Historical MG0-MG5 object lacks required Study/CellType/CellType1 metadata.")
}
allowed_MG_states <- paste0("MG", 0:5)
rows <- list()

for (i in seq_len(nrow(spec))) {
  x <- spec[i, , drop = FALSE]
  source_obj <- qread(x$source_path[[1]])
  stored_mg <- qread(x$mg_path[[1]])
  if (!"CellType" %in% colnames(source_obj@meta.data)) stop(x$study[[1]], ": source CellType missing")

  selected_labels <- x$selected_labels[[1]]
  source_cells <- colnames(source_obj)[as.character(source_obj$CellType) %in% selected_labels]
  source_label_counts <- table(as.character(source_obj$CellType)[as.character(source_obj$CellType) %in% selected_labels])
  final_cells <- colnames(final_obj)[as.character(final_obj$Study) == x$study[[1]]]
  final_metadata <- final_obj@meta.data[final_cells, , drop = FALSE]
  mg0_mg5_counts <- table(as.character(final_metadata$CellType1))

  rows[[i]] <- data.frame(
    study = x$study[[1]],
    selection_metadata_column = "CellType",
    selected_source_labels = paste(selected_labels, collapse = "; "),
    source_cells_with_selected_labels = length(source_cells),
    stored_mg_object_cells = ncol(stored_mg),
    final_MG0_MG5_object_cells = length(final_cells),
    source_selection_equals_stored_mg = setequal(source_cells, colnames(stored_mg)),
    stored_mg_equals_final_study_cells = setequal(colnames(stored_mg), final_cells),
    final_study_CellType_within_selected_labels = all(as.character(final_metadata$CellType) %in% selected_labels),
    final_study_CellType1_all_MG0_to_MG5 = all(as.character(final_metadata$CellType1) %in% allowed_MG_states),
    source_label_count_detail = paste(paste0(names(source_label_counts), "=", as.integer(source_label_counts)), collapse = "; "),
    final_MG0_to_MG5_count_detail = paste(paste0(names(mg0_mg5_counts), "=", as.integer(mg0_mg5_counts)), collapse = "; "),
    stringsAsFactors = FALSE
  )
  rm(source_obj, stored_mg)
  gc()
}

audit <- do.call(rbind, rows)

sink(report_path)
cat("Historical six-study MG0-MG5 input trace (read-only)\n")
cat("Final object: ", final_path, "\n\n", sep = "")
print(audit, row.names = FALSE, right = FALSE)
cat("\nInterpretation rule:\n")
cat("A study has an exact historical chain only when both barcode-set equality columns are TRUE, final CellType remains within the selection labels, and final CellType1 is exclusively MG0-MG5.\n")
cat("This tests only the original project's CellType selection field. It does not claim that semantically related labels in other source metadata columns were selected.\n")
sink()
