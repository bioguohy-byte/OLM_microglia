#!/usr/bin/env Rscript

# Read-only preflight for the supplied integrated Human/Mouse/Monkey object.
# This script identifies the existing species and annotation columns before any
# marker calculation, module scoring, or re-labelling is attempted.

options(stringsAsFactors = FALSE)

input_qs <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/02.cross_species_integrated_noref_Human_Mouse_Monkey.qs")
project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/OLM_MG0_top50_cross_species_20260915")

for (one_dir in c("scripts", "results", "figures", "logs", "tmp", ".cache", ".config")) {
  dir.create(file.path(project_dir, one_dir), recursive = TRUE, showWarnings = FALSE)
}
Sys.setenv(
  R_LIBS = file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_libs"),
  R_HISTFILE = file.path(project_dir, "logs", ".Rhistory"),
  XDG_CACHE_HOME = file.path(project_dir, ".cache"),
  XDG_CONFIG_HOME = file.path(project_dir, ".config"),
  R_USER_CACHE_DIR = file.path(project_dir, ".cache", "R"),
  TMPDIR = file.path(project_dir, "tmp")
)
.libPaths(c(file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(SeuratObject)
})

if (!file.exists(input_qs)) stop("Input QS object is absent: ", input_qs)

cross_species <- qread(input_qs)
if (!inherits(cross_species, "Seurat")) stop("The supplied QS object is not a Seurat object. Class: ", paste(class(cross_species), collapse = ", "))

metadata <- cross_species[[]]
object_summary <- data.frame(
  input_qs = input_qs,
  object_class = paste(class(cross_species), collapse = "; "),
  cells = ncol(cross_species),
  features = nrow(cross_species),
  assays = paste(Assays(cross_species), collapse = "; "),
  reductions = paste(Reductions(cross_species), collapse = "; "),
  active_assay = DefaultAssay(cross_species),
  active_identity_levels = paste(utils::head(levels(Idents(cross_species)), 50), collapse = "; "),
  stringsAsFactors = FALSE
)
write.table(object_summary, file.path(project_dir, "results", "00_object_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

metadata_summary <- data.frame(
  metadata_column = colnames(metadata),
  class = vapply(metadata, function(x) paste(class(x), collapse = "; "), character(1)),
  non_missing_cells = vapply(metadata, function(x) sum(!is.na(x)), numeric(1)),
  unique_values = vapply(metadata, function(x) length(unique(x[!is.na(x)])), numeric(1)),
  example_values = vapply(metadata, function(x) paste(utils::head(sort(unique(as.character(x[!is.na(x)]))), 10), collapse = " | "), character(1)),
  stringsAsFactors = FALSE
)
write.table(metadata_summary, file.path(project_dir, "results", "00_metadata_column_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

candidate_pattern <- "species|organism|human|mouse|monkey|cell.?type|annotation|celltype|cluster|ident|mg0|mg_0|ol[mi]"
candidate_columns <- metadata_summary$metadata_column[grepl(candidate_pattern, metadata_summary$metadata_column, ignore.case = TRUE, perl = TRUE)]
candidate_value_table <- data.frame()
for (one_column in candidate_columns) {
  one_values <- as.character(metadata[[one_column]])
  one_table <- sort(table(one_values, useNA = "ifany"), decreasing = TRUE)
  candidate_value_table <- rbind(
    candidate_value_table,
    data.frame(
      metadata_column = one_column,
      value = names(one_table),
      cells = as.integer(one_table),
      stringsAsFactors = FALSE
    )
  )
}
write.table(candidate_value_table, file.path(project_dir, "results", "00_candidate_species_and_annotation_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

identity_table <- sort(table(as.character(Idents(cross_species)), useNA = "ifany"), decreasing = TRUE)
write.table(
  data.frame(active_identity = names(identity_table), cells = as.integer(identity_table), stringsAsFactors = FALSE),
  file.path(project_dir, "results", "00_active_identity_values.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

feature_preview <- data.frame(gene = utils::head(rownames(cross_species), 100), stringsAsFactors = FALSE)
write.table(feature_preview, file.path(project_dir, "results", "00_feature_name_preview.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(project_dir, "logs", "00_preflight_sessionInfo.txt"))
writeLines(format(Sys.time(), tz = "Asia/Shanghai"), file.path(project_dir, "logs", "00_preflight_completed_at.txt"))
message("Preflight complete: ", ncol(cross_species), " cells; ", nrow(cross_species), " features; candidate annotation columns: ", paste(candidate_columns, collapse = ", "))
