#!/usr/bin/env Rscript

# Inspect the locally installed, versioned babelgene Human-to-Mouse mapping
# table format before applying it to the Human MG0 top-50 signature.

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/three_slice_ST_MG0_MG5_dual_module_score_20260916")
dir.create(file.path(project_dir, "results", "03_Human_MG0_top50_strict_mouse_ortholog_module_score", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library"), file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_libs"), .libPaths()))

suppressPackageStartupMessages(library(babelgene))

test_human_symbols <- c("RNF220", "ERBB3", "KIAA1755")
test_mapping <- babelgene::orthologs(genes = test_human_symbols, species = "mouse", human = TRUE)
write.table(
  as.data.frame(test_mapping),
  file.path(project_dir, "results", "03_Human_MG0_top50_strict_mouse_ortholog_module_score", "tables", "babelgene_test_human_to_mouse_mapping.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
capture.output(str(test_mapping), file = file.path(project_dir, "logs", "03a_check_babelgene_mapping_format_structure.txt"))
capture.output(sessionInfo(), file = file.path(project_dir, "logs", "03a_check_babelgene_mapping_format_sessionInfo.txt"))
