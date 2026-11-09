#!/usr/bin/env Rscript

# Spatial microglia: Human MG0 top-50 converted strictly to mouse orthologs.
# The three supplied MERFISH slice objects are mouse.  Human MG0 top-50 genes
# are therefore mapped to mouse orthologs before scoring.
# Only unambiguous Human-to-Mouse orthologs supported by at least two HCOP
# source databases are retained; their mouse MGI symbols are then intersected
# with each native MERFISH RNA panel.
# No input object is overwritten.

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/three_slice_ST_MG0_MG5_dual_module_score_20260916")
result_dir <- file.path(project_dir, "results", "03_Human_MG0_top50_strict_mouse_ortholog_module_score")
figure_dir <- file.path(result_dir, "figures")
table_dir <- file.path(result_dir, "tables")
object_dir <- file.path(result_dir, "objects")
log_dir <- file.path(project_dir, "logs")
tmp_dir <- file.path(project_dir, "tmp")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir, XDG_CACHE_HOME = file.path(project_dir, "cache"))
dir.create(Sys.getenv("XDG_CACHE_HOME"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library"), file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(babelgene)
  library(SCP)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260916)

slice_paths <- c(
  s2r1 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r1_microglia_ST_MG0_MG5_from_v2.rds"),
  s2r2 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r2_microglia_ST_MG0_MG5_from_v2.rds"),
  s2r3 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r3_microglia_ST_MG0_MG5_from_v2.rds")
)
human_top50_tsv <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/OLM_MG0_top50_cross_species_20260915/results/01_human_MG0_top50_score_cross_species_OLM/tables/Human_MG0_FindMarkers_top50_reference_signature.tsv")
target_order <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")
st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35", "ST_MG1" = "#4DBBD5", "ST_MG2" = "#00A087",
  "ST_MG3" = "#3C5488", "ST_MG4" = "#F39B7F", "ST_MG5" = "#7E6148"
)

human_top50 <- read.delim(human_top50_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (!"gene" %in% colnames(human_top50) || nrow(human_top50) != 50L) {
  stop("The Human MG0 reference file must have exactly 50 rows and a gene column.")
}
human_top50_genes <- unique(as.character(human_top50$gene))
if (length(human_top50_genes) != 50L) stop("Human MG0 top-50 genes are not unique.")

# babelgene 22.9 bundles the HGNC Comparison of Orthology Predictions (HCOP)
# database locally. Therefore mapping does not depend on a live service. To
# form the strict signature, require: a non-empty mouse MGI symbol; support
# from at least two source databases; and exactly one unique candidate mouse
# symbol per human reference gene.
ortholog_raw <- as.data.frame(babelgene::orthologs(genes = human_top50_genes, species = "mouse", human = TRUE), stringsAsFactors = FALSE)
ortholog_raw <- ortholog_raw[ortholog_raw$human_symbol %in% human_top50_genes, , drop = FALSE]
ortholog_raw$symbol <- as.character(ortholog_raw$symbol)
ortholog_raw$support_n <- as.numeric(ortholog_raw$support_n)
strict_candidate_rows <- ortholog_raw[
  !is.na(ortholog_raw$symbol) & ortholog_raw$symbol != "" & ortholog_raw$support_n >= 2,
  , drop = FALSE
]
strict_candidate_rows <- strict_candidate_rows[!duplicated(strict_candidate_rows[, c("human_symbol", "symbol")]), , drop = FALSE]

# A human gene is included only if it has one unique strict mouse MGI symbol.
# Genes with zero or multiple candidates remain in the audit but are excluded.
strict_mouse_symbol_count <- table(strict_candidate_rows$human_symbol)
strict_mouse_symbol_count <- strict_mouse_symbol_count[human_top50_genes]
strict_mouse_symbol_count[is.na(strict_mouse_symbol_count)] <- 0L
strict_mouse_symbol_count <- as.integer(strict_mouse_symbol_count)
names(strict_mouse_symbol_count) <- human_top50_genes

ortholog_audit <- human_top50
ortholog_audit$human_gene <- as.character(ortholog_audit$gene)
ortholog_audit$strict_mouse_symbol_count <- strict_mouse_symbol_count[ortholog_audit$human_gene]
ortholog_audit$mouse_mgi_symbol <- vapply(ortholog_audit$human_gene, function(gene_name) {
  matching_symbols <- unique(strict_candidate_rows$symbol[strict_candidate_rows$human_symbol == gene_name])
  if (length(matching_symbols) == 1L) matching_symbols else NA_character_
}, character(1))
ortholog_audit$strict_mapping_status <- ifelse(
  ortholog_audit$strict_mouse_symbol_count == 1L,
  "retained_unique_mouse_ortholog_support_n_ge_2",
  ifelse(ortholog_audit$strict_mouse_symbol_count == 0L, "excluded_no_mouse_ortholog_with_support_n_ge_2", "excluded_multiple_mouse_ortholog_candidates")
)
ortholog_audit$raw_babelgene_rows <- vapply(ortholog_audit$human_gene, function(gene_name) sum(ortholog_raw$human_symbol == gene_name), integer(1))
ortholog_audit$babelgene_candidate_mouse_symbols <- vapply(ortholog_audit$human_gene, function(gene_name) {
  paste(sort(unique(ortholog_raw$symbol[ortholog_raw$human_symbol == gene_name & !is.na(ortholog_raw$symbol) & ortholog_raw$symbol != ""])), collapse = ";")
}, character(1))
ortholog_audit$babelgene_candidate_support_n <- vapply(ortholog_audit$human_gene, function(gene_name) {
  paste(sort(unique(ortholog_raw$support_n[ortholog_raw$human_symbol == gene_name & !is.na(ortholog_raw$symbol) & ortholog_raw$symbol != ""])), collapse = ";")
}, character(1))
write.table(ortholog_raw, file.path(table_dir, "Human_MG0_top50_babelgene_all_human_to_mouse_homolog_rows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ortholog_audit, file.path(table_dir, "Human_MG0_top50_strict_human_to_mouse_ortholog_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

strict_mouse_signature <- unique(ortholog_audit$mouse_mgi_symbol[ortholog_audit$strict_mapping_status == "retained_unique_mouse_ortholog_support_n_ge_2"])
if (length(strict_mouse_signature) == 0L) stop("No unique Human-to-Mouse orthologs supported by at least two sources were available from the Human MG0 top-50 signature.")
write.table(
  data.frame(mouse_mgi_symbol = strict_mouse_signature, stringsAsFactors = FALSE),
  file.path(table_dir, "Human_MG0_top50_strict_mouse_ortholog_signature.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

all_slice_intersection <- data.frame()
all_slice_cluster_summary <- data.frame()
all_slice_comparisons <- data.frame()

for (slice_id in names(slice_paths)) {
  slice_object <- readRDS(slice_paths[[slice_id]])
  DefaultAssay(slice_object) <- "RNA"
  if (!"CellType" %in% colnames(slice_object[[]])) stop(slice_id, " lacks CellType.")
  slice_object$plot_celltype <- factor(as.character(slice_object$CellType), levels = target_order)
  if (any(is.na(slice_object$plot_celltype))) stop(slice_id, " contains a CellType outside ST_MG0--ST_MG5.")
  if (!all(target_order %in% levels(droplevels(slice_object$plot_celltype)))) stop(slice_id, " lacks one or more ST_MG groups.")

  slice_RNA_genes <- rownames(slice_object[["RNA"]])
  slice_mouse_signature_genes <- intersect(strict_mouse_signature, slice_RNA_genes)
  if (length(slice_mouse_signature_genes) == 0L) stop(slice_id, " has no strict Human-to-Mouse ortholog signature genes in the native MERFISH RNA panel.")
  slice_intersection_table <- merge(
    ortholog_audit[, c("human_gene", "mouse_mgi_symbol", "strict_mapping_status")],
    data.frame(mouse_mgi_symbol = strict_mouse_signature, stringsAsFactors = FALSE),
    by = "mouse_mgi_symbol", all.y = TRUE, sort = FALSE
  )
  slice_intersection_table$slice <- slice_id
  slice_intersection_table$present_in_slice_RNA <- slice_intersection_table$mouse_mgi_symbol %in% slice_RNA_genes
  slice_intersection_table <- slice_intersection_table[, c("slice", "human_gene", "mouse_mgi_symbol", "strict_mapping_status", "present_in_slice_RNA")]
  all_slice_intersection <- rbind(all_slice_intersection, slice_intersection_table)

  slice_object <- AddModuleScore(
    object = slice_object,
    assay = "RNA",
    features = list(slice_mouse_signature_genes),
    ctrl = 10,
    name = "Human_MG0_top50_mouse_ortholog",
    seed = 20260916
  )
  score_column_from_AddModuleScore <- "Human_MG0_top50_mouse_ortholog1"
  if (!score_column_from_AddModuleScore %in% colnames(slice_object[[]])) stop(slice_id, " did not produce ", score_column_from_AddModuleScore, ".")
  score_column <- "Human_MG0_top50_mouse_ortholog_score"
  slice_object[[score_column]] <- slice_object[[score_column_from_AddModuleScore]][, 1]

  comparison_groups <- levels(slice_object$plot_celltype)
  my_comparisons <- lapply(comparison_groups[-1], function(x) c(x, comparison_groups[1]))
  slice_violin <- FeatureStatPlot(
    srt = slice_object,
    assay = "RNA",
    slot = "data",
    group.by = "plot_celltype",
    stat.by = score_column,
    add_box = FALSE,
    comparisons = my_comparisons,
    bg.by = "plot_celltype",
    palcolor = st_mg_colors_npg,
    bg_palcolor = st_mg_colors_npg,
    legend.position = "right",
    xlab = "",
    theme_args = list(base_size = 10)
  )
  slice_violin <- slice_violin & theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.margin = ggplot2::margin(t = 16, r = 4, b = 4, l = 4, unit = "pt")
  )
  slice_violin <- slice_violin + patchwork::plot_annotation(title = paste0(slice_id, ": Human MG0 top-50 strict mouse-ortholog module score"))
  ggsave(file.path(figure_dir, paste0(slice_id, "_Human_MG0_top50_strict_mouse_ortholog_violin.png")), slice_violin, width = 8.5, height = 6.2, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, paste0(slice_id, "_Human_MG0_top50_strict_mouse_ortholog_violin.pdf")), slice_violin, width = 8.5, height = 6.2, device = grDevices::cairo_pdf, bg = "white")

  score_values <- as.numeric(slice_object[[score_column]][, 1])
  celltype_values <- as.character(slice_object$plot_celltype)
  slice_cluster_summary <- do.call(rbind, lapply(comparison_groups, function(group_name) {
    group_scores <- score_values[celltype_values == group_name]
    data.frame(
      slice = slice_id,
      cluster = group_name,
      cells = length(group_scores),
      mean_Human_MG0_top50_mouse_ortholog_score = mean(group_scores),
      median_Human_MG0_top50_mouse_ortholog_score = median(group_scores),
      sd_Human_MG0_top50_mouse_ortholog_score = sd(group_scores),
      stringsAsFactors = FALSE
    )
  }))
  all_slice_cluster_summary <- rbind(all_slice_cluster_summary, slice_cluster_summary)

  slice_comparison_statistics <- do.call(rbind, lapply(comparison_groups[-1], function(comparison_cluster) {
    mg0_scores <- score_values[celltype_values == "ST_MG0"]
    other_scores <- score_values[celltype_values == comparison_cluster]
    test_result <- wilcox.test(other_scores, mg0_scores, alternative = "two.sided", exact = FALSE)
    data.frame(
      slice = slice_id,
      reference_cluster = "ST_MG0",
      comparison_cluster = comparison_cluster,
      ST_MG0_cells = length(mg0_scores),
      comparison_cells = length(other_scores),
      ST_MG0_mean_Human_MG0_top50_mouse_ortholog_score = mean(mg0_scores),
      comparison_mean_Human_MG0_top50_mouse_ortholog_score = mean(other_scores),
      wilcoxon_W = unname(test_result$statistic),
      p_value_unadjusted = test_result$p.value,
      stringsAsFactors = FALSE
    )
  }))
  slice_comparison_statistics$p_value_BH_within_slice <- p.adjust(slice_comparison_statistics$p_value_unadjusted, method = "BH")
  slice_comparison_statistics$plot_significance <- ifelse(slice_comparison_statistics$p_value_unadjusted <= 1e-4, "****", ifelse(slice_comparison_statistics$p_value_unadjusted <= 1e-3, "***", ifelse(slice_comparison_statistics$p_value_unadjusted <= 1e-2, "**", ifelse(slice_comparison_statistics$p_value_unadjusted <= 0.05, "*", "ns"))))
  all_slice_comparisons <- rbind(all_slice_comparisons, slice_comparison_statistics)

  slice_cell_audit <- data.frame(
    cell_barcode = colnames(slice_object),
    slice = slice_id,
    CellType = as.character(slice_object$CellType),
    plot_celltype = as.character(slice_object$plot_celltype),
    Human_MG0_top50_mouse_ortholog_score = score_values,
    stringsAsFactors = FALSE
  )
  write.table(slice_cell_audit, file.path(table_dir, paste0(slice_id, "_Human_MG0_top50_strict_mouse_ortholog_per_cell.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(slice_object, file.path(object_dir, paste0(slice_id, "_ST_MG0_MG5_Human_MG0_top50_strict_mouse_ortholog_scored.rds")), compress = "gzip")
}

write.table(all_slice_intersection, file.path(table_dir, "Human_MG0_top50_strict_mouse_ortholog_intersection_by_slice.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_slice_cluster_summary, file.path(table_dir, "Human_MG0_top50_strict_mouse_ortholog_summary_by_slice_and_cluster.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_slice_comparisons, file.path(table_dir, "Human_MG0_top50_strict_mouse_ortholog_Wilcoxon_ST_MG0_vs_other_clusters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parameter_record <- data.frame(
  item = c("reference_script", "human_marker_reference", "ortholog_source", "ortholog_rule", "spatial_species", "score_method", "score_assay", "control_genes", "score_column", "grouping", "reference_group", "violin_method", "test"),
  value = c(
    file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/scripts/spatial_module.r"),
    human_top50_tsv,
    paste0("babelgene ", packageVersion("babelgene"), "; bundled HGNC Comparison of Orthology Predictions (HCOP) data"),
    "Retain a non-empty mouse MGI symbol with support_n >= 2; require exactly one unique mouse candidate per human gene; then require exact match to the native MERFISH RNA feature.",
    "Mouse",
    "Seurat::AddModuleScore", "RNA", "10", score_column,
    "CellType (ST_MG0--ST_MG5)", "ST_MG0", "SCP::FeatureStatPlot", "two-sided Wilcoxon; raw p displayed as stars, BH p saved per slice"
  ),
  stringsAsFactors = FALSE
)
write.table(parameter_record, file.path(result_dir, "03_Human_MG0_top50_strict_mouse_ortholog_parameter_record.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "03_Human_MG0_top50_strict_mouse_ortholog_module_score_sessionInfo.txt"))
writeLines(as.character(Sys.time()), file.path(log_dir, "03_Human_MG0_top50_strict_mouse_ortholog_module_score_completed_at.txt"))
message("Human MG0 top-50 strict Human-to-Mouse-ortholog module-score branch completed.")
