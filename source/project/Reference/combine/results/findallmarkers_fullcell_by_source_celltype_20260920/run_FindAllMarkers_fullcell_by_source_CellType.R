#!/usr/bin/env Rscript

# Per-study FindAllMarkers outputs for the exact full-cell objects used in the
# CellType QC violin plots. Historical MG barcodes are matched to the final
# human reference, replacing their source CellType for analysis only:
# MG0 -> OLM; MG1-MG5 retained. Source objects are never saved or overwritten.

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
out_dir <- file.path(project_root, "combine/results/findallmarkers_fullcell_by_source_celltype_20260920")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

final_mg_path <- file.path(project_root, "combine/data/processed_data/mg_cca_meta_v2.qs")
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

marker_parameters <- list(
  assay = "RNA",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
mg_levels <- c("OLM", "MG1", "MG2", "MG3", "MG4", "MG5")

final_mg <- qread(final_mg_path)
final_meta <- final_mg[[]]
if (!all(c("Study", "CellType1") %in% colnames(final_meta))) {
  stop("The final MG reference lacks Study and/or CellType1 metadata.")
}
if (!all(final_meta$CellType1 %in% paste0("MG", 0:5))) {
  stop("The final MG reference contains CellType1 outside MG0-MG5.")
}

write_tsv <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

all_marker_list <- list()
mapping_audit_list <- list()
group_count_list <- list()
superseded_dir <- file.path(out_dir, "superseded_averageexpression_files")
old_expression_files <- list.files(
  out_dir,
  pattern = "_FindAllMarkers_average_normalized_expression_matrix\\.tsv$",
  full.names = TRUE
)
if (length(old_expression_files) > 0L) {
  dir.create(superseded_dir, recursive = TRUE, showWarnings = FALSE)
  file.rename(old_expression_files, file.path(superseded_dir, basename(old_expression_files)))
}

for (i in seq_len(nrow(study_manifest))) {
  study <- study_manifest$Study[[i]]
  historical_study <- study_manifest$historical_Study[[i]]
  source_path <- study_manifest$source_path[[i]]
  message("Processing ", study)

  source_object <- qread(source_path)
  source_meta <- source_object[[]]
  if (!all(c("CellType", "nCount_RNA", "nFeature_RNA") %in% colnames(source_meta))) {
    stop(study, " is missing CellType, nCount_RNA, and/or nFeature_RNA metadata.")
  }
  if (!"RNA" %in% Assays(source_object)) stop(study, " does not contain an RNA assay.")
  tryCatch(
    GetAssayData(source_object, assay = "RNA", layer = "data"),
    error = function(e) stop(study, " RNA normalized data are unavailable: ", conditionMessage(e))
  )

  historical_cells <- rownames(final_meta)[final_meta$Study == historical_study]
  if (length(historical_cells) == 0L) stop("No historical MG cells were found for ", study)
  if (!all(historical_cells %in% rownames(source_meta))) {
    missing_cells <- setdiff(historical_cells, rownames(source_meta))
    stop(study, " historical MG cells are absent from the source object (first: ", missing_cells[[1]], ")")
  }

  marker_group <- as.character(source_meta$CellType)
  marker_group[is.na(marker_group) | marker_group == ""] <- "Unannotated"
  historical_states <- as.character(final_meta[historical_cells, "CellType1", drop = TRUE])
  marker_group[match(historical_cells, rownames(source_meta))] <- ifelse(
    historical_states == "MG0", "OLM", historical_states
  )
  non_mg_levels <- sort(setdiff(unique(marker_group), mg_levels))
  marker_group <- factor(marker_group, levels = c(mg_levels[mg_levels %in% marker_group], non_mg_levels))
  names(marker_group) <- rownames(source_meta)
  marker_group <- marker_group[colnames(source_object)]
  source_object$marker_group <- marker_group
  DefaultAssay(source_object) <- "RNA"
  Idents(source_object) <- "marker_group"

  group_counts <- as.data.frame(table(marker_group), stringsAsFactors = FALSE)
  colnames(group_counts) <- c("marker_group", "n_cells")
  group_counts$Study <- study
  group_counts <- group_counts[, c("Study", "marker_group", "n_cells")]
  group_count_list[[study]] <- group_counts

  marker_file <- file.path(out_dir, paste0(study, "_FindAllMarkers_by_source_CellType.tsv"))
  if (file.exists(marker_file)) {
    markers <- read.delim(marker_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    markers <- FindAllMarkers(
      source_object,
      assay = marker_parameters$assay,
      test.use = marker_parameters$test.use,
      only.pos = marker_parameters$only.pos,
      min.pct = marker_parameters$min.pct,
      logfc.threshold = marker_parameters$logfc.threshold,
      verbose = FALSE
    )
    if (!"gene" %in% colnames(markers)) markers$gene <- rownames(markers)
    markers$Study <- study
    markers$analysis_group <- as.character(markers$cluster)
    markers <- markers[, c("Study", "analysis_group", setdiff(colnames(markers), c("Study", "analysis_group"))), drop = FALSE]
    write_tsv(markers, marker_file)
  }
  all_marker_list[[study]] <- markers

  marker_genes <- sort(unique(markers$gene))
  if (length(marker_genes) == 0L) {
    stop(study, " has no positive markers under the historical thresholds.")
  }
  fc_column <- intersect(c("avg_log2FC", "avg_logFC"), colnames(markers))
  if (length(fc_column) == 0L) stop(study, " FindAllMarkers output has no log-fold-change column.")
  fc_matrix <- matrix(NA_real_, nrow = length(marker_genes), ncol = length(levels(marker_group)),
                      dimnames = list(marker_genes, levels(marker_group)))
  for (group in levels(marker_group)) {
    group_markers <- markers[markers$analysis_group == group, c("gene", fc_column[[1]]), drop = FALSE]
    if (nrow(group_markers) == 0L) next
    group_fc <- tapply(group_markers[[fc_column[[1]]]], group_markers$gene, max)
    fc_matrix[names(group_fc), group] <- as.numeric(group_fc)
  }
  fc_output <- data.frame(gene = rownames(fc_matrix), as.data.frame(fc_matrix, check.names = FALSE), check.names = FALSE)
  write_tsv(fc_output, file.path(out_dir, paste0(study, "_FindAllMarkers_", fc_column[[1]], "_matrix.tsv")))

  source_data <- GetAssayData(source_object, assay = "RNA", layer = "data")
  available_marker_genes <- intersect(marker_genes, rownames(source_data))
  cells_by_group <- split(colnames(source_object), as.character(marker_group))
  average_expression <- do.call(cbind, lapply(levels(marker_group), function(group) {
    Matrix::rowMeans(source_data[available_marker_genes, cells_by_group[[group]], drop = FALSE])
  }))
  rownames(average_expression) <- available_marker_genes
  colnames(average_expression) <- levels(marker_group)
  expression_output <- data.frame(
    gene = rownames(average_expression),
    as.data.frame(as.matrix(average_expression), check.names = FALSE),
    check.names = FALSE
  )
  write_tsv(expression_output, file.path(out_dir, paste0(study, "_FindAllMarkers_mean_RNA_data_layer_expression_matrix.tsv")))

  mapping_audit_list[[study]] <- data.frame(
    Study = study,
    historical_CellType1 = paste0("MG", 0:5),
    analysis_group = mg_levels,
    n_cells = as.integer(table(factor(historical_states, levels = paste0("MG", 0:5)))),
    stringsAsFactors = FALSE
  )
  rm(source_object, source_meta, markers, fc_matrix, average_expression)
  gc(verbose = FALSE)
}

all_markers <- do.call(rbind, all_marker_list)
write_tsv(all_markers, file.path(out_dir, "All_studies_FindAllMarkers_by_source_CellType.tsv"))
write_tsv(do.call(rbind, group_count_list), file.path(out_dir, "analysis_group_cell_counts.tsv"))
write_tsv(do.call(rbind, mapping_audit_list), file.path(out_dir, "historical_MG0_MG5_mapping_audit.tsv"))

writeLines(c(
  "FindAllMarkers provenance",
  "",
  "Input: six original full-cell Seurat objects used for the CellType QC violin plots.",
  "Group labels: original CellType, except exact historical MG barcodes matched to mg_cca_meta_v2.qs.",
  "Historical mapping: MG0 is written as OLM; MG1-MG5 retain their historical names.",
  "FindAllMarkers parameters: RNA assay; default Wilcoxon; only.pos=TRUE; min.pct=0.25; logfc.threshold=0.25.",
  "Each study is analyzed separately. No source Seurat object is overwritten.",
  "Mean RNA data-layer expression matrices contain the union of positive FindAllMarkers genes for that study.",
  "Means are calculated directly from the stored RNA/data layer without Seurat exponentiation or back-transformation."
), file.path(out_dir, "README.txt"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Finished all studies. Output directory: ", out_dir)
