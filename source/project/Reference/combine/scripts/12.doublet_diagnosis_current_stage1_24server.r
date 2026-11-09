#!/usr/bin/env Rscript

# Stage 3, part 1: current reanalysis doublet diagnosis.
#
# Input is the Stage-2 object saved after repeated basic QC but BEFORE the
# Stage-2 DoubletFinder removal.  It retains raw integer RNA counts and the
# uncorrected diagnostic UMAP.  The existing DoubletFinder calls are therefore
# the calls for exactly the same input cells used below for scDblFinder and
# Scrublet.  No cell is removed by this diagnostic workflow.

rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
invisible(gc())

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(qs)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(scDblFinder)
  library(BiocParallel)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260907)

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
stage2_root <- file.path(project_root, "combine/results/full_count_reanalysis/02_microglia_integration")
result_root <- file.path(project_root, "combine/results/full_count_reanalysis/03_current_MG0_doublet_diagnosis")
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
object_dir <- file.path(result_root, "object")
runtime_dir <- file.path(result_root, "runtime")
scrublet_input_dir <- file.path(result_root, "scrublet_input")
for (d in c(figure_dir, table_dir, object_dir, runtime_dir, scrublet_input_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

diagnostic_file <- file.path(stage2_root, "object", "selected_Microglia_QC_before_doublet_removal.qs")
reference_file <- file.path(stage2_root, "object", "integrated_Microglia_processed.qs")
mg0_selection_file <- file.path(project_root, "combine/results/full_count_reanalysis/manual_selections/integrated_Microglia_MG0_cluster_selection_TEMPLATE.csv")
df_parameters_file <- file.path(stage2_root, "table", "selected_Microglia_DoubletFinder_parameters.csv")

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
         dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
         device = grDevices::pdf, bg = "white", limitsize = FALSE)
}
get_counts <- function(obj) {
  x <- tryCatch(SeuratObject::LayerData(obj[["RNA"]], layer = "counts"), error = function(e) NULL)
  if (is.null(x) || !length(x)) x <- tryCatch(Seurat::GetAssayData(obj, assay = "RNA", layer = "counts"), error = function(e) NULL)
  if (is.null(x) || !length(x)) x <- tryCatch(Seurat::GetAssayData(obj, assay = "RNA", slot = "counts"), error = function(e) NULL)
  if (is.null(x) || !length(x)) stop("RNA count layer is unavailable.")
  if (!inherits(x, "dgCMatrix")) x <- as(x, "dgCMatrix")
  if (length(x@x) && any(abs(x@x - round(x@x)) > 1e-8)) stop("RNA counts are not integer-valued.")
  x
}

if (any(!file.exists(c(diagnostic_file, reference_file, mg0_selection_file, df_parameters_file)))) {
  stop("One or more required Stage-2 inputs are missing.")
}

selection <- read.csv(mg0_selection_file, check.names = FALSE, stringsAsFactors = FALSE)
selection$selected <- tolower(trimws(as.character(selection$select_as_MG0))) %in% c("true", "t", "1", "yes", "y")
selected_clusters <- as.character(selection$integrated_cluster[selection$selected])
if (length(selected_clusters) != 1L || selected_clusters != "2") {
  stop("This stage requires the recorded manual selection: exactly integrated cluster 2 as MG0.")
}

message("Loading QC-passed pre-removal diagnostic object and final singlet reference ...")
diagnostic <- qs::qread(diagnostic_file)
reference <- qs::qread(reference_file)
counts <- get_counts(diagnostic)
if (ncol(diagnostic) != ncol(counts)) stop("Diagnostic object and count matrix disagree.")
if (!all(c("DoubletFinder_score", "DoubletFinder_call", "DF_partition", "diagnostic_cluster") %in% colnames(diagnostic[[]]))) {
  stop("The Stage-2 diagnostic object is missing required DoubletFinder/cluster metadata.")
}
if (!"integrated_cluster" %in% colnames(reference[[]])) stop("Final reference has no integrated_cluster metadata.")

# The manual cluster-2 definition was made on post-DoubletFinder singlets.
# Anchor transfer extends that definition to the pre-removal diagnostic cells,
# including cells called doublets, without changing the manual reference label.
transfer_features <- intersect(Seurat::VariableFeatures(reference), rownames(diagnostic))
if (length(transfer_features) < 1000L) stop("Fewer than 1,000 shared reference variable features for label transfer.")
reference <- Seurat::NormalizeData(reference, verbose = FALSE)
diagnostic <- Seurat::NormalizeData(diagnostic, verbose = FALSE)
transfer_dims <- seq_len(min(30L, ncol(Seurat::Embeddings(reference, "pca"))))
anchors <- Seurat::FindTransferAnchors(
  reference = reference, query = diagnostic, normalization.method = "LogNormalize",
  reference.reduction = "pca", features = transfer_features, dims = transfer_dims,
  verbose = FALSE
)
predictions <- Seurat::TransferData(
  anchorset = anchors, refdata = as.character(reference$integrated_cluster),
  dims = transfer_dims, verbose = FALSE
)
diagnostic$transferred_integrated_cluster <- as.character(predictions$predicted.id)
score_column <- paste0("prediction.score.", selected_clusters)
diagnostic$MG0_transfer_score <- if (score_column %in% colnames(predictions)) {
  as.numeric(predictions[[score_column]])
} else {
  NA_real_
}
diagnostic$MG0_status <- factor(
  ifelse(diagnostic$transferred_integrated_cluster == selected_clusters,
         "MG0 (manual cluster 2 transfer)", "Other predicted cluster"),
  levels = c("Other predicted cluster", "MG0 (manual cluster 2 transfer)")
)

df_parameters <- read.csv(df_parameters_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("DF_partition", "expected_doublet_rate") %in% colnames(df_parameters))) {
  stop("DoubletFinder parameter table is incomplete.")
}
partition_cells <- split(colnames(diagnostic), as.character(diagnostic$DF_partition))
partition_cells <- partition_cells[order(names(partition_cells))]
if (any(lengths(partition_cells) < 300L)) stop("A saved DoubletFinder partition has fewer than 300 cells.")
partition_rate <- setNames(df_parameters$expected_doublet_rate, as.character(df_parameters$DF_partition))
if (any(!names(partition_cells) %in% names(partition_rate))) stop("A DoubletFinder partition has no recorded expected rate.")

run_scdblfinder_partition <- function(cells, partition_name, expected_rate) {
  partition_counts <- counts[, cells, drop = FALSE]
  partition_counts <- partition_counts[Matrix::rowSums(partition_counts) > 0, , drop = FALSE]
  partition_md <- diagnostic[[]][cells, , drop = FALSE]
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = partition_counts))
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(partition_md)
  sce <- scDblFinder::scDblFinder(
    sce, clusters = as.character(partition_md$diagnostic_cluster), dbr = expected_rate,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE), verbose = TRUE
  )
  threshold <- tryCatch(as.numeric(S4Vectors::metadata(sce)$scDblFinder.threshold), error = function(e) NA_real_)
  data.frame(
    Cell = colnames(sce), DF_partition = partition_name,
    scDblFinder_score = as.numeric(SummarizedExperiment::colData(sce)$scDblFinder.score),
    scDblFinder_call = ifelse(tolower(as.character(SummarizedExperiment::colData(sce)$scDblFinder.class)) == "doublet", "Doublet", "Singlet"),
    scDblFinder_threshold = threshold,
    expected_doublet_rate = expected_rate,
    n_partition_cells = length(cells), n_detected_features = nrow(partition_counts),
    stringsAsFactors = FALSE
  )
}

message("Running scDblFinder in the same technical partitions used by DoubletFinder ...")
sc_checkpoint_dir <- file.path(runtime_dir, "scDblFinder_partitions")
dir.create(sc_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
sc_calls_by_partition <- lapply(names(partition_cells), function(partition_name) {
  checkpoint <- file.path(sc_checkpoint_dir, paste0(safe_name(partition_name), "_scDblFinder_calls.csv"))
  if (file.exists(checkpoint)) return(read.csv(checkpoint, stringsAsFactors = FALSE, check.names = FALSE))
  result <- run_scdblfinder_partition(partition_cells[[partition_name]], partition_name, partition_rate[[partition_name]])
  write.csv(result, checkpoint, row.names = FALSE)
  result
})
sc_calls <- do.call(rbind, sc_calls_by_partition)
if (anyDuplicated(sc_calls$Cell) || !setequal(sc_calls$Cell, colnames(diagnostic))) stop("scDblFinder per-cell calls are incomplete or duplicated.")
sc_calls <- sc_calls[match(colnames(diagnostic), sc_calls$Cell), , drop = FALSE]
diagnostic$scDblFinder_score <- sc_calls$scDblFinder_score
diagnostic$scDblFinder_call <- factor(sc_calls$scDblFinder_call, levels = c("Singlet", "Doublet"))
diagnostic$scDblFinder_threshold <- sc_calls$scDblFinder_threshold

message("Exporting sparse raw-count Scrublet inputs in the same technical partitions ...")
scrublet_manifest <- do.call(rbind, lapply(names(partition_cells), function(partition_name) {
  cells <- partition_cells[[partition_name]]
  part_dir <- file.path(scrublet_input_dir, safe_name(partition_name))
  dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
  part_counts <- counts[, cells, drop = FALSE]
  part_counts <- part_counts[Matrix::rowSums(part_counts) > 0, , drop = FALSE]
  Matrix::writeMM(part_counts, file.path(part_dir, "counts_genes_by_cells.mtx"))
  write.table(data.frame(cell_id = colnames(part_counts)), file.path(part_dir, "cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(gene_id = rownames(part_counts)), file.path(part_dir, "genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  data.frame(
    DF_partition = partition_name, safe_partition = safe_name(partition_name),
    Study = paste(sort(unique(as.character(diagnostic$Study[cells]))), collapse = ";"),
    n_cells = ncol(part_counts), n_genes = nrow(part_counts),
    expected_doublet_rate = partition_rate[[partition_name]], stringsAsFactors = FALSE
  )
}))
write.csv(scrublet_manifest, file.path(scrublet_input_dir, "scrublet_manifest.csv"), row.names = FALSE)

umap <- Seurat::Embeddings(diagnostic, "umap_diagnostic")
if (ncol(umap) < 2L) stop("Diagnostic UMAP coordinates are unavailable.")
meta <- diagnostic[[]]
meta$Cell <- rownames(meta)
meta$UMAP_1 <- umap[meta$Cell, 1]
meta$UMAP_2 <- umap[meta$Cell, 2]
write.csv(meta, file.path(table_dir, "per_cell_current_MG0_three_method_partial_calls.csv"), row.names = FALSE)
write.csv(sc_calls, file.path(table_dir, "scDblFinder_per_cell_calls.csv"), row.names = FALSE)
write.csv(scrublet_manifest, file.path(table_dir, "Scrublet_input_manifest.csv"), row.names = FALSE)

# Current labels for all six diagnostic maps are based on transferred final
# cluster labels.  This preserves the user's manual cluster-2 decision while
# permitting doublet calls to be displayed before any method removes a cell.
p_cluster <- Seurat::DimPlot(diagnostic, reduction = "umap_diagnostic", group.by = "diagnostic_cluster", label = TRUE, repel = TRUE) +
  ggtitle("QC-passed selected Microglia: diagnostic UMAP clusters")
p_mg0 <- Seurat::DimPlot(diagnostic, reduction = "umap_diagnostic", group.by = "MG0_status", cols = c("Other predicted cluster" = "#d9d9d9", "MG0 (manual cluster 2 transfer)" = "#c51b7d")) +
  ggtitle("QC-passed selected Microglia: MG0 transferred from manual cluster 2")
p_df <- Seurat::DimPlot(diagnostic, reduction = "umap_diagnostic", group.by = "DoubletFinder_call", cols = c(Singlet = "#bdbdbd", Doublet = "#d7301f")) +
  ggtitle("QC-passed selected Microglia: DoubletFinder call")
for (item in list(Diagnostic_UMAP_clusters = p_cluster, Diagnostic_UMAP_MG0_transfer = p_mg0, Diagnostic_UMAP_DoubletFinder = p_df)) {
  plot_name <- names(item)
  save_plot_both(item[[1]], file.path(figure_dir, plot_name), 10, 7)
  print(item[[1]])
}

qs::qsave(diagnostic, file.path(object_dir, "current_MG0_QC_before_doublet_removal_three_method_partial.qs"), preset = "high")
writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "stage1_R_sessionInfo.txt"))
writeLines(c(
  "Scope: the six manually selected Microglia populations after repeated basic QC and before Stage-2 DoubletFinder removal (16,931 cells).",
  "Raw counts: RNA counts are sparse integer values. Existing normalized values, old source UMAPs and old CellType1 labels are not used.",
  "DoubletFinder: calls were computed during Stage 2 on this exact QC-passed object, within recorded technical partitions; they are reused rather than recalculated.",
  "scDblFinder: newly run from the same raw counts and technical partitions as the recorded DoubletFinder calls.",
  "Scrublet: sparse count matrices are exported per identical technical partition for the foreground RStudio terminal run in stage 2.",
  "MG0: manual final integrated cluster 2 is transferred to pre-removal cells by Seurat anchor transfer. transferred_integrated_cluster and MG0_transfer_score are retained per cell.",
  "No doublet call deletes or relabels a cell in this diagnostic stage."
), file.path(runtime_dir, "stage1_method_manifest.txt"))
message("Current MG0 doublet diagnosis stage 1 complete. Run 12.doublet_diagnosis_current_stage2_scrublet_24server.py in the RStudio Terminal.")
