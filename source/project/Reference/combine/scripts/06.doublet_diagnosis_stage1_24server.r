#!/usr/bin/env Rscript

# MG0--MG5 doublet diagnosis, stage 1: raw count QC, diagnostic Seurat UMAP,
# DoubletFinder/scDblFinder calls, and Scrublet input export.

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(BiocParallel)
  library(Matrix)
  library(ggplot2)
  library(qs)
})
set.seed(20260902)

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
input_file <- file.path(project_root, "combine/data/processed_data/mg_cca_meta_v2.qs")
result_root <- file.path(project_root, "combine/results/MG0_MG5_doublet_diagnosis")
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
checkpoint_dir <- file.path(project_root, "combine/data/processed_data/MG0_MG5_doublet_diagnosis")
scrublet_input_dir <- file.path(result_root, "scrublet_input")
runtime_dir <- file.path(result_root, "runtime")
for (d in c(result_root, figure_dir, table_dir, checkpoint_dir, scrublet_input_dir, runtime_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = grDevices::pdf, bg = "white", limitsize = FALSE)
}

# Each study gets empirical conservative QC thresholds. Low-library boundaries are
# its 1st percentile with a minimal floor of 200; mt is its 99th percentile
# constrained to 10--20%. No upper UMI/gene cutoff is applied before doublet calling.
make_qc_thresholds <- function(md) c(
  min_nCount_RNA = max(200, floor(as.numeric(quantile(md$nCount_RNA, 0.01, na.rm = TRUE)))),
  min_nFeature_RNA = max(200, floor(as.numeric(quantile(md$nFeature_RNA, 0.01, na.rm = TRUE)))),
  max_percent_mt = min(20, max(10, ceiling(as.numeric(quantile(md$percent.mt, 0.99, na.rm = TRUE)))))
)

run_doubletfinder <- function(seu, dims_use, expected_rate) {
  n_exp_poi <- max(1L, round(expected_rate * ncol(seu)))
  message("DoubletFinder: parameter sweep")
  sweep <- DoubletFinder::paramSweep(seu, PCs = dims_use, sct = FALSE)
  message("DoubletFinder: summarizing sweep")
  sweep_stats <- DoubletFinder::summarizeSweep(sweep, GT = FALSE)
  message("DoubletFinder: selecting pK")
  pk_table <- DoubletFinder::find.pK(sweep_stats)
  pk_values <- suppressWarnings(as.numeric(as.character(pk_table$pK)))
  best_pk <- pk_values[which.max(pk_table$BCmetric)]
  if (!is.finite(best_pk)) stop("DoubletFinder pK optimisation returned no numeric pK.")
  message("DoubletFinder: estimating homotypic fraction")
  homotypic_fraction <- DoubletFinder::modelHomotypic(as.character(seu$diagnostic_cluster))
  n_exp_adj <- max(1L, round(n_exp_poi * (1 - homotypic_fraction)))
  message("DoubletFinder: final classification")
  seu <- DoubletFinder::doubletFinder(seu, PCs = dims_use, pN = 0.25, pK = best_pk,
                                      nExp = n_exp_adj, reuse.pANN = NULL, sct = FALSE)
  md_names <- colnames(seu[[]])
  class_col <- grep("^DF.classifications_", md_names, value = TRUE)
  score_col <- grep("^pANN_", md_names, value = TRUE)
  if (length(class_col) != 1L || length(score_col) != 1L) stop("Cannot uniquely identify DoubletFinder output columns.")
  seu$DoubletFinder_score <- seu[[score_col]][, 1]
  seu$DoubletFinder_call <- factor(ifelse(grepl("Doublet", seu[[class_col]][, 1]), "Doublet", "Singlet"), levels = c("Singlet", "Doublet"))
  list(seu = seu, diagnostics = data.frame(pN = 0.25, pK = best_pk, nExp_raw = n_exp_poi,
       nExp_homotypic_adjusted = n_exp_adj, homotypic_fraction = homotypic_fraction))
}

run_scdblfinder <- function(counts, metadata, clusters, expected_rate) {
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(metadata)
  sce <- scDblFinder::scDblFinder(sce, clusters = as.character(clusters), dbr = expected_rate,
                                  BPPARAM = BiocParallel::SerialParam(progressbar = FALSE), verbose = TRUE)
  list(score = as.numeric(SummarizedExperiment::colData(sce)$scDblFinder.score),
       call = as.character(SummarizedExperiment::colData(sce)$scDblFinder.class),
       threshold = as.numeric(S4Vectors::metadata(sce)$scDblFinder.threshold))
}

message("Loading MG0--MG5 object and recalculating QC directly from integer RNA counts ...")
mg <- qs::qread(input_file)
counts_all <- mg[["RNA"]]@counts
if (!inherits(counts_all, "dgCMatrix")) counts_all <- as(counts_all, "dgCMatrix")
if (any(abs(counts_all@x - round(counts_all@x)) > 1e-8)) stop("RNA counts are not integer-valued.")
meta_all <- mg[[]]
mt_genes <- grep("^MT-", rownames(counts_all), value = TRUE)
ribo_genes <- grep("^RP[SL]", rownames(counts_all), value = TRUE)
meta_all$nCount_RNA <- Matrix::colSums(counts_all)
meta_all$nFeature_RNA <- Matrix::colSums(counts_all > 0)
meta_all$percent.mt <- if (length(mt_genes)) Matrix::colSums(counts_all[mt_genes, , drop = FALSE]) / pmax(meta_all$nCount_RNA, 1) * 100 else 0
meta_all$percent.ribo <- if (length(ribo_genes)) Matrix::colSums(counts_all[ribo_genes, , drop = FALSE]) / pmax(meta_all$nCount_RNA, 1) * 100 else 0
missing_cols <- setdiff(c("Study", "Sample", "CellType1"), colnames(meta_all))
if (length(missing_cols)) stop("Missing metadata: ", paste(missing_cols, collapse = ", "))
study_order <- c("Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
                 "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024")
studies <- intersect(study_order, unique(as.character(meta_all$Study)))
if (length(studies) != 6L) stop("Expected six studies; found ", paste(studies, collapse = ", "))

threshold_rows <- list()
qc_cell_rows <- list()
scrublet_manifest_rows <- list()
df_parameter_rows <- list()
stage1_metadata_rows <- list()

for (study in studies) {
  message("\n===== ", study, " =====")
  study_cells <- rownames(meta_all)[as.character(meta_all$Study) == study]
  study_counts <- counts_all[, study_cells, drop = FALSE]
  study_meta <- meta_all[study_cells, , drop = FALSE]
  thresholds <- make_qc_thresholds(study_meta)
  study_meta$passes_basic_QC <- with(study_meta,
    nCount_RNA >= thresholds[["min_nCount_RNA"]] & nFeature_RNA >= thresholds[["min_nFeature_RNA"]] & percent.mt <= thresholds[["max_percent_mt"]])
  threshold_rows[[study]] <- data.frame(Study = study, n_before = nrow(study_meta),
    n_after = sum(study_meta$passes_basic_QC), n_removed = sum(!study_meta$passes_basic_QC),
    min_nCount_RNA = thresholds[["min_nCount_RNA"]], min_nFeature_RNA = thresholds[["min_nFeature_RNA"]],
    max_percent_mt = thresholds[["max_percent_mt"]])
  qc_cell_rows[[study]] <- cbind(Cell = rownames(study_meta), Study = study, study_meta[, c(
    "Sample", "CellType1", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "passes_basic_QC"), drop = FALSE])
  keep <- which(study_meta$passes_basic_QC)
  if (length(keep) < 500L) stop(study, " has fewer than 500 cells after basic QC.")
  study_counts <- study_counts[, keep, drop = FALSE]
  study_meta <- study_meta[keep, , drop = FALSE]
  expected_rate <- min(0.20, 0.008 * ncol(study_counts) / 1000)
  study_meta$expected_doublet_rate <- expected_rate

  seu <- Seurat::CreateSeuratObject(counts = study_counts, project = study, meta.data = study_meta)
  seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  n_hvg <- min(2000L, nrow(seu) - 1L)
  seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = n_hvg, verbose = FALSE)
  hvg <- Seurat::VariableFeatures(seu)
  seu <- Seurat::ScaleData(seu, features = hvg, verbose = FALSE)
  n_pcs <- min(30L, length(hvg) - 1L, ncol(seu) - 1L)
  seu <- Seurat::RunPCA(seu, features = hvg, npcs = n_pcs, verbose = FALSE, seed.use = 20260902)
  dims_use <- seq_len(min(30L, ncol(Seurat::Embeddings(seu, "pca"))))
  seu <- Seurat::FindNeighbors(seu, dims = dims_use, verbose = FALSE)
  seu <- Seurat::FindClusters(seu, resolution = 0.4, verbose = FALSE)
  seu$diagnostic_cluster <- as.character(seu$seurat_clusters)
  seu <- Seurat::RunUMAP(seu, dims = dims_use, seed.use = 20260902, verbose = FALSE)
  study_safe <- safe_name(study)
  qs::qsave(seu, file.path(checkpoint_dir, paste0("MG_doublet_pre_DF_", study_safe, ".qs")), preset = "high")

  df_out <- tryCatch(run_doubletfinder(seu, dims_use, expected_rate), error = function(e) e)
  if (inherits(df_out, "error")) stop(study, ": DoubletFinder failed: ", conditionMessage(df_out))
  seu <- df_out$seu
  df_parameter_rows[[study]] <- cbind(Study = study, expected_doublet_rate = expected_rate, df_out$diagnostics)
  sc_out <- tryCatch(run_scdblfinder(study_counts, seu[[]], seu$diagnostic_cluster, expected_rate), error = function(e) e)
  if (inherits(sc_out, "error")) stop(study, ": scDblFinder failed: ", conditionMessage(sc_out))
  seu$scDblFinder_score <- sc_out$score
  seu$scDblFinder_call <- factor(ifelse(tolower(sc_out$call) == "doublet", "Doublet", "Singlet"), levels = c("Singlet", "Doublet"))
  seu$scDblFinder_threshold <- sc_out$threshold

  study_scrublet_dir <- file.path(scrublet_input_dir, study_safe)
  dir.create(study_scrublet_dir, recursive = TRUE, showWarnings = FALSE)
  Matrix::writeMM(study_counts, file.path(study_scrublet_dir, "counts_genes_by_cells.mtx"))
  utils::write.table(data.frame(cell_id = colnames(study_counts)), file.path(study_scrublet_dir, "cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(gene_id = rownames(study_counts)), file.path(study_scrublet_dir, "genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  scrublet_manifest_rows[[study]] <- data.frame(Study = study, safe_study = study_safe, n_cells = ncol(study_counts),
    n_genes = nrow(study_counts), expected_doublet_rate = expected_rate)
  qs::qsave(seu, file.path(checkpoint_dir, paste0("MG_doublet_stage1_", study_safe, ".qs")), preset = "high")
  current_meta <- seu[[]]
  current_meta$Cell <- rownames(current_meta)
  umap <- Seurat::Embeddings(seu, "umap")
  current_meta$UMAP_1 <- umap[rownames(current_meta), 1]
  current_meta$UMAP_2 <- umap[rownames(current_meta), 2]
  stage1_metadata_rows[[study]] <- current_meta[, c("Cell", "Study", "Sample", "CellType1", "nCount_RNA", "nFeature_RNA",
    "percent.mt", "percent.ribo", "diagnostic_cluster", "expected_doublet_rate", "DoubletFinder_score", "DoubletFinder_call",
    "scDblFinder_score", "scDblFinder_call", "UMAP_1", "UMAP_2"), drop = FALSE]
  rm(seu, study_counts, study_meta, current_meta, umap)
  gc(verbose = FALSE)
}

threshold_table <- do.call(rbind, threshold_rows)
qc_cell_table <- do.call(rbind, qc_cell_rows)
scrublet_manifest <- do.call(rbind, scrublet_manifest_rows)
df_parameter_table <- do.call(rbind, df_parameter_rows)
stage1_metadata <- do.call(rbind, stage1_metadata_rows)
utils::write.csv(threshold_table, file.path(table_dir, "basic_QC_thresholds_and_retention.csv"), row.names = FALSE)
utils::write.csv(qc_cell_table, file.path(table_dir, "basic_QC_per_cell.csv"), row.names = FALSE)
utils::write.csv(scrublet_manifest, file.path(scrublet_input_dir, "scrublet_manifest.csv"), row.names = FALSE)
utils::write.csv(df_parameter_table, file.path(table_dir, "DoubletFinder_parameter_selection.csv"), row.names = FALSE)
utils::write.csv(stage1_metadata, file.path(table_dir, "stage1_cell_metadata_before_scrublet.csv"), row.names = FALSE)

qc_long <- rbind(data.frame(Study = qc_cell_table$Study, metric = "nCount_RNA", value = qc_cell_table$nCount_RNA),
  data.frame(Study = qc_cell_table$Study, metric = "nFeature_RNA", value = qc_cell_table$nFeature_RNA),
  data.frame(Study = qc_cell_table$Study, metric = "percent.mt", value = qc_cell_table$percent.mt))
threshold_long <- rbind(data.frame(Study = threshold_table$Study, metric = "nCount_RNA", threshold = threshold_table$min_nCount_RNA),
  data.frame(Study = threshold_table$Study, metric = "nFeature_RNA", threshold = threshold_table$min_nFeature_RNA),
  data.frame(Study = threshold_table$Study, metric = "percent.mt", threshold = threshold_table$max_percent_mt))
p_qc <- ggplot(qc_long, aes(x = value)) + geom_histogram(bins = 60, fill = "#9ecae1", colour = "white", linewidth = 0.1) +
  geom_vline(data = threshold_long, aes(xintercept = threshold), colour = "#d7301f", linewidth = 0.45) +
  facet_grid(metric ~ Study, scales = "free") +
  labs(title = "Data-driven basic QC thresholds before doublet diagnosis",
       subtitle = "Red line: lower nCount/nFeature boundary or upper mitochondrial boundary", x = NULL, y = "Cell count") +
  theme_bw(base_size = 10) + theme(strip.text = element_text(face = "bold"))
save_plot_both(p_qc, file.path(figure_dir, "Basic_QC_thresholds_six_studies"), 15, 8)
# Keep the diagnostic QC figure visible in RStudio's Plots history as well as
# saving its reproducible PNG and PDF files.
print(p_qc)

writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "stage1_R_sessionInfo.txt"))
writeLines(c(
  "Scope: six study-specific microglia subsets retained in mg_cca_meta_v2.qs.",
  "Input: integer RNA@counts only; existing normalized, SCT, integrated and UMAP slots are not used.",
  "Basic QC: empirical 1st percentile lower thresholds with floors of 200 UMIs/200 genes; mitochondrial upper threshold is the 99th percentile constrained to 10--20%.",
  "No upper UMI/gene filter was applied, so potential doublets remain available to the callers.",
  "Expected doublet prior for all methods: min(0.20, 0.008 * recovered cells / 1000).",
  "Limitation: post-MG-subset calls primarily find within-MG doublets; heterotypic MG--non-MG doublets may be underestimated."
), file.path(runtime_dir, "stage1_method_manifest.txt"))
message("Stage 1 complete. Run Scrublet stage 2, then stage 3 plotting.")
