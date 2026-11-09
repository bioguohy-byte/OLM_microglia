#!/usr/bin/env Rscript

# Artificial-doublet visualisation, part 3 of 3.
#
# For each study, this script creates a standalone Seurat object containing
# real pre-removal MG00-diagnostic cells plus three method-specific artificial
# profile sets.  It performs only LogNormalize/HVG/Scale/PCA/UMAP for visual
# inspection.  Artificial profiles are never batch-corrected with real cells
# and never used to reassign a real cell's doublet call or identity.

rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
invisible(gc())

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(qs)
  library(ggplot2)
})

set.seed(20260908)

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
combine_root <- file.path(project_root, "combine")
current_object_file <- file.path(
  combine_root,
  "results/full_count_reanalysis/03_current_MG00_doublet_diagnosis/object/current_MG00_QC_before_doublet_removal_all_three_calls.qs"
)
result_root <- file.path(combine_root, "results/artificial_doublet_visualization")
input_dir <- file.path(result_root, "input")
scrublet_dir <- file.path(result_root, "scrublet_output")
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
object_dir <- file.path(result_root, "object")
runtime_dir <- file.path(result_root, "runtime")
for (directory in c(figure_dir, table_dir, object_dir, runtime_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

study_order <- c(
  "Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
  "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"
)
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width = 8, height = 6) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
    dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
    device = grDevices::pdf, bg = "white", limitsize = FALSE)
}
get_counts <- function(object) {
  counts <- tryCatch(SeuratObject::LayerData(object[["RNA"]], layer = "counts"), error = function(e) NULL)
  if (is.null(counts) || !length(counts)) counts <- tryCatch(Seurat::GetAssayData(object, assay = "RNA", layer = "counts"), error = function(e) NULL)
  if (is.null(counts) || !length(counts)) counts <- tryCatch(Seurat::GetAssayData(object, assay = "RNA", slot = "counts"), error = function(e) NULL)
  if (is.null(counts) || !length(counts)) stop("The current MG00 object has no RNA count layer.")
  as(counts, "dgCMatrix")
}
read_sparse_bundle <- function(folder, matrix_file, cell_file, gene_file) {
  matrix_file <- file.path(folder, matrix_file)
  cell_file <- file.path(folder, cell_file)
  gene_file <- file.path(folder, gene_file)
  if (any(!file.exists(c(matrix_file, cell_file, gene_file)))) {
    stop("Missing sparse artificial-profile bundle in ", folder)
  }
  counts <- as(Matrix::readMM(matrix_file), "dgCMatrix")
  cells <- read.delim(cell_file, stringsAsFactors = FALSE, check.names = FALSE)$cell_id
  genes <- read.delim(gene_file, stringsAsFactors = FALSE, check.names = FALSE)$gene_id
  if (!identical(dim(counts), c(length(genes), length(cells)))) stop("Sparse matrix dimensions disagree with manifests in ", folder)
  rownames(counts) <- genes
  colnames(counts) <- cells
  counts
}
align_to_feature_union <- function(counts, feature_union) {
  if (identical(rownames(counts), feature_union)) return(as(counts, "dgCMatrix"))
  feature_index <- match(rownames(counts), feature_union)
  if (anyNA(feature_index) || anyDuplicated(feature_index)) stop("Artificial profile genes cannot be aligned to the real-cell feature union.")
  triplet <- as(counts, "dgTMatrix")
  output <- Matrix::sparseMatrix(
    i = feature_index[triplet@i + 1L], j = triplet@j + 1L, x = triplet@x,
    dims = c(length(feature_union), ncol(counts)),
    dimnames = list(feature_union, colnames(counts))
  )
  as(output, "dgCMatrix")
}
plot_base <- function(dat, title, subtitle = NULL) {
  ggplot(dat, aes(Artificial_UMAP_1, Artificial_UMAP_2)) +
    coord_equal() +
    labs(title = title, subtitle = subtitle, x = "Artificial-doublet UMAP 1", y = "Artificial-doublet UMAP 2") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      plot.margin = margin(t = 12, r = 14, b = 10, l = 14),
      legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 8)
    )
}
standardise_parent_provenance <- function(dat, method, parent_information) {
  if (!"parent_1" %in% colnames(dat)) dat$parent_1 <- NA_character_
  if (!"parent_2" %in% colnames(dat)) dat$parent_2 <- NA_character_
  if (!"scDblFinder_origin_cluster" %in% colnames(dat)) dat$scDblFinder_origin_cluster <- NA_character_
  dat$method <- method
  dat$parent_information <- parent_information
  dat[, c(
    "artificial_cell_id", "Study", "DF_partition", "method", "parent_1", "parent_2",
    "scDblFinder_origin_cluster", "parent_information", "construction"
  ), drop = FALSE]
}
make_artificial_source_plot <- function(dat, study) {
  dat <- dat[order(dat$profile_origin), , drop = FALSE]
  plot_base(dat, study, "Real cells and method-specific artificial profiles") +
    geom_point(aes(colour = profile_origin), size = 0.42, alpha = 0.80) +
    scale_colour_manual(values = c(
      "Real cell" = "#BDBDBD", "DoubletFinder artificial" = "#E64B35",
      "scDblFinder artificial" = "#4DBBD5", "Scrublet artificial" = "#00A087"
    )) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.2, alpha = 1)))
}
make_real_overlay_plot <- function(dat, field, title, positive_label, positive_colour, study) {
  background <- dat[dat$profile_origin != "Real cell", , drop = FALSE]
  real <- dat[dat$profile_origin == "Real cell", , drop = FALSE]
  real$overlay <- ifelse(as.character(real[[field]]) == positive_label, positive_label, "Other real cells")
  real$overlay <- factor(real$overlay, levels = c("Other real cells", positive_label))
  real <- real[order(real$overlay), , drop = FALSE]
  plot_base(dat, study, title) +
    geom_point(data = background, colour = "#F0F0F0", size = 0.32, alpha = 0.35) +
    geom_point(data = real, aes(colour = overlay), size = 0.52, alpha = 0.90) +
    scale_colour_manual(values = c("Other real cells" = "#BDBDBD", setNames(positive_colour, positive_label))) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
}

if (any(!file.exists(c(current_object_file, file.path(input_dir, "partition_manifest.csv"), file.path(table_dir, "real_cell_metadata_and_existing_calls.csv"))))) {
  stop("Stage-1 inputs are missing; run stage 1 then the Scrublet stage first.")
}
if (!dir.exists(scrublet_dir)) stop("Scrublet artificial output directory is missing; run stage 2 first.")

message("Loading real pre-removal MG00 cells and artificial profile manifests ...")
current <- qs::qread(current_object_file)
counts_all <- get_counts(current)
meta <- read.csv(file.path(table_dir, "real_cell_metadata_and_existing_calls.csv"), stringsAsFactors = FALSE, check.names = FALSE)
rownames(meta) <- meta$Cell
if (!identical(colnames(counts_all), meta$Cell)) stop("Current RNA count columns and saved real-cell metadata disagree.")
manifest <- read.csv(file.path(input_dir, "partition_manifest.csv"), stringsAsFactors = FALSE, check.names = FALSE)
if (!setequal(study_order, unique(as.character(manifest$Study)))) stop("Partition manifest lacks one or more expected studies.")

all_summary <- list()
all_parent_rows <- list()
all_plot_stems <- character()
object_index <- list()
for (study in study_order) {
  study_meta <- meta[meta$Study == study, , drop = FALSE]
  real_cells <- study_meta$Cell
  real_counts <- counts_all[, real_cells, drop = FALSE]
  study_manifest <- manifest[manifest$Study == study, , drop = FALSE]
  study_manifest <- study_manifest[order(study_manifest$DF_partition), , drop = FALSE]

  df_blocks <- list()
  sc_blocks <- list()
  scrublet_blocks <- list()
  df_parent_blocks <- list()
  sc_origin_blocks <- list()
  scrublet_parent_blocks <- list()
  for (row_id in seq_len(nrow(study_manifest))) {
    item <- study_manifest[row_id, , drop = FALSE]
    part_dir <- file.path(input_dir, item$safe_partition)
    scrublet_part_dir <- file.path(scrublet_dir, item$safe_partition)
    df_blocks[[row_id]] <- align_to_feature_union(read_sparse_bundle(
      part_dir, "DoubletFinder_artificial_counts_genes_by_cells.mtx",
      "DoubletFinder_artificial_cells.tsv", "DoubletFinder_artificial_genes.tsv"
    ), rownames(real_counts))
    sc_blocks[[row_id]] <- align_to_feature_union(read_sparse_bundle(
      part_dir, "scDblFinder_artificial_counts_genes_by_cells.mtx",
      "scDblFinder_artificial_cells.tsv", "scDblFinder_artificial_genes.tsv"
    ), rownames(real_counts))
    scrublet_blocks[[row_id]] <- align_to_feature_union(read_sparse_bundle(
      scrublet_part_dir, "Scrublet_artificial_counts_genes_by_cells.mtx",
      "Scrublet_artificial_cells.tsv", "Scrublet_artificial_genes.tsv"
    ), rownames(real_counts))
    df_parent_blocks[[row_id]] <- read.csv(file.path(part_dir, "DoubletFinder_artificial_parent_pairs.csv"), stringsAsFactors = FALSE, check.names = FALSE)
    sc_origin_blocks[[row_id]] <- read.csv(file.path(part_dir, "scDblFinder_artificial_origin_clusters.csv"), stringsAsFactors = FALSE, check.names = FALSE)
    scrublet_parent_blocks[[row_id]] <- read.csv(file.path(scrublet_part_dir, "Scrublet_artificial_parent_pairs.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  }
  df_counts <- do.call(cbind, df_blocks)
  sc_counts <- do.call(cbind, sc_blocks)
  scrublet_counts <- do.call(cbind, scrublet_blocks)
  combined_counts <- cbind(real_counts, df_counts, sc_counts, scrublet_counts)
  if (anyDuplicated(colnames(combined_counts))) stop(study, ": a real/artificial cell ID is duplicated.")

  # Build a concise, explicit provenance table rather than carrying every
  # author metadata column into the artificial-profile Seurat object.
  real_provenance <- data.frame(
    profile_origin = "Real cell", Study = study, DF_partition = study_meta$DF_partition,
    historical_MG0 = study_meta$historical_MG0,
    historical_label_available = study_meta$historical_label_available,
    MG00_status = study_meta$MG00_status,
    DoubletFinder_call = study_meta$DoubletFinder_call,
    scDblFinder_call = study_meta$scDblFinder_call,
    Scrublet_call = study_meta$Scrublet_call,
    stringsAsFactors = FALSE, row.names = real_cells
  )
  artificial_provenance <- function(cell_ids, origin) data.frame(
    # `data.frame()` uses a one-row data frame when every supplied value is a
    # scalar.  Make the leading field explicitly cell-length so its row names
    # are aligned one-to-one with every artificial profile.
    profile_origin = rep(origin, length(cell_ids)), Study = study, DF_partition = NA_character_,
    historical_MG0 = NA, historical_label_available = NA, MG00_status = NA_character_,
    DoubletFinder_call = NA_character_, scDblFinder_call = NA_character_, Scrublet_call = NA_character_,
    stringsAsFactors = FALSE, row.names = cell_ids
  )
  combined_meta <- rbind(
    real_provenance,
    artificial_provenance(colnames(df_counts), "DoubletFinder artificial"),
    artificial_provenance(colnames(sc_counts), "scDblFinder artificial"),
    artificial_provenance(colnames(scrublet_counts), "Scrublet artificial")
  )
  combined_meta$profile_origin <- factor(combined_meta$profile_origin, levels = c(
    "Real cell", "DoubletFinder artificial", "scDblFinder artificial", "Scrublet artificial"
  ))
  if (!identical(rownames(combined_meta), colnames(combined_counts))) stop(study, ": combined provenance is not aligned to counts.")

  # Select variable features only on real cells. This ensures that artificial
  # profiles can be projected and visualised without determining the feature
  # panel. No Harmony or other batch correction is applied to artificial data.
  real_for_hvg <- Seurat::CreateSeuratObject(real_counts, project = safe_name(study))
  real_for_hvg <- Seurat::NormalizeData(real_for_hvg, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  real_for_hvg <- Seurat::FindVariableFeatures(real_for_hvg, selection.method = "vst",
    nfeatures = min(3000L, nrow(real_for_hvg) - 1L), verbose = FALSE)
  features <- Seurat::VariableFeatures(real_for_hvg)
  if (length(features) < 500L) stop(study, ": fewer than 500 real-cell variable features are available.")
  artificial_object <- Seurat::CreateSeuratObject(combined_counts, meta.data = combined_meta, project = safe_name(study))
  artificial_object <- Seurat::NormalizeData(artificial_object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  artificial_object <- Seurat::ScaleData(artificial_object, features = features, verbose = FALSE)
  n_pcs <- min(30L, length(features) - 1L, ncol(artificial_object) - 1L)
  artificial_object <- Seurat::RunPCA(artificial_object, features = features, npcs = n_pcs, verbose = FALSE)
  artificial_object <- Seurat::RunUMAP(artificial_object, dims = seq_len(n_pcs), reduction.name = "umap_artificial_doublet",
    seed.use = 20260908, verbose = FALSE)
  umap <- Seurat::Embeddings(artificial_object, "umap_artificial_doublet")
  plot_data <- artificial_object[[]]
  plot_data$Cell <- rownames(plot_data)
  plot_data$Artificial_UMAP_1 <- umap[plot_data$Cell, 1]
  plot_data$Artificial_UMAP_2 <- umap[plot_data$Cell, 2]

  plots <- list(
    artificial_sources = make_artificial_source_plot(plot_data, study),
    DoubletFinder_calls = make_real_overlay_plot(plot_data, "DoubletFinder_call", "DoubletFinder calls on real cells", "Doublet", "#D7301F", study),
    scDblFinder_calls = make_real_overlay_plot(plot_data, "scDblFinder_call", "scDblFinder calls on real cells", "Doublet", "#D7301F", study),
    Scrublet_calls = make_real_overlay_plot(plot_data, "Scrublet_call", "Scrublet calls on real cells", "Doublet", "#D7301F", study),
    original_MG0 = make_real_overlay_plot(plot_data, "historical_MG0", "Historical original MG0 on real cells", "TRUE", "#2C7FB8", study),
    MG00 = make_real_overlay_plot(plot_data, "MG00_status", "MG00 on real cells", "MG00 (manual cluster 2 transfer)", "#C51B7D", study)
  )
  for (plot_name in names(plots)) {
    stem <- file.path(figure_dir, paste0("UMAP_", safe_name(study), "_", plot_name))
    save_plot_both(plots[[plot_name]], stem)
    print(plots[[plot_name]])
    all_plot_stems <- c(all_plot_stems, stem)
  }

  # Preserve a full, study-specific R object: raw artificial and real matrices,
  # origin metadata, and the visualisation PCA/UMAP reductions are all retained.
  object_file <- file.path(object_dir, paste0(safe_name(study), "_real_and_artificial_doublet_UMAP.qs"))
  qs::qsave(artificial_object, object_file, preset = "high")
  write.csv(plot_data, file.path(table_dir, paste0(safe_name(study), "_real_and_artificial_UMAP_metadata.csv")), row.names = FALSE)
  object_index[[study]] <- object_file

  all_summary[[study]] <- data.frame(
    Study = study,
    profile_origin = c("Real cell", "DoubletFinder artificial", "scDblFinder artificial", "Scrublet artificial"),
    n_profiles = c(ncol(real_counts), ncol(df_counts), ncol(sc_counts), ncol(scrublet_counts)),
    stringsAsFactors = FALSE
  )
  all_parent_rows[[study]] <- rbind(
    standardise_parent_provenance(do.call(rbind, df_parent_blocks), "DoubletFinder", "individual parent IDs"),
    standardise_parent_provenance(do.call(rbind, sc_origin_blocks), "scDblFinder", "official helper origin cluster; individual IDs unavailable"),
    standardise_parent_provenance(do.call(rbind, scrublet_parent_blocks), "Scrublet", "individual parent IDs")
  )
  message("Built and saved artificial-doublet UMAP object for ", study, " (", ncol(artificial_object), " total profiles).")
  rm(real_counts, df_counts, sc_counts, scrublet_counts, combined_counts, artificial_object, real_for_hvg, plot_data)
  invisible(gc())
}

profile_summary <- do.call(rbind, all_summary)
profile_summary$profile_origin <- factor(profile_summary$profile_origin,
  levels = c("Real cell", "DoubletFinder artificial", "scDblFinder artificial", "Scrublet artificial"))
write.csv(profile_summary, file.path(table_dir, "artificial_profile_counts_by_study_method.csv"), row.names = FALSE)
all_parent_data <- do.call(rbind, all_parent_rows)
write.csv(all_parent_data, file.path(table_dir, "artificial_profile_parent_provenance.csv"), row.names = FALSE)
qs::qsave(object_index, file.path(object_dir, "per_study_artificial_doublet_UMAP_object_index.qs"), preset = "high")

p_counts <- ggplot(profile_summary, aes(Study, n_profiles, fill = profile_origin)) +
  geom_col(position = position_dodge(width = 0.80), width = 0.72) +
  scale_fill_manual(values = c(
    "Real cell" = "#BDBDBD", "DoubletFinder artificial" = "#E64B35",
    "scDblFinder artificial" = "#4DBBD5", "Scrublet artificial" = "#00A087"
  )) +
  labs(title = "Profiles per artificial-doublet UMAP", x = NULL, y = "Number of profiles", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
save_plot_both(p_counts, file.path(figure_dir, "Artificial_profile_counts_by_study_method"), 10.5, 6)
print(p_counts)
all_plot_stems <- c(all_plot_stems, file.path(figure_dir, "Artificial_profile_counts_by_study_method"))

real_summary <- meta
real_summary$historical_MG0_label <- ifelse(real_summary$historical_MG0, "Historical MG0", "Other real cell")
real_summary$MG00_label <- ifelse(real_summary$MG00_status == "MG00 (manual cluster 2 transfer)", "MG00", "Other real cell")
historical_label_counts <- aggregate(Cell ~ Study + historical_MG0_label, data = real_summary, FUN = length)
names(historical_label_counts) <- c("Study", "label", "n_real_cells")
historical_label_counts$definition <- "Historical label"
MG00_label_counts <- aggregate(Cell ~ Study + MG00_label, data = real_summary, FUN = length)
names(MG00_label_counts) <- c("Study", "label", "n_real_cells")
MG00_label_counts$definition <- "MG00 reanalysis label"
label_counts <- rbind(historical_label_counts, MG00_label_counts)
p_labels <- ggplot(label_counts, aes(Study, n_real_cells, fill = label)) +
  geom_col(position = "stack", width = 0.72) +
  facet_wrap(~definition, nrow = 1) +
  scale_fill_manual(values = c("Other real cell" = "#BDBDBD", "Historical MG0" = "#2C7FB8", "MG00" = "#C51B7D")) +
  labs(title = "Real-cell labels projected onto artificial-doublet UMAPs", x = NULL, y = "Number of real cells", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
save_plot_both(p_labels, file.path(figure_dir, "Real_cell_historical_MG0_and_MG00_counts"), 10.5, 6)
print(p_labels)
all_plot_stems <- c(all_plot_stems, file.path(figure_dir, "Real_cell_historical_MG0_and_MG00_counts"))

pair_check <- data.frame(
  figure_stem = all_plot_stems,
  PNG_exists = file.exists(paste0(all_plot_stems, ".png")),
  PDF_exists = file.exists(paste0(all_plot_stems, ".pdf")),
  stringsAsFactors = FALSE
)
write.csv(pair_check, file.path(table_dir, "PNG_PDF_pair_check.csv"), row.names = FALSE)
if (any(!pair_check$PNG_exists | !pair_check$PDF_exists)) stop("At least one artificial-doublet figure lacks a PNG or PDF counterpart.")

writeLines(c(
  "Objects: one study-specific Seurat object contains real cells plus DoubletFinder, scDblFinder, and Scrublet artificial profiles, their raw matrix columns, provenance metadata, PCA and UMAP.",
  "Real cells: current MG00 pre-removal diagnostic universe. Existing real-cell calls are displayed but never recalculated or changed here.",
  "Artificial construction: DoubletFinder profile mean with pN=0.25; scDblFinder official getArtificialDoublets helper with propRandom=0; Scrublet official simulate_doublets with raw-count addition and ratio 1/3.",
  "Features: selected only from the real cells of each study; artificial profiles did not select variable genes.",
  "Embedding: LogNormalize, ScaleData, PCA and UMAP separately per study. No Harmony/batch correction and no artificial-versus-real filtering were applied.",
  "Interpretation: artificial profiles visualise the simulated neighborhood geometry. They are diagnostic references rather than observed cells or proof of a real-cell label."
), file.path(runtime_dir, "stage3_method_manifest.txt"))
writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "stage3_R_sessionInfo.txt"))
message("Artificial-doublet visualisation is complete: ", figure_dir)
