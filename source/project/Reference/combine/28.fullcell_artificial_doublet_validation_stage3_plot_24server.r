#!/usr/bin/env Rscript

# Corrected artificial-doublet validation, stage 3.
#
# Scope: all real profiles are the complete QC-passed real-cell universe in
# each study. A fixed PCA/UMAP reference is learned from those real cells
# only. Method-specific artificial profiles are projected into that fixed
# reference; they never change the real-cell manifold, QC status, labels,
# clusters, or doublet calls.

rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
invisible(gc())
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(qs)
  library(ggplot2)
  library(uwot)
})

set.seed(20260913)
project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
combine_root <- file.path(project_root, "combine")
diagnosis_root <- file.path(combine_root, "results/fullcell_historical_MG0_MG5_doublet_diagnosis_v2")
result_root <- file.path(combine_root, "results/fullcell_artificial_doublet_validation_v2")
input_dir <- file.path(result_root, "input")
scrublet_dir <- file.path(result_root, "scrublet_output")
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
object_dir <- file.path(result_root, "object")
runtime_dir <- file.path(result_root, "runtime")
for (directory in c(figure_dir, table_dir, object_dir, runtime_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)

study_order <- c("Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
                 "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024")
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width = 7.5, height = 5) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = grDevices::pdf, bg = "white", limitsize = FALSE)
}
get_layer <- function(object, layer) {
  SeuratObject::LayerData(object, assay = Seurat::DefaultAssay(object), layer = layer)
}
read_bundle <- function(folder, prefix) {
  files <- c(file.path(folder, paste0(prefix, "_counts_genes_by_cells.mtx")),
             file.path(folder, paste0(prefix, "_cells.tsv")),
             file.path(folder, paste0(prefix, "_genes.tsv")))
  if (any(!file.exists(files))) stop("Missing artificial-profile bundle: ", paste(files, collapse = "; "))
  counts <- as(Matrix::readMM(files[[1]]), "dgCMatrix")
  cells <- read.delim(files[[2]], stringsAsFactors = FALSE, check.names = FALSE)$cell_id
  genes <- read.delim(files[[3]], stringsAsFactors = FALSE, check.names = FALSE)$gene_id
  if (!identical(dim(counts), c(length(genes), length(cells)))) stop("Sparse artificial matrix dimensions disagree with its manifest.")
  rownames(counts) <- genes; colnames(counts) <- cells
  counts
}
restrict_to_features <- function(counts, features) {
  hit <- match(features, rownames(counts))
  if (all(!is.na(hit))) return(counts[hit, , drop = FALSE])
  present <- which(!is.na(hit))
  if (!length(present)) stop("Artificial profile bundle has no reference PCA features.")
  out <- Matrix::sparseMatrix(i = integer(), j = integer(), dims = c(length(features), ncol(counts)),
                              dimnames = list(features, colnames(counts)))
  out[present, ] <- counts[hit[present], , drop = FALSE]
  as(out, "dgCMatrix")
}
normalise_log_sparse <- function(counts, scale_factor = 10000) {
  library_sizes <- pmax(Matrix::colSums(counts), 1)
  normalised <- counts %*% Matrix::Diagonal(x = scale_factor / library_sizes)
  normalised@x <- log1p(normalised@x)
  as(normalised, "dgCMatrix")
}
reference_scaling <- function(real_data, features) {
  x <- real_data[features, , drop = FALSE]
  n <- ncol(x)
  center <- Matrix::rowMeans(x)
  variance <- (Matrix::rowSums(x ^ 2) - n * center ^ 2) / max(n - 1, 1)
  scale <- sqrt(pmax(variance, 1e-8))
  list(center = as.numeric(center), scale = as.numeric(scale))
}
project_counts_to_pca <- function(counts, features, loadings, center, scale) {
  selected <- restrict_to_features(counts, features)
  normalised <- normalise_log_sparse(selected)
  scaled_loadings <- sweep(loadings, 1, scale, "/")
  coordinates <- as.matrix(Matrix::t(normalised) %*% scaled_loadings)
  offset <- colSums(loadings * (center / scale))
  sweep(coordinates, 2, offset, "-")
}
plot_base <- function(data, title, xlim, ylim) {
  ggplot(data, aes(Artificial_UMAP_1, Artificial_UMAP_2)) +
    coord_equal(xlim = xlim, ylim = ylim) +
    labs(title = title,
         subtitle = "Fixed UMAP from all QC-passed real cells; artificial profiles are projected by the fixed real-cell PCA/UMAP model",
         x = "Full-cell reference UMAP 1", y = "Full-cell reference UMAP 2") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = .5),
          plot.subtitle = element_text(hjust = .5, size = 8.5),
          legend.title = element_blank(), legend.position = "bottom")
}
make_artificial_plot <- function(data, method, title, xlim, ylim) {
  artificial_name <- paste0(method, " artificial")
  data <- data[data$profile_origin %in% c("Real cell", artificial_name), , drop = FALSE]
  data$overlay <- ifelse(data$profile_origin == artificial_name, artificial_name, "All QC-passed real cells")
  data$overlay <- factor(data$overlay, levels = c("All QC-passed real cells", artificial_name))
  data <- data[order(data$overlay), , drop = FALSE]
  colour <- c("DoubletFinder" = "#E64B35", "scDblFinder" = "#4DBBD5", "Scrublet" = "#00A087")[[method]]
  plot_base(data, title, xlim, ylim) +
    geom_point(aes(colour = overlay), size = .32, alpha = .72) +
    scale_colour_manual(values = c("All QC-passed real cells" = "#D9D9D9", setNames(colour, artificial_name))) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
}
make_real_call_plot <- function(data, method, title, xlim, ylim) {
  real <- data[data$profile_origin == "Real cell", , drop = FALSE]
  field <- paste0(method, "_call")
  real$overlay <- factor(as.character(real[[field]]), levels = c("Singlet", "Doublet"))
  real <- real[order(real$overlay), , drop = FALSE]
  plot_base(real, title, xlim, ylim) +
    geom_point(aes(colour = overlay), size = .36, alpha = .80) +
    scale_colour_manual(values = c(Singlet = "#BDBDBD", Doublet = "#7A1FA2")) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
}
make_real_label_plot <- function(data, field, positive, colour, title, xlim, ylim) {
  real <- data[data$profile_origin == "Real cell", , drop = FALSE]
  real$overlay <- ifelse(as.logical(real[[field]]), positive, "Other QC-passed real cells")
  real$overlay <- factor(real$overlay, levels = c("Other QC-passed real cells", positive))
  real <- real[order(real$overlay), , drop = FALSE]
  plot_base(real, title, xlim, ylim) +
    geom_point(aes(colour = overlay), size = .36, alpha = .80) +
    scale_colour_manual(values = c("Other QC-passed real cells" = "#D9D9D9", setNames(colour, positive))) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
}

manifest <- read.csv(file.path(input_dir, "partition_manifest.csv"), stringsAsFactors = FALSE, check.names = FALSE)
all_summaries <- list(); all_metadata <- list(); all_provenance <- list(); all_stems <- character(); object_index <- list()

for (study in study_order) {
  safe <- safe_name(study)
  source_file <- file.path(diagnosis_root, "object", paste0(safe, "_all_QC_passed_fullcell_all_three_calls.qs"))
  if (!file.exists(source_file)) stop("Missing corrected complete real-cell object for ", study)
  object <- qs::qread(source_file)
  real_meta <- object[[]]
  real_data <- get_layer(object, "data")
  real_pca <- Seurat::Embeddings(object, "pca")
  all_loadings <- Seurat::Loadings(object, "pca")
  npcs <- min(30L, ncol(real_pca), ncol(all_loadings))
  features <- rownames(all_loadings)[seq_len(min(3000L, nrow(all_loadings)))]
  loadings <- all_loadings[features, seq_len(npcs), drop = FALSE]
  scaling <- reference_scaling(real_data, features)
  rows <- manifest[manifest$Study == study, , drop = FALSE]
  if (!nrow(rows)) stop("No technical partitions available for ", study)

  # Learn UMAP once from the complete real-cell PCA space, retaining the model
  # only to transform the artificial profiles afterwards.
  set.seed(20260913)
  umap_fit <- uwot::umap(real_pca[, seq_len(npcs), drop = FALSE], n_neighbors = 30,
                         n_components = 2, metric = "euclidean", n_epochs = 200,
                         ret_model = TRUE, n_threads = 1, verbose = FALSE)
  real_umap <- umap_fit$embedding
  colnames(real_umap) <- c("Artificial_UMAP_1", "Artificial_UMAP_2")
  real_plot <- data.frame(Cell = rownames(real_meta), profile_origin = "Real cell", Study = study,
                          DF_partition = as.character(real_meta$DF_partition),
                          historical_MG0 = as.logical(real_meta$historical_MG0), MG00 = as.logical(real_meta$MG00),
                          DoubletFinder_call = as.character(real_meta$DoubletFinder_call),
                          scDblFinder_call = as.character(real_meta$scDblFinder_call),
                          Scrublet_call = as.character(real_meta$Scrublet_call),
                          Artificial_UMAP_1 = real_umap[, 1], Artificial_UMAP_2 = real_umap[, 2],
                          stringsAsFactors = FALSE, row.names = rownames(real_meta))

  # Project all technical partitions for one method in one transform call.
  # This retains every artificial profile while avoiding 91 separate UMAP
  # transform initialisations.
  artificial_plot_blocks <- list(); provenance_blocks <- list(); profile_count <- c("Real cell" = nrow(real_plot))
  for (method in c("DoubletFinder", "scDblFinder", "Scrublet")) {
    pca_blocks <- vector("list", nrow(rows)); cell_blocks <- vector("list", nrow(rows)); partition_blocks <- vector("list", nrow(rows))
    for (i in seq_len(nrow(rows))) {
      row <- rows[i, , drop = FALSE]
      input_folder <- file.path(input_dir, row$safe_partition)
      source_folder <- if (method == "Scrublet") file.path(scrublet_dir, row$safe_partition) else input_folder
      counts <- read_bundle(source_folder, paste0(method, "_artificial"))
      pca_blocks[[i]] <- project_counts_to_pca(counts, features, loadings, scaling$center, scaling$scale)
      cell_blocks[[i]] <- colnames(counts)
      partition_blocks[[i]] <- rep(as.character(row$DF_partition), ncol(counts))
      if (method %in% c("DoubletFinder", "Scrublet")) {
        pairs <- read.csv(file.path(source_folder, paste0(method, "_artificial_parent_pairs.csv")), stringsAsFactors = FALSE, check.names = FALSE)
        pairs$method <- method
        provenance_blocks[[length(provenance_blocks) + 1L]] <- pairs
      }
      rm(counts); invisible(gc(FALSE))
    }
    pca_projection <- do.call(rbind, pca_blocks)
    umap_projection <- uwot::umap_transform(pca_projection, umap_fit, n_threads = 1)
    cells <- unlist(cell_blocks, use.names = FALSE); partitions <- unlist(partition_blocks, use.names = FALSE)
    artificial_plot_blocks[[method]] <- data.frame(
      Cell = cells, profile_origin = paste0(method, " artificial"), Study = study,
      DF_partition = partitions, historical_MG0 = NA, MG00 = NA,
      DoubletFinder_call = NA, scDblFinder_call = NA, Scrublet_call = NA,
      Artificial_UMAP_1 = umap_projection[, 1], Artificial_UMAP_2 = umap_projection[, 2],
      stringsAsFactors = FALSE, row.names = cells)
    profile_count[paste0(method, " artificial")] <- length(cells)
    rm(pca_blocks, cell_blocks, partition_blocks, pca_projection, umap_projection, cells, partitions); invisible(gc())
  }
  artificial_plot <- do.call(rbind, artificial_plot_blocks)
  plot_data <- rbind(real_plot, artificial_plot)
  xlim <- range(plot_data$Artificial_UMAP_1, finite = TRUE)
  ylim <- range(plot_data$Artificial_UMAP_2, finite = TRUE)
  plots <- list(
    DoubletFinder_artificial = make_artificial_plot(plot_data, "DoubletFinder", paste0(study, ": DoubletFinder artificial profiles"), xlim, ylim),
    DoubletFinder_real_calls = make_real_call_plot(plot_data, "DoubletFinder", paste0(study, ": DoubletFinder-labelled real doublets"), xlim, ylim),
    scDblFinder_artificial = make_artificial_plot(plot_data, "scDblFinder", paste0(study, ": scDblFinder artificial profiles"), xlim, ylim),
    scDblFinder_real_calls = make_real_call_plot(plot_data, "scDblFinder", paste0(study, ": scDblFinder-labelled real doublets"), xlim, ylim),
    Scrublet_artificial = make_artificial_plot(plot_data, "Scrublet", paste0(study, ": Scrublet artificial profiles"), xlim, ylim),
    Scrublet_real_calls = make_real_call_plot(plot_data, "Scrublet", paste0(study, ": Scrublet-labelled real doublets"), xlim, ylim),
    historical_MG0 = make_real_label_plot(plot_data, "historical_MG0", "MG0", "#2C7FB8", paste0(study, ": MG0 (original MG0--MG5 label)"), xlim, ylim),
    MG00 = make_real_label_plot(plot_data, "MG00", "MG00", "#C51B7D", paste0(study, ": MG00"), xlim, ylim)
  )
  for (name in names(plots)) {
    stem <- file.path(figure_dir, paste0("UMAP_", safe, "_", name))
    save_plot_both(plots[[name]], stem); print(plots[[name]]); all_stems <- c(all_stems, stem)
  }
  audit_object <- list(
    study = study, real_cell_source_qs = source_file,
    fixed_real_pca_reduction = "pca", fixed_real_umap_method = "uwot::umap on all QC-passed real-cell PCA coordinates",
    projection_method = "library-size normalisation + fixed PCA loadings + uwot::umap_transform",
    reference_features = features, profile_metadata = plot_data,
    parent_pair_files = list(input_dir = input_dir, scrublet_dir = scrublet_dir)
  )
  audit_path <- file.path(object_dir, paste0(safe, "_all_QC_passed_real_and_method_artificial_projected_reference.qs"))
  qs::qsave(audit_object, audit_path, preset = "high")
  write.csv(plot_data, file.path(table_dir, paste0(safe, "_all_QC_passed_real_and_method_artificial_projected_reference_metadata.csv")), row.names = FALSE)
  object_index[[study]] <- audit_path
  all_summaries[[study]] <- data.frame(Study = study, profile_origin = names(profile_count), n_profiles = unname(profile_count), stringsAsFactors = FALSE)
  all_metadata[[study]] <- plot_data[, c("Cell", "profile_origin", "Study", "DF_partition", "historical_MG0", "MG00", "DoubletFinder_call", "scDblFinder_call", "Scrublet_call", "Artificial_UMAP_1", "Artificial_UMAP_2"), drop = FALSE]
  if (length(provenance_blocks)) all_provenance[[study]] <- do.call(rbind, provenance_blocks)
  message("Saved fixed full-cell reference and projected artificial profiles for ", study, " with ", nrow(plot_data), " profiles.")
  rm(object, real_meta, real_data, real_pca, all_loadings, loadings, scaling, umap_fit, real_umap, real_plot,
     artificial_plot_blocks, artificial_plot, provenance_blocks, plot_data, plots, audit_object); invisible(gc())
}

summary <- do.call(rbind, all_summaries)
metadata <- do.call(rbind, all_metadata)
pair_check <- data.frame(figure_stem = all_stems, PNG_exists = file.exists(paste0(all_stems, ".png")), PDF_exists = file.exists(paste0(all_stems, ".pdf")))
if (any(!pair_check$PNG_exists | !pair_check$PDF_exists)) stop("A full-cell artificial-validation figure lacks PNG or PDF.")
write.csv(summary, file.path(table_dir, "profile_counts_by_study_method.csv"), row.names = FALSE)
write.csv(metadata, file.path(table_dir, "all_studies_fullcell_artificial_validation_projected_reference_metadata.csv"), row.names = FALSE)
if (length(all_provenance)) write.csv(do.call(rbind, all_provenance), file.path(table_dir, "artificial_parent_pairs_auditable_subset.csv"), row.names = FALSE)
write.csv(pair_check, file.path(table_dir, "PNG_PDF_pair_check.csv"), row.names = FALSE)
qs::qsave(object_index, file.path(object_dir, "per_study_fullcell_artificial_validation_object_index.qs"), preset = "high")
writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "stage3_R_sessionInfo.txt"))
writeLines(c(
  "Scope correction: every real profile is a QC-passed cell from the complete study object. MG0 and MG00 are display-only overlays.",
  "Each study learns a UMAP from all QC-passed real-cell PCA coordinates only; every artificial profile is subsequently transformed into that fixed reference.",
  "Artificial and caller-labelled real doublets are intentionally shown in separate complete figures; matched pairs share coordinates and axes.",
  "Artificial profiles are visual references only and never change real-cell QC, labels, clustering or calls."
), file.path(runtime_dir, "stage3_scope_manifest.txt"))
message("Corrected full-cell artificial-validation stage 3 complete: ", figure_dir)
