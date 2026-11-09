#!/usr/bin/env Rscript

# OL–MG artificial-doublet reference, stage 3
# Compare real MG00 and MG0 with controlled artificial OL–MG,
# OL–OL and MG–MG profiles.  Artificial profiles are descriptive references:
# no artificial profile participates in real-cell QC, doublet calling, labels,
# or biological differential-expression conclusions.

local({
  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
    library(qs)
    library(ggplot2)
    library(patchwork)
    library(RANN)
  })
  set.seed(20260909)
  project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
  combine_root <- file.path(project_root, "combine")
  stage2_root <- file.path(combine_root, "results/OL_MG_artificial_reference/02_artificial_profile_construction")
  input_object_dir <- file.path(combine_root, "results/full_count_reanalysis/01_per_study/object")
  result_root <- file.path(combine_root, "results/OL_MG_artificial_reference/03_comparison")
  figure_dir <- file.path(result_root, "figure")
  table_dir <- file.path(result_root, "table")
  runtime_dir <- file.path(result_root, "runtime")
  de_dir <- file.path(table_dir, "real_MG00_vs_MG0_DE")
  for (d in c(figure_dir, table_dir, runtime_dir, de_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  study_order <- c("Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021", "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024")
  safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  ol_markers <- c("PLP1", "MOBP", "MAG", "CNP", "CLDN11", "SOX10", "MYRF", "OPALIN")
  mg_markers <- c("P2RY12", "TMEM119", "CX3CR1", "GPR34", "CSF1R", "SPI1", "C1QA", "C1QB", "C1QC", "AIF1", "TYROBP", "FCER1G")
  if (!requireNamespace("RANN", quietly = TRUE)) stop("RANN is required for nearest artificial-profile distance calculations.")

  save_plot_both <- function(plot, stem, width, height) {
    ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 320, bg = "white", limitsize = FALSE)
    ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = grDevices::pdf, bg = "white", limitsize = FALSE)
    if (!file.exists(paste0(stem, ".png")) || !file.exists(paste0(stem, ".pdf"))) stop("Incomplete PNG/PDF pair: ", stem)
  }
  get_counts <- function(object) {
    out <- tryCatch(SeuratObject::LayerData(object, assay = "RNA", layer = "counts"), error = function(e) NULL)
    if (is.null(out)) out <- tryCatch(SeuratObject::GetAssayData(object, assay = "RNA", layer = "counts"), error = function(e) NULL)
    if (is.null(out)) out <- SeuratObject::GetAssayData(object, assay = "RNA", slot = "counts")
    if (!inherits(out, "sparseMatrix")) out <- Matrix::Matrix(out, sparse = TRUE)
    out
  }
  assert_integer <- function(x, label) if (length(x@x) && any(abs(x@x - round(x@x)) > 1e-8)) stop(label, " has non-integer count values.")
  panel_score <- function(counts, markers) {
    present <- intersect(markers, rownames(counts))
    if (length(present) < 3L) stop("Too few panel markers available: ", paste(markers, collapse = ","))
    subset <- counts[present, , drop = FALSE]
    lib <- Matrix::colSums(counts)
    normalized <- log1p(subset %*% Matrix::Diagonal(x = 10000 / pmax(lib, 1)))
    list(score = as.numeric(Matrix::colMeans(normalized)), present = present)
  }
  marker_dot_data <- function(counts, groups, markers) {
    markers <- intersect(markers, rownames(counts))
    out <- list()
    for (g in levels(groups)) {
      cells <- which(groups == g)
      if (!length(cells)) next
      sub <- counts[markers, cells, drop = FALSE]
      norm <- log1p(sub %*% Matrix::Diagonal(x = 10000 / pmax(Matrix::colSums(counts[, cells, drop = FALSE]), 1)))
      out[[g]] <- data.frame(group = g, gene = markers, avg_expression = as.numeric(Matrix::rowMeans(norm)), pct_expressed = as.numeric(Matrix::rowMeans(sub > 0)), stringsAsFactors = FALSE)
    }
    do.call(rbind, out)
  }
  plot_dot <- function(dot_data, title) {
    ggplot(dot_data, aes(gene, group)) + geom_point(aes(size = pct_expressed, colour = avg_expression)) +
      scale_colour_gradient(low = "grey90", high = "#54278F") + scale_size(range = c(0.5, 8)) +
      labs(title = title, x = NULL, y = NULL, colour = "Average\nexpression", size = "Percent\nexpressed") +
      theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, face = "bold"))
  }
  make_group <- function(meta) {
    group <- rep("Other real QC-passed cell", nrow(meta))
    art <- meta$cell_origin == "artificial_profile"
    group[art] <- paste(meta$profile_class[art], meta$construction[art], sep = " | ")
    group[!art & meta$MG0 %in% TRUE] <- "MG0"
    group[!art & meta$MG00 %in% TRUE] <- "MG00"
    factor(group, levels = c("Other real QC-passed cell", "MG0", "MG00", "OL_MG | raw_count_sum", "OL_MG | parent_depth_matched", "OL_OL | raw_count_sum", "OL_OL | parent_depth_matched", "MG_MG | raw_count_sum", "MG_MG | parent_depth_matched"))
  }
  nearest_olmg <- function(meta, pcs) {
    distance <- rep(NA_real_, nrow(meta))
    primary_art <- which(meta$cell_origin == "artificial_profile" & meta$profile_class == "OL_MG" & meta$construction == "parent_depth_matched")
    real <- which(meta$cell_origin == "real_QC_passed")
    for (batch in unique(as.character(meta$batch_id[real]))) {
      rb <- real[as.character(meta$batch_id[real]) == batch]
      ab <- primary_art[as.character(meta$batch_id[primary_art]) == batch]
      if (!length(rb) || !length(ab)) next
      distance[rb] <- RANN::nn2(data = as.matrix(pcs[ab, , drop = FALSE]), query = as.matrix(pcs[rb, , drop = FALSE]), k = 1L)$nn.dists[, 1]
    }
    distance
  }
  match_background <- function(dat, target_group) {
    targets <- which(dat$comparison_group == target_group & dat$cell_origin == "real_QC_passed")
    back <- which(
      dat$MG_parent %in% TRUE &
        dat$cell_origin == "real_QC_passed" &
        !(dat$MG00 %in% TRUE) &
        !(dat$MG0 %in% TRUE)
    )
    out <- list()
    target_dat <- dat[targets, , drop = FALSE]
    for (study in unique(target_dat$Study)) for (batch in unique(target_dat$batch_id[target_dat$Study == study])) {
      t_idx <- targets[dat$Study[targets] == study & dat$batch_id[targets] == batch]
      b_idx <- back[dat$Study[back] == study & dat$batch_id[back] == batch]
      if (!length(t_idx) || !length(b_idx)) next
      q <- as.matrix(dat[t_idx, c("nCount_RNA", "nFeature_RNA")])
      b <- as.matrix(dat[b_idx, c("nCount_RNA", "nFeature_RNA")])
      scale_center <- apply(rbind(q, b), 2, median, na.rm = TRUE)
      scale_spread <- apply(rbind(q, b), 2, mad, na.rm = TRUE); scale_spread[!is.finite(scale_spread) | scale_spread == 0] <- 1
      q <- sweep(sweep(q, 2, scale_center, "-"), 2, scale_spread, "/")
      b <- sweep(sweep(b, 2, scale_center, "-"), 2, scale_spread, "/")
      nearest <- RANN::nn2(data = b, query = q, k = 1L)$nn.idx[, 1]
      out[[length(out) + 1L]] <- data.frame(Study = study, batch_id = batch, target_group = target_group,
        target_Cell = dat$Cell[t_idx], matched_MG_Cell = dat$Cell[b_idx[nearest]], stringsAsFactors = FALSE)
    }
    if (length(out)) do.call(rbind, out) else data.frame(Study = character(), batch_id = character(), target_group = character(), target_Cell = character(), matched_MG_Cell = character())
  }

  all_metrics <- list(); all_external <- list(); all_effects <- list(); all_de <- list(); all_figures <- list(); all_mapping <- list()
  for (study in study_order) {
    safe <- safe_name(study)
    message("\n===== OL–MG comparison: ", study, " =====")
    object_file <- file.path(stage2_root, "object", paste0(safe, "_real_cells_plus_artificial_profiles.qs"))
    qc_file <- file.path(input_object_dir, paste0(safe, "_QC_before_doublet_removal.qs"))
    if (!file.exists(object_file) || !file.exists(qc_file)) stop(study, ": missing stage-2 object or QC source.")
    object <- qs::qread(object_file); qc <- qs::qread(qc_file)
    counts <- get_counts(object); assert_integer(counts, study)
    meta <- object[[]]; meta$Cell <- rownames(meta)
    if (!all(c("cell_origin", "profile_class", "construction", "MG0", "MG00", "MG_parent", "batch_id") %in% names(meta))) stop(study, ": required stage-2 metadata missing.")
    meta$Study <- study
    meta$nCount_RNA <- Matrix::colSums(counts); meta$nFeature_RNA <- Matrix::colSums(counts > 0)
    qc_meta <- qc[[]]
    for (metric in c("percent.mt", "percent.ribo")) meta[[metric]] <- if (metric %in% names(qc_meta)) qc_meta[meta$Cell, metric] else NA_real_
    ol <- panel_score(counts, ol_markers); mg <- panel_score(counts, mg_markers)
    meta$OL_independent_score <- ol$score; meta$MG_independent_score <- mg$score; meta$dual_lineage_score <- pmin(meta$OL_independent_score, meta$MG_independent_score)
    meta$comparison_group <- make_group(meta)
    pcs <- Seurat::Embeddings(object, "pca_real_reference_projection"); pcs <- pcs[meta$Cell, seq_len(min(20L, ncol(pcs))), drop = FALSE]
    meta$nearest_artificial_OL_MG_distance <- nearest_olmg(meta, pcs)
    umap <- as.data.frame(Seurat::Embeddings(object, "umap_real_reference_projection")[meta$Cell, , drop = FALSE]); names(umap) <- c("UMAP_1", "UMAP_2")
    dat <- cbind(meta, umap)
    chosen <- dat$comparison_group != "Other real QC-passed cell"
    group_colours <- c("Other real QC-passed cell" = "#D9D9D9", "MG0" = "#2C7FB8", "MG00" = "#C51B7D", "OL_MG | raw_count_sum" = "#D7301F", "OL_MG | parent_depth_matched" = "#FC8D59", "OL_OL | raw_count_sum" = "#31A354", "OL_OL | parent_depth_matched" = "#74C476", "MG_MG | raw_count_sum" = "#756BB1", "MG_MG | parent_depth_matched" = "#9E9AC8")
    p_umap <- ggplot(dat[order(chosen), ], aes(UMAP_1, UMAP_2, colour = comparison_group)) + geom_point(size = 0.28, alpha = 0.8) + coord_equal() + scale_colour_manual(values = group_colours, drop = FALSE) + labs(title = paste0(study, ": real cells and artificial OL–MG references"), colour = NULL) + theme_classic(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom", legend.text = element_text(size = 7))
    p_score <- (ggplot(dat, aes(UMAP_1, UMAP_2, colour = OL_independent_score)) + geom_point(size = 0.25) + scale_colour_viridis_c() + coord_equal() + labs(title = "Independent OL score") + theme_classic(base_size = 10)) + (ggplot(dat, aes(UMAP_1, UMAP_2, colour = MG_independent_score)) + geom_point(size = 0.25) + scale_colour_viridis_c() + coord_equal() + labs(title = "Independent MG score") + theme_classic(base_size = 10)) + plot_annotation(title = paste0(study, ": independent lineage-score UMAPs"))
    dot_groups <- factor(as.character(dat$comparison_group), levels = setdiff(levels(dat$comparison_group), "Other real QC-passed cell"))
    dot_data <- rbind(marker_dot_data(counts, dot_groups, ol_markers), marker_dot_data(counts, dot_groups, mg_markers))
    p_dot <- plot_dot(dot_data, paste0(study, ": independent OL/MG marker DotPlot"))
    p_distance <- ggplot(dat[dat$comparison_group != "Other real QC-passed cell", ], aes(comparison_group, nearest_artificial_OL_MG_distance, fill = comparison_group)) + geom_boxplot(outlier.size = 0.3, na.rm = TRUE) + scale_fill_manual(values = group_colours) + labs(title = paste0(study, ": distance to depth-matched artificial OL–MG"), x = NULL, y = "Nearest PCA-space distance") + theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
    qcd <- dat[dat$cell_origin == "real_QC_passed" & dat$comparison_group %in% c("MG0", "MG00"), ]
    qcd_long <- reshape(qcd[, c("comparison_group", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")], varying = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), v.names = "value", timevar = "metric", times = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), direction = "long")
    p_qc <- ggplot(qcd_long, aes(comparison_group, value, fill = comparison_group)) + geom_violin(scale = "width", na.rm = TRUE) + geom_boxplot(width = 0.12, outlier.size = 0.2, na.rm = TRUE) + facet_wrap(~metric, scales = "free_y", ncol = 4) + scale_fill_manual(values = group_colours) + labs(title = paste0(study, ": real-cell QC metrics"), x = NULL, y = NULL) + theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
    plots <- list(reference_UMAP = list(p_umap, 10, 7), independent_score_UMAP = list(p_score, 12, 6), independent_marker_DotPlot = list(p_dot, 13, 6), artificial_distance = list(p_distance, 12, 6), real_cell_QC = list(p_qc, 12, 5))
    for (nm in names(plots)) { stem <- file.path(figure_dir, paste0(safe, "_", nm)); save_plot_both(plots[[nm]][[1]], stem, plots[[nm]][[2]], plots[[nm]][[3]]); print(plots[[nm]][[1]]); all_figures[[length(all_figures) + 1L]] <- data.frame(Study = study, figure = nm, stem = stem, stringsAsFactors = FALSE) }
    real <- dat[dat$cell_origin == "real_QC_passed" & dat$comparison_group %in% c("MG0", "MG00"), ]
    for (method in c("DoubletFinder_call", "scDblFinder_call", "Scrublet_call")) if (method %in% names(real)) {
      groups_for_method <- split(real, interaction(real$Study, real$comparison_group, drop = TRUE))
      temp <- do.call(rbind, lapply(groups_for_method, function(x) data.frame(
        Study = x$Study[[1]], comparison_group = as.character(x$comparison_group[[1]]),
        n_available = sum(!is.na(x[[method]])), n_positive = sum(x[[method]] == "Doublet", na.rm = TRUE),
        method = method, stringsAsFactors = FALSE
      )))
      if (nrow(temp)) all_external[[length(all_external) + 1L]] <- temp
    }
    all_metrics[[study]] <- dat
    for (target in c("MG00", "MG0")) {
      pairs <- match_background(dat, target)
      if (nrow(pairs)) {
        target_rows <- dat[match(pairs$target_Cell, dat$Cell), ]; bg_rows <- dat[match(pairs$matched_MG_Cell, dat$Cell), ]
        for (metric in c("OL_independent_score", "MG_independent_score", "dual_lineage_score", "nearest_artificial_OL_MG_distance")) {
          delta <- target_rows[[metric]] - bg_rows[[metric]]; n <- sum(is.finite(delta)); if (!n) next
          all_effects[[length(all_effects) + 1L]] <- data.frame(Study = study, target_group = target, metric = metric, n_matched = n, effect_target_minus_matched_MG = mean(delta, na.rm = TRUE), ci_low = mean(delta, na.rm = TRUE) - 1.96 * stats::sd(delta, na.rm = TRUE) / sqrt(n), ci_high = mean(delta, na.rm = TRUE) + 1.96 * stats::sd(delta, na.rm = TRUE) / sqrt(n), stringsAsFactors = FALSE)
        }
        write.csv(pairs, file.path(table_dir, paste0(safe, "_", safe_name(target), "_to_MG_QC_match.csv")), row.names = FALSE)
      }
    }
    # DE is confined to real cells; overlapping MG00/MG0 cells are excluded from
    # the MG0 comparator to preserve a two-group contrast.
    de_cells <- dat$Cell[dat$cell_origin == "real_QC_passed" & (dat$MG00 %in% TRUE | (dat$MG0 %in% TRUE & !dat$MG00 %in% TRUE))]
    de_group <- ifelse(dat$MG00[match(de_cells, dat$Cell)] %in% TRUE, "MG00", "MG0_only")
    if (sum(de_group == "MG00") >= 10L && sum(de_group == "MG0_only") >= 10L) {
      de_object <- Seurat::CreateSeuratObject(counts = counts[, de_cells, drop = FALSE], project = study); de_object$comparison_group <- de_group
      de_object <- Seurat::NormalizeData(de_object, verbose = FALSE); Seurat::Idents(de_object) <- "comparison_group"
      de <- Seurat::FindMarkers(de_object, ident.1 = "MG00", ident.2 = "MG0_only", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0, verbose = FALSE)
      de$gene <- rownames(de); de$Study <- study; write.csv(de, file.path(de_dir, paste0(safe, "_real_MG00_vs_MG0_only.csv")), row.names = FALSE); all_de[[study]] <- de
    }
    rm(object, qc, counts, meta, dat, real); gc(verbose = FALSE)
  }
  metrics <- do.call(rbind, all_metrics); write.csv(metrics, file.path(table_dir, "per_cell_comparison_metrics.csv"), row.names = FALSE)
  external <- if (length(all_external)) do.call(rbind, all_external) else data.frame(); if (nrow(external)) { external$doublet_rate <- external$n_positive / external$n_available; write.csv(external, file.path(table_dir, "external_doublet_call_rates_real_MG00_MG0.csv"), row.names = FALSE) }
  effects <- if (length(all_effects)) do.call(rbind, all_effects) else data.frame(); if (nrow(effects)) write.csv(effects, file.path(table_dir, "study_batch_QC_matched_effects.csv"), row.names = FALSE)
  p_summary <- ggplot(metrics[metrics$comparison_group != "Other real QC-passed cell", ], aes(comparison_group, dual_lineage_score, fill = comparison_group)) + geom_boxplot(outlier.size = 0.18) + facet_wrap(~Study, scales = "free_x") + labs(title = "Six-study comparison: independent dual-lineage score", x = NULL, y = "min(OL score, MG score)") + theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot_both(p_summary, file.path(figure_dir, "six_study_dual_lineage_score_summary"), 15, 9); print(p_summary)
  if (nrow(external)) { p_external <- ggplot(external, aes(comparison_group, doublet_rate, fill = method)) + geom_col(position = "dodge") + facet_wrap(~Study, scales = "free_x") + labs(title = "Existing doublet calls: real MG0 and MG00", x = NULL, y = "Positive rate") + theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_text(hjust = 0.5, face = "bold")); save_plot_both(p_external, file.path(figure_dir, "six_study_external_doublet_call_rates"), 15, 8); print(p_external) }
  if (nrow(effects)) { p_forest <- ggplot(effects, aes(effect_target_minus_matched_MG, interaction(Study, target_group), colour = target_group)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") + geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2) + geom_point(size = 2) + facet_wrap(~metric, scales = "free_x") + labs(title = "Study/batch and QC-matched effects versus MG background", x = "Target minus matched MG", y = NULL) + theme_classic(base_size = 9) + theme(plot.title = element_text(hjust = 0.5, face = "bold")); save_plot_both(p_forest, file.path(figure_dir, "six_study_QC_matched_effect_forest"), 15, 9); print(p_forest) }
  figure_files <- do.call(rbind, all_figures); figure_files$png_exists <- file.exists(paste0(figure_files$stem, ".png")); figure_files$pdf_exists <- file.exists(paste0(figure_files$stem, ".pdf")); write.csv(figure_files, file.path(table_dir, "per_study_figure_manifest.csv"), row.names = FALSE)
  write.csv(data.frame(key = c("script", "run_timestamp", "OL_panel", "MG_panel", "distance_definition", "MG0_populations", "external_call_use", "DE_scope"), value = c("25.OL_MG_artificial_doublet_comparison_24server.r", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), paste(ol_markers, collapse = ";"), paste(mg_markers, collapse = ";"), "nearest parent-depth-matched artificial OL_MG in same-study/same-batch PCA reference space", "all QC-passed MG0 plus DoubletFinder-singlet sensitivity via DoubletFinder_call", "existing DoubletFinder/scDblFinder/Scrublet calls are descriptive external validation only", "real MG00 versus real MG0-only within each study; artificial profiles excluded"), stringsAsFactors = FALSE), file.path(runtime_dir, "stage3_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "sessionInfo.txt")); message("\nStage 3 complete. Output: ", result_root)
})
