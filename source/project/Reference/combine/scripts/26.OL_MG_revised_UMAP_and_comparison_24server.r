#!/usr/bin/env Rscript

# OL–MG artificial-reference revision.
#
# This is the authoritative terminology/display revision for the OL–MG branch.
# It does not rebuild counts, QC, doublet calls, parent pairings, PCA, or UMAP.
# Instead, it upgrades the six analysis objects from `canonical_MG_parent` to
# `MG_parent`, adds a derived `MG_other` comparison label, then redraws all
# required evidence using fixed UMAP coordinates and one highlighted category
# per figure.
#
# Key audit correction:
#   The original-pipeline 16.2% MG0 rate came from the original MG0–MG5 DoubletFinder
#   run on all 3,174 MG0 cells.  It cannot be compared to the
#   previously plotted post-selection reanalysis DoubletFinder rate.  This
#   script reports the two universes separately and uses scDblFinder/Scrublet
#   only as non-circular external callers for MG00.

local({
  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
    library(qs)
    library(ggplot2)
    library(RANN)
  })

  set.seed(20260910)
  project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
  combine_root <- file.path(project_root, "combine")
  legacy_stage2_root <- file.path(combine_root, "results/OL_MG_artificial_reference/02_artificial_profile_construction")
  result_root <- file.path(combine_root, "results/OL_MG_artificial_reference/04_revised_MG_terminology_and_comparison")
  object_dir <- file.path(result_root, "object")
  figure_dir <- file.path(result_root, "figure")
  table_dir <- file.path(result_root, "table")
  runtime_dir <- file.path(result_root, "runtime")
  for (dir_path in c(object_dir, figure_dir, table_dir, runtime_dir)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

  study_order <- c("Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021", "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024")
  safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  is_true <- function(x) !is.na(x) & trimws(tolower(as.character(x))) %in% c("true", "t", "1", "yes", "y")
  ol_markers <- c("PLP1", "MOBP", "MAG", "CNP", "CLDN11", "SOX10", "MYRF", "OPALIN")
  mg_markers <- c("P2RY12", "TMEM119", "CX3CR1", "GPR34", "CSF1R", "SPI1", "C1QA", "C1QB", "C1QC", "AIF1", "TYROBP", "FCER1G")
  if (!requireNamespace("RANN", quietly = TRUE)) stop("RANN is required for nearest artificial OL–MG distances.")

  save_plot_both <- function(plot, stem, width, height) {
    ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 320, bg = "white", limitsize = FALSE)
    # Cairo embeds UTF-8 glyphs in the PDF output. The base pdf() device on
    # this server writes the data correctly but cannot encode Chinese labels.
    ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
    if (!file.exists(paste0(stem, ".png")) || !file.exists(paste0(stem, ".pdf"))) stop("Incomplete PNG/PDF pair: ", stem)
  }
  get_counts <- function(object) {
    out <- tryCatch(SeuratObject::LayerData(object, assay = "RNA", layer = "counts"), error = function(e) NULL)
    if (is.null(out)) out <- tryCatch(SeuratObject::GetAssayData(object, assay = "RNA", layer = "counts"), error = function(e) NULL)
    if (is.null(out)) out <- SeuratObject::GetAssayData(object, assay = "RNA", slot = "counts")
    if (!inherits(out, "sparseMatrix")) out <- Matrix::Matrix(out, sparse = TRUE)
    out
  }
  assert_integer_counts <- function(x, label) {
    if (length(x@x) && any(abs(x@x - round(x@x)) > 1e-8)) stop(label, " contains non-integer counts.")
  }
  panel_score <- function(counts, markers) {
    present <- intersect(markers, rownames(counts))
    if (length(present) < 3L) stop("Too few panel markers available: ", paste(markers, collapse = ", "))
    sub <- counts[present, , drop = FALSE]
    lib <- Matrix::colSums(counts)
    norm <- log1p(sub %*% Matrix::Diagonal(x = 10000 / pmax(lib, 1)))
    list(score = as.numeric(Matrix::colMeans(norm)), present = present)
  }
  nearest_artificial_olmg <- function(meta, pcs) {
    output <- rep(NA_real_, nrow(meta))
    is_reference <- meta$cell_origin == "artificial_profile" & meta$profile_class == "OL_MG" & meta$construction == "parent_depth_matched"
    is_real <- meta$cell_origin == "real_QC_passed"
    for (batch in sort(unique(as.character(meta$batch_id[is_real])))) {
      real_rows <- which(is_real & as.character(meta$batch_id) == batch)
      reference_rows <- which(is_reference & as.character(meta$batch_id) == batch)
      if (!length(real_rows) || !length(reference_rows)) next
      output[real_rows] <- RANN::nn2(
        data = as.matrix(pcs[reference_rows, , drop = FALSE]),
        query = as.matrix(pcs[real_rows, , drop = FALSE]),
        k = 1L
      )$nn.dists[, 1L]
    }
    output
  }
  fixed_limits <- function(dat) {
    xr <- range(dat$UMAP_1, na.rm = TRUE); yr <- range(dat$UMAP_2, na.rm = TRUE)
    xp <- max(diff(xr) * 0.04, 0.25); yp <- max(diff(yr) * 0.04, 0.25)
    list(x = xr + c(-xp, xp), y = yr + c(-yp, yp))
  }
  single_highlight_umap <- function(dat, flag, highlight_label, colour, title, limits) {
    background <- dat[dat$cell_origin == "real_QC_passed", , drop = FALSE]
    foreground <- dat[flag, , drop = FALSE]
    if (!nrow(foreground)) stop(title, ": no cells/profile to highlight.")
    ggplot() +
      geom_point(data = background, aes(UMAP_1, UMAP_2), colour = "#D9D9D9", size = 0.25, alpha = 0.78) +
      geom_point(data = foreground, aes(UMAP_1, UMAP_2), colour = colour, size = 0.34, alpha = 0.95) +
      coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
      labs(
        title = title,
        subtitle = paste0("Single-label staining on the fixed real-cell reference UMAP; ", highlight_label, " (n=", nrow(foreground), ")"),
        x = "Reference UMAP 1", y = "Reference UMAP 2"
      ) +
      theme_classic(base_size = 10) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 8))
  }
  # Deliberately plain method diagrams.  They are analysis figures, saved in
  # both PNG and PDF, and can be placed in the companion presentation without
  # changing any biological result.
  simple_logic_plot <- function(title, subtitle, boxes, links, xlim = c(0, 10), ylim = c(0, 10)) {
    ggplot() +
      geom_rect(
        data = boxes,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE, fill = "#F7F7F7", colour = "#4D4D4D", linewidth = 0.45
      ) +
      geom_text(
        data = boxes,
        aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
        inherit.aes = FALSE, size = 4.0, lineheight = 0.95
      ) +
      geom_segment(
        data = links,
        aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE, linewidth = 0.45, colour = "#4D4D4D",
        arrow = grid::arrow(length = grid::unit(2.0, "mm"), type = "closed")
      ) +
      coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE, clip = "off") +
      labs(title = title, subtitle = subtitle) +
      theme_void(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 9, margin = margin(b = 8)),
        plot.margin = margin(12, 16, 12, 16)
      )
  }
  upgrade_mg_terminology <- function(object, study) {
    meta <- object@meta.data
    if (!"MG_parent" %in% names(meta)) {
      if (!"canonical_MG_parent" %in% names(meta)) stop(study, ": legacy object has neither MG_parent nor canonical_MG_parent.")
      meta$MG_parent <- as.logical(meta$canonical_MG_parent)
    }
    meta$canonical_MG_parent <- NULL
    for (field in intersect(c("parent_1_type", "parent_2_type"), names(meta))) {
      meta[[field]] <- gsub("canonical_MG", "MG", as.character(meta[[field]]), fixed = TRUE)
    }
    object@meta.data <- meta
    object
  }
  # MG0 and MG00 are complete labels, not mutually exclusive comparison bins:
  # a shared barcode remains in both label-specific summaries.  The reference
  # group is the reviewer-confirmed MG population after removing both labels.
  comparison_rows <- function(dat) {
    real <- dat$cell_origin == "real_QC_passed"
    mg0 <- dat$MG0 %in% TRUE
    mg00 <- dat$MG00 %in% TRUE
    specs <- list(
      "OL" = real & dat$mature_OL_parent %in% TRUE,
      "MG" = real & dat$MG_parent %in% TRUE,
      "MG other" = real & dat$MG_other %in% TRUE,
      "MG0" = real & mg0,
      "MG00" = real & mg00,
      "Artificial OL–MG | raw count sum" = dat$cell_origin == "artificial_profile" & dat$profile_class == "OL_MG" & dat$construction == "raw_count_sum",
      "Artificial OL–MG | parent-depth matched" = dat$cell_origin == "artificial_profile" & dat$profile_class == "OL_MG" & dat$construction == "parent_depth_matched"
    )
    out <- lapply(names(specs), function(label) {
      x <- dat[specs[[label]], , drop = FALSE]
      x$comparison_group <- label
      x
    })
    do.call(rbind, out)
  }

  historical_file <- file.path(combine_root, "results/MG0_MG5_doublet_diagnosis/table/six_study_cell_metadata_with_all_doublet_calls.csv")
  if (!file.exists(historical_file)) stop("Original MG0–MG5 table is missing: ", historical_file)
  historical <- read.csv(historical_file, check.names = FALSE, stringsAsFactors = FALSE)
  required_historical <- c("Study", "Cell", "CellType1", "DoubletFinder_call")
  if (length(setdiff(required_historical, names(historical)))) stop("Original MG0 table missing required columns.")
  historical$key <- paste(historical$Study, historical$Cell, sep = "||")
  if (anyDuplicated(historical$key)) stop("Historical Study/barcode keys are duplicated.")

  group_levels <- c(
    "OL", "MG", "MG other", "MG0", "MG00",
    "Artificial OL–MG | raw count sum", "Artificial OL–MG | parent-depth matched"
  )
  group_colours <- c(
    "OL" = "#D7301F", "MG" = "#006D2C", "MG other" = "#1B9E77",
    "MG0" = "#2C7FB8", "MG00" = "#C51B7D",
    "Artificial OL–MG | raw count sum" = "#FC8D59", "Artificial OL–MG | parent-depth matched" = "#F03B20"
  )
  # Every comparison plot used in the presentation is generated once per
  # study and once per metric.  The cross-study faceted summaries below are
  # retained only as compact supplementary overviews, never as the source for
  # a presentation slide.  This keeps every displayed analysis result a
  # complete, uncropped R-generated single figure.
  single_comparison_boxplot <- function(metric_dat, metric_column, y_label, title, subtitle) {
    metric_dat$comparison_group <- factor(as.character(metric_dat$comparison_group), levels = group_levels)
    primary_dat <- metric_dat[metric_dat$comparison_group == "Artificial OL–MG | parent-depth matched", , drop = FALSE]
    ggplot(metric_dat, aes_string(x = "comparison_group", y = metric_column, fill = "comparison_group")) +
      geom_boxplot(outlier.size = 0.22, na.rm = TRUE) +
      geom_boxplot(data = primary_dat, fill = NA, colour = "#B2182B", linewidth = 1.0, outlier.shape = NA, na.rm = TRUE) +
      scale_fill_manual(values = group_colours, drop = FALSE) +
      labs(title = title, subtitle = subtitle, x = NULL, y = y_label) +
      theme_classic(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 8)
      )
  }
  all_metrics <- list(); all_depth <- list(); all_external <- list(); all_current_df <- list(); all_mapping <- list(); figure_rows <- list()

  # 1) Overall analysis logic.
  p_logic_overall <- simple_logic_plot(
    "OL–MG 人工双细胞参照与 MG0/MG00 比较：总体逻辑",
    "所有比较均在研究内完成；人工 profile 不参与真实细胞的 QC、聚类或标签重定义。",
    data.frame(
      xmin = 1, xmax = 9, ymin = c(8.3, 6.5, 4.7, 2.9, 1.1), ymax = c(9.4, 7.6, 5.8, 4.0, 2.2),
      label = c(
        "用户确认成熟 OL 簇与 MG 簇\n每研究、每技术 batch 内选取 DoubletFinder singlet 父细胞",
        "OL 整数 count + MG 整数 count\n构建人工 OL–MG 原始 profile",
        "下采样至两父细胞平均 UMI\n深度匹配人工 OL–MG 作为主参照",
        "真实细胞先拟合 PCA / UMAP\n人工 profile 仅投射到固定参考空间",
        "比较 OL、MG、MG other、MG0、MG00\n并与人工 OL–MG（原始相加与深度匹配）参照比较\n双谱系分数、PCA 最近距离、外部双细胞调用"
      ), stringsAsFactors = FALSE
    ),
    data.frame(x = 5, y = c(8.3, 6.5, 4.7, 2.9), xend = 5, yend = c(7.6, 5.8, 4.0, 2.2))
  )
  overall_stem <- file.path(figure_dir, "method_logic_overall_OL_MG_artificial_reference")
  save_plot_both(p_logic_overall, overall_stem, 9.5, 8); print(p_logic_overall)
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = "method", figure = "logic_overall", stem = overall_stem, stringsAsFactors = FALSE)

  # 2) Count construction and depth matching logic.
  p_logic_depth <- simple_logic_plot(
    "人工 OL–MG profile：原始 count 相加与深度匹配",
    "两种 profile 均保留整数 count；主分析使用深度匹配版以减少总 UMI 翻倍的影响。",
    data.frame(
      xmin = c(0.4, 2.7, 5.0, 7.3), xmax = c(2.1, 4.4, 6.7, 9.6), ymin = 3.7, ymax = 6.4,
      label = c(
        "OL 父细胞\n整数向量 cOL\nUMIOL",
        "MG 父细胞\n整数向量 cMG\nUMIMG",
        "原始相加 profile\ncOL + cMG\n总 UMI = UMIOL + UMIMG",
        "深度匹配 profile（主参照）\n对相加 profile 多项分布下采样\n目标 UMI = floor[(UMIOL + UMIMG) / 2]"
      ), stringsAsFactors = FALSE
    ),
    data.frame(x = c(2.1, 4.4, 6.7), y = c(5.1, 5.1, 5.1), xend = c(2.7, 5.0, 7.3), yend = c(5.1, 5.1, 5.1)),
    xlim = c(0, 10), ylim = c(2.3, 7.6)
  )
  depth_logic_stem <- file.path(figure_dir, "method_logic_artificial_OL_MG_depth_matching")
  save_plot_both(p_logic_depth, depth_logic_stem, 12, 5.5); print(p_logic_depth)
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = "method", figure = "logic_artificial_depth_matching", stem = depth_logic_stem, stringsAsFactors = FALSE)

  # 3) Nearest-reference PCA distance logic.
  p_logic_distance <- simple_logic_plot(
    "PCA 最近距离：真实细胞与主人工 OL–MG 参照的相似性",
    "距离在每研究、同一 batch、固定真实细胞 PCA 的前 20 个维度中计算；UMAP 只用于可视化。",
    data.frame(
      xmin = 1, xmax = 9, ymin = c(7.8, 5.9, 4.0, 2.1), ymax = c(8.9, 7.0, 5.1, 3.2),
      label = c(
        "真实 QC 通过细胞：标准化、缩放、PCA、UMAP\n得到固定的真实细胞参考空间",
        "深度匹配人工 OL–MG：使用同一基因空间\n投射到已拟合 PCA / UMAP，不重新拟合真实细胞流形",
        "对每个真实细胞 i：在同一 batch 的人工 profiles 中\n寻找前 20 个 PCA 维度的最近邻 j",
        "d(i) = minj ||PCi - PCj||\n距离较小表示更接近人工 OL–MG；人工 profile 自身距离为 0"
      ), stringsAsFactors = FALSE
    ),
    data.frame(x = 5, y = c(7.8, 5.9, 4.0), xend = 5, yend = c(7.0, 5.1, 3.2))
  )
  distance_logic_stem <- file.path(figure_dir, "method_logic_PCA_nearest_artificial_OL_MG_distance")
  save_plot_both(p_logic_distance, distance_logic_stem, 10, 7.5); print(p_logic_distance)
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = "method", figure = "logic_PCA_nearest_artificial_OL_MG_distance", stem = distance_logic_stem, stringsAsFactors = FALSE)

  # 4) Independent OL/MG scoring logic.  The two panels do not reuse the six
  # markers that selected MG00, avoiding a circular score construction.
  p_logic_score <- simple_logic_plot(
    "独立 OL score、MG score 与双谱系分数",
    "每个 score 使用原始整数 count 的 library-size 标准化表达；双谱系分数要求两个谱系信号同时存在。",
    data.frame(
      xmin = 0.8, xmax = 9.2,
      ymin = c(8.0, 6.25, 3.25, 1.15), ymax = c(9.0, 7.25, 5.45, 2.35),
      label = c(
        "每个细胞 j 的整数 count cgj；L_j = 总 UMI",
        "对每个基因标准化：xgj = log1p(10,000 × cgj / L_j)",
        "OL score = OL panel 的平均 xgj\nOL panel：PLP1、MOBP、MAG、CNP、CLDN11、SOX10、MYRF、OPALIN\n\nMG score = MG panel 的平均 xgj\nMG panel：P2RY12、TMEM119、CX3CR1、GPR34、CSF1R、SPI1、C1QA/B/C、AIF1、TYROBP、FCER1G",
        "双谱系分数 = min(OL score, MG score)\n只有 OL 和 MG 两种独立信号均较高时，双谱系分数才较高"
      ), stringsAsFactors = FALSE
    ),
    data.frame(x = 5, y = c(8.0, 6.25, 3.25), xend = 5, yend = c(7.25, 5.45, 2.35)),
    xlim = c(0, 10), ylim = c(0.5, 9.6)
  )
  score_logic_stem <- file.path(figure_dir, "method_logic_independent_OL_MG_and_dual_scores")
  save_plot_both(p_logic_score, score_logic_stem, 11, 8); print(p_logic_score)
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = "method", figure = "logic_independent_OL_MG_and_dual_scores", stem = score_logic_stem, stringsAsFactors = FALSE)

  for (study in study_order) {
    safe <- safe_name(study)
    message("\n===== Revised OL–MG display and audit: ", study, " =====")
    input_file <- file.path(legacy_stage2_root, "object", paste0(safe, "_real_cells_plus_artificial_profiles.qs"))
    if (!file.exists(input_file)) stop(study, ": missing legacy artificial-reference object.")
    object <- upgrade_mg_terminology(qs::qread(input_file), study)
    counts <- get_counts(object); assert_integer_counts(counts, study)
    meta <- object[[]]
    meta$Cell <- rownames(meta)
    if (!"MG0" %in% names(meta)) {
      if (!"historical_MG0" %in% names(meta)) stop(study, ": object has neither MG0 nor historical_MG0 metadata.")
      meta$MG0 <- as.logical(meta$historical_MG0)
    }
    meta$historical_MG0 <- NULL
    required_meta <- c("cell_origin", "profile_class", "construction", "MG0", "MG00", "mature_OL_parent", "MG_parent", "batch_id", "original_barcode", "DoubletFinder_call", "scDblFinder_call", "Scrublet_call")
    if (length(setdiff(required_meta, names(meta)))) stop(study, ": upgraded object lacks required metadata: ", paste(setdiff(required_meta, names(meta)), collapse = ", "))

    key <- paste(study, as.character(meta$original_barcode), sep = "||")
    hidx <- match(key, historical$key)
    meta$MG0 <- !is.na(hidx) & historical$CellType1[hidx] == "MG0"
    meta$original_pipeline_DoubletFinder_call <- ifelse(!is.na(hidx), as.character(historical$DoubletFinder_call[hidx]), NA_character_)
    meta$MG0_only <- meta$MG0 %in% TRUE & !(meta$MG00 %in% TRUE)
    meta$MG00_only <- !(meta$MG0 %in% TRUE) & meta$MG00 %in% TRUE
    meta$MG0_and_MG00 <- meta$MG0 %in% TRUE & meta$MG00 %in% TRUE
    meta$MG_other <- meta$cell_origin == "real_QC_passed" & meta$MG_parent %in% TRUE & !(meta$MG0 %in% TRUE) & !(meta$MG00 %in% TRUE)

    # Keep terminology and the new MG_other label in memory while scores and
    # distances are calculated.  The authoritative object is saved below with
    # these derived fields included; counts, projections and pairings remain
    # unchanged.
    object@meta.data <- meta[rownames(object@meta.data), setdiff(names(meta), "Cell"), drop = FALSE]

    ol_score <- panel_score(counts, ol_markers)$score
    mg_score <- panel_score(counts, mg_markers)$score
    meta$OL_independent_score <- ol_score
    meta$MG_independent_score <- mg_score
    meta$dual_lineage_score <- pmin(ol_score, mg_score)
    pcs <- Seurat::Embeddings(object, "pca_real_reference_projection")
    pcs <- pcs[meta$Cell, seq_len(min(20L, ncol(pcs))), drop = FALSE]
    meta$nearest_artificial_OL_MG_distance <- nearest_artificial_olmg(meta, pcs)
    object@meta.data <- meta[rownames(object@meta.data), setdiff(names(meta), "Cell"), drop = FALSE]
    qs::qsave(object, file.path(object_dir, paste0(safe, "_real_cells_plus_artificial_profiles_MG_terminology.qs")), preset = "high")
    umap <- as.data.frame(Seurat::Embeddings(object, "umap_real_reference_projection")[meta$Cell, , drop = FALSE])
    names(umap) <- c("UMAP_1", "UMAP_2")
    dat <- cbind(meta, umap)
    limits <- fixed_limits(dat)

    # Every figure below uses the same stored UMAP coordinates and limits, with
    # exactly one highlighted label.  This is the project-wide default.
    umap_specs <- list(
      list(id = "MG0", flag = dat$cell_origin == "real_QC_passed" & dat$MG0 %in% TRUE, label = "MG0 (all mapped cells; overlap with MG00 is retained)", colour = "#2C7FB8"),
      list(id = "MG00", flag = dat$cell_origin == "real_QC_passed" & dat$MG00 %in% TRUE, label = "MG00", colour = "#C51B7D"),
      list(id = "Mature_OL", flag = dat$cell_origin == "real_QC_passed" & dat$mature_OL_parent %in% TRUE, label = "reviewer-confirmed mature OL parent", colour = "#D7301F"),
      list(id = "MG", flag = dat$cell_origin == "real_QC_passed" & dat$MG_parent %in% TRUE, label = "reviewer-confirmed MG parent", colour = "#1B9E77"),
      list(id = "Artificial_OL_MG_parent_depth_matched", flag = dat$cell_origin == "artificial_profile" & dat$profile_class == "OL_MG" & dat$construction == "parent_depth_matched", label = "primary artificial OL–MG reference: parent-depth-matched", colour = "#F03B20"),
      list(id = "Artificial_OL_MG_raw_count_sum", flag = dat$cell_origin == "artificial_profile" & dat$profile_class == "OL_MG" & dat$construction == "raw_count_sum", label = "artificial OL–MG: raw-count-sum sensitivity reference", colour = "#FC8D59")
    )
    for (spec in umap_specs) {
      plot <- single_highlight_umap(dat, spec$flag, spec$label, spec$colour, paste0(study, ": ", spec$label), limits)
      stem <- file.path(figure_dir, paste0(safe, "_single_label_UMAP_", spec$id))
      save_plot_both(plot, stem, 7.5, 6.2); print(plot)
      figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = study, figure = paste0("single_label_UMAP_", spec$id), stem = stem, stringsAsFactors = FALSE)
    }

    comparison_dat <- comparison_rows(dat)
    comparison_dat$comparison_group <- factor(comparison_dat$comparison_group, levels = group_levels)
    if (!nrow(comparison_dat)) stop(study, ": no OL–MG comparison cells/profiles.")
    all_metrics[[study]] <- comparison_dat

    # Per-study, per-metric single figures.  These are the presentation
    # sources; each is a full R plot, never a crop or a panel extracted from a
    # six-study facet figure.
    score_specs <- list(
      list(id = "independent_OL_score", column = "OL_independent_score", y = "Independent OL score", label = "independent OL score"),
      list(id = "independent_MG_score", column = "MG_independent_score", y = "Independent MG score", label = "independent MG score"),
      list(id = "dual_lineage_score", column = "dual_lineage_score", y = "Dual-lineage score: min(OL score, MG score)", label = "dual-lineage score")
    )
    for (score_spec in score_specs) {
      score_plot <- single_comparison_boxplot(
        comparison_dat, score_spec$column, score_spec$y,
        paste0(study, ": ", score_spec$label),
        "OL, MG, MG other, MG0, MG00, and both artificial OL–MG constructions. Red outline = primary parent-depth-matched artificial OL–MG reference."
      )
      score_stem <- file.path(figure_dir, paste0(safe, "_", score_spec$id))
      save_plot_both(score_plot, score_stem, 10, 6.6); print(score_plot)
      figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = study, figure = score_spec$id, stem = score_stem, stringsAsFactors = FALSE)
    }

    # Distances are defined for real cells only.  A depth-matched artificial
    # profile's self-distance is identically zero, so it is shown as a labelled
    # reference point rather than a misleading boxplot distribution.
    distance_dat <- comparison_dat[comparison_dat$cell_origin == "real_QC_passed" & comparison_dat$comparison_group %in% group_levels[seq_len(5L)], , drop = FALSE]
    distance_dat$distance_group <- factor(as.character(distance_dat$comparison_group), levels = group_levels[seq_len(5L)])
    zero_reference <- data.frame(
      Study = study,
      distance_group = factor("Artificial OL–MG self-reference", levels = c(group_levels[seq_len(5L)], "Artificial OL–MG self-reference")),
      nearest_artificial_OL_MG_distance = 0
    )
    distance_dat$distance_group <- factor(as.character(distance_dat$distance_group), levels = levels(zero_reference$distance_group))
    p_distance <- ggplot(distance_dat, aes(distance_group, nearest_artificial_OL_MG_distance, fill = distance_group)) +
      geom_boxplot(outlier.size = 0.25, na.rm = TRUE) +
      geom_hline(yintercept = 0, linetype = 2, colour = "#F03B20") +
      geom_point(data = zero_reference, aes(distance_group, nearest_artificial_OL_MG_distance), inherit.aes = FALSE, shape = 8, size = 3, colour = "#F03B20") +
      scale_fill_manual(values = c(group_colours[group_levels[seq_len(5L)]], "Artificial OL–MG self-reference" = "white"), drop = FALSE) +
      scale_x_discrete(labels = c(
        "OL" = "OL", "MG" = "MG", "MG other" = "MG other", "MG0" = "MG0", "MG00" = "MG00",
        "Artificial OL–MG self-reference" = "Artificial reference\nself-distance = 0"
      )) +
      labs(
        title = paste0(study, ": PCA distance to the primary artificial OL–MG reference"),
        subtitle = "Fixed 20-PC space, same batch. Lower distance = closer to parent-depth-matched artificial OL–MG; red star = artificial self-distance 0.",
        x = NULL, y = "Nearest PCA-space distance"
      ) +
      theme_classic(base_size = 9) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 7))
    distance_stem <- file.path(figure_dir, paste0(safe, "_distance_to_primary_artificial_OL_MG"))
    save_plot_both(p_distance, distance_stem, 10, 7.2); print(p_distance)
    figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = study, figure = "distance_to_primary_artificial_OL_MG", stem = distance_stem, stringsAsFactors = FALSE)

    depth_rows <- dat[dat$cell_origin == "artificial_profile" & dat$profile_class == "OL_MG", c("Study", "construction", "nCount_RNA", "parent_1_umi", "parent_2_umi", "raw_profile_umi", "depth_matched_target_umi")]
    depth_rows$construction <- factor(depth_rows$construction, levels = c("raw_count_sum", "parent_depth_matched"), labels = c("Raw count sum", "Parent-depth matched (primary)"))
    all_depth[[study]] <- depth_rows
    p_depth_study <- ggplot(depth_rows, aes(construction, nCount_RNA, fill = construction)) +
      geom_boxplot(outlier.size = 0.24) +
      scale_fill_manual(values = c("Raw count sum" = "#FC8D59", "Parent-depth matched (primary)" = "#F03B20")) +
      labs(
        title = paste0(study, ": artificial OL–MG depth audit"),
        subtitle = "Raw count sum = OL parent count + MG parent count. Primary depth-matched profile = multinomial downsample of the summed count to the mean parent UMI depth.",
        x = NULL, y = "Artificial profile total UMI"
      ) +
      theme_classic(base_size = 10) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 8))
    depth_study_stem <- file.path(figure_dir, paste0(safe, "_artificial_OL_MG_depth_audit"))
    save_plot_both(p_depth_study, depth_study_stem, 9, 6); print(p_depth_study)
    figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = study, figure = "artificial_OL_MG_depth_audit", stem = depth_study_stem, stringsAsFactors = FALSE)

    # Non-circular external methods.  Reanalysis DoubletFinder is deliberately
    # excluded here because MG00 was selected after that exact call removed
    # doublets; it is reported in a separate conditioning audit below.
    ext_groups <- c("MG", "MG other", "MG0", "MG00")
    external_study_rows <- list()
    for (method in c("scDblFinder_call", "Scrublet_call")) {
      for (group in ext_groups) {
        current <- comparison_dat[comparison_dat$cell_origin == "real_QC_passed" & comparison_dat$comparison_group == group, , drop = FALSE]
        available <- !is.na(current[[method]])
        external_row <- data.frame(
          Study = study, group = group, method = method, n_available = sum(available),
          n_doublet = sum(current[[method]][available] == "Doublet"),
          doublet_rate = if (sum(available)) sum(current[[method]][available] == "Doublet") / sum(available) else NA_real_,
          stringsAsFactors = FALSE
        )
        all_external[[length(all_external) + 1L]] <- external_row
        external_study_rows[[length(external_study_rows) + 1L]] <- external_row
      }
    }
    external_study <- do.call(rbind, external_study_rows)
    external_study$group <- factor(external_study$group, levels = ext_groups)
    p_external_study <- ggplot(external_study, aes(group, doublet_rate, fill = method)) +
      geom_col(position = position_dodge(width = 0.78), na.rm = TRUE) +
      geom_text(aes(label = ifelse(is.na(doublet_rate), "NA", sprintf("n=%d", n_available))), position = position_dodge(width = 0.78), vjust = -0.18, size = 2.6, na.rm = TRUE) +
      scale_fill_manual(values = c(scDblFinder_call = "#00A651", Scrublet_call = "#5B8FF9")) +
      scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = expansion(mult = c(0, 0.18))) +
      labs(title = paste0(study, ": non-circular external doublet evidence"), subtitle = "scDblFinder and Scrublet only. MG00 was selected after reanalysis DoubletFinder removal, so that call is not independent validation.", x = NULL, y = "Positive rate") +
      theme_classic(base_size = 10) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 8))
    external_study_stem <- file.path(figure_dir, paste0(safe, "_external_scDblFinder_Scrublet_non_circular_rates"))
    save_plot_both(p_external_study, external_study_stem, 9, 6); print(p_external_study)
    figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = study, figure = "external_scDblFinder_Scrublet_non_circular_rates", stem = external_study_stem, stringsAsFactors = FALSE)

    conditioning_study_rows <- list()
    for (group in ext_groups) {
      current <- comparison_dat[comparison_dat$cell_origin == "real_QC_passed" & comparison_dat$comparison_group == group, , drop = FALSE]
      conditioning_row <- data.frame(
        Study = study, group = group, n_available = sum(!is.na(current$DoubletFinder_call)),
        n_doublet = sum(current$DoubletFinder_call == "Doublet", na.rm = TRUE),
        current_reanalysis_DoubletFinder_rate = if (sum(!is.na(current$DoubletFinder_call))) sum(current$DoubletFinder_call == "Doublet", na.rm = TRUE) / sum(!is.na(current$DoubletFinder_call)) else NA_real_,
        stringsAsFactors = FALSE
      )
      all_current_df[[length(all_current_df) + 1L]] <- conditioning_row
      conditioning_study_rows[[length(conditioning_study_rows) + 1L]] <- conditioning_row
    }
    conditioning_study <- do.call(rbind, conditioning_study_rows)
    conditioning_study$group <- factor(conditioning_study$group, levels = ext_groups)
    p_conditioning_study <- ggplot(conditioning_study, aes(group, current_reanalysis_DoubletFinder_rate, fill = group)) +
      geom_col(na.rm = TRUE) +
      geom_text(aes(label = ifelse(is.na(current_reanalysis_DoubletFinder_rate), "NA", sprintf("%.1f%%\\n(n=%d)", 100 * current_reanalysis_DoubletFinder_rate, n_available))), vjust = -0.18, size = 2.6, na.rm = TRUE) +
      scale_fill_manual(values = group_colours[ext_groups]) +
      scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = expansion(mult = c(0, 0.20))) +
      labs(title = paste0(study, ": DoubletFinder conditioning audit"), subtitle = "MG00 was assigned after this exact reanalysis DoubletFinder removal. Its low/zero rate is structural and is not independent MG00 validation.", x = NULL, y = "Reanalysis DoubletFinder positive rate") +
      theme_classic(base_size = 10) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 8))
    conditioning_study_stem <- file.path(figure_dir, paste0(safe, "_DoubletFinder_conditioning_audit"))
    save_plot_both(p_conditioning_study, conditioning_study_stem, 9, 6); print(p_conditioning_study)
    figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = study, figure = "DoubletFinder_conditioning_audit", stem = conditioning_study_stem, stringsAsFactors = FALSE)
    historical_mg0_total <- historical[historical$Study == study & historical$CellType1 == "MG0", , drop = FALSE]
    all_mapping[[length(all_mapping) + 1L]] <- data.frame(
      Study = study,
      MG0_original_total = nrow(historical_mg0_total),
      MG0_original_DoubletFinder_doublets = sum(historical_mg0_total$DoubletFinder_call == "Doublet"),
      MG0_original_DoubletFinder_rate = mean(historical_mg0_total$DoubletFinder_call == "Doublet"),
      MG0_mapped_to_reanalysis = sum(dat$cell_origin == "real_QC_passed" & dat$MG0 %in% TRUE),
      MG00_mapped_to_reanalysis = sum(dat$cell_origin == "real_QC_passed" & dat$MG00 %in% TRUE),
      MG0_and_MG00_audit_only = sum(dat$cell_origin == "real_QC_passed" & dat$MG0_and_MG00 %in% TRUE),
      MG_other = sum(dat$cell_origin == "real_QC_passed" & dat$MG_other %in% TRUE),
      stringsAsFactors = FALSE
    )
    rm(object, counts, meta, pcs, umap, dat, comparison_dat, distance_dat); gc(verbose = FALSE)
  }

  metrics <- do.call(rbind, all_metrics)
  depth <- do.call(rbind, all_depth)
  external <- do.call(rbind, all_external)
  current_df <- do.call(rbind, all_current_df)
  mapping <- do.call(rbind, all_mapping)
  write.csv(metrics, file.path(table_dir, "revised_per_cell_metrics.csv"), row.names = FALSE)
  write.csv(depth, file.path(table_dir, "artificial_OL_MG_depth_audit.csv"), row.names = FALSE)
  write.csv(external, file.path(table_dir, "external_scDblFinder_Scrublet_non_circular_rates.csv"), row.names = FALSE)
  write.csv(current_df, file.path(table_dir, "current_reanalysis_DoubletFinder_conditioning_audit.csv"), row.names = FALSE)
  write.csv(mapping, file.path(table_dir, "MG0_denominator_and_mapping_audit.csv"), row.names = FALSE)

  primary_group <- "Artificial OL–MG | parent-depth matched"
  # Six-study view of the PCA-distance result.  Artificial profiles have a
  # defined self-distance of zero; they are shown as stars rather than mixed
  # into the real-cell boxplots.
  distance_real_groups <- c("OL", "MG", "MG other", "MG0", "MG00")
  distance_plot_dat <- metrics[metrics$comparison_group %in% distance_real_groups, , drop = FALSE]
  distance_zero <- data.frame(
    Study = factor(study_order, levels = study_order),
    comparison_group = factor(rep(primary_group, length(study_order)), levels = group_levels),
    nearest_artificial_OL_MG_distance = 0
  )
  p_distance_summary <- ggplot(distance_plot_dat, aes(comparison_group, nearest_artificial_OL_MG_distance, fill = comparison_group)) +
    geom_boxplot(outlier.size = 0.13, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = 2, colour = "#F03B20") +
    geom_point(data = distance_zero, aes(comparison_group, nearest_artificial_OL_MG_distance), inherit.aes = FALSE, shape = 8, size = 3.0, colour = "#F03B20") +
    scale_fill_manual(values = c(group_colours[distance_real_groups], "Artificial OL–MG | parent-depth matched" = "white"), drop = FALSE) +
    facet_wrap(~Study, scales = "free_y") +
    labs(
      title = "Six-study PCA distance to the primary artificial OL–MG reference",
      subtitle = "For every real cell: nearest same-batch parent-depth-matched artificial OL–MG in the first 20 fixed reference PCs. Lower = more similar. Red star at 0 = artificial OL–MG self-distance. MG other = MG excluding MG0 and MG00; MG0/MG00 retain shared cells.",
      x = NULL, y = "Nearest PCA-space distance"
    ) +
    theme_classic(base_size = 8.5) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 7))
  distance_summary_stem <- file.path(figure_dir, "six_study_PCA_distance_to_primary_artificial_OL_MG_with_MGother")
  save_plot_both(p_distance_summary, distance_summary_stem, 16, 9); print(p_distance_summary)
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = "six-study", figure = "PCA_distance_to_primary_artificial_OL_MG_with_MGother", stem = distance_summary_stem, stringsAsFactors = FALSE)

  score_long <- rbind(
    data.frame(metrics, score_name = "OL score", score_value = metrics$OL_independent_score, stringsAsFactors = FALSE),
    data.frame(metrics, score_name = "MG score", score_value = metrics$MG_independent_score, stringsAsFactors = FALSE)
  )
  p_independent_scores <- ggplot(score_long, aes(comparison_group, score_value, fill = comparison_group)) +
    geom_boxplot(outlier.size = 0.12, na.rm = TRUE) +
    geom_boxplot(data = score_long[score_long$comparison_group == primary_group, , drop = FALSE], fill = NA, colour = "#B2182B", linewidth = 0.9, outlier.shape = NA) +
    scale_fill_manual(values = group_colours, drop = FALSE) +
    facet_grid(score_name ~ Study, scales = "free_y") +
    labs(
      title = "Six-study independent OL and MG scores",
      subtitle = "MG other = confirmed MG after excluding MG0 and MG00. MG0 and MG00 each retain shared cells. Red outline = primary parent-depth-matched artificial OL–MG reference.",
      x = NULL, y = "Mean log-normalized expression across the indicated marker panel"
    ) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 6.8))
  independent_stem <- file.path(figure_dir, "six_study_independent_OL_and_MG_scores_with_MGother")
  save_plot_both(p_independent_scores, independent_stem, 18, 10); print(p_independent_scores)
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = "six-study", figure = "independent_OL_and_MG_scores_with_MGother", stem = independent_stem, stringsAsFactors = FALSE)

  p_dual <- ggplot(metrics, aes(comparison_group, dual_lineage_score, fill = comparison_group)) +
    geom_boxplot(outlier.size = 0.15, na.rm = TRUE) +
    geom_boxplot(data = metrics[metrics$comparison_group == primary_group, , drop = FALSE], fill = NA, colour = "#B2182B", linewidth = 0.9, outlier.shape = NA) +
    scale_fill_manual(values = group_colours, drop = FALSE) +
    facet_wrap(~Study, scales = "free_x") +
    labs(
      title = "Six-study dual-lineage comparison: OL, MG, MG other, MG0, MG00 and artificial OL–MG",
      subtitle = "MG other = confirmed MG after excluding MG0 and MG00. MG0 and MG00 each retain shared cells. Red outline = primary OL–MG reference (parent-depth-matched).",
      x = NULL, y = "min(independent OL score, independent MG score)"
    ) +
    theme_classic(base_size = 8.5) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 7))
  dual_stem <- file.path(figure_dir, "six_study_dual_lineage_score_with_MGother_and_primary_reference")
  save_plot_both(p_dual, dual_stem, 16, 9); print(p_dual)
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = "six-study", figure = "dual_lineage_score_with_MGother_and_primary_reference", stem = dual_stem, stringsAsFactors = FALSE)

  p_depth <- ggplot(depth, aes(construction, nCount_RNA, fill = construction)) +
    geom_boxplot(outlier.size = 0.18) +
    facet_wrap(~Study, scales = "free_y") +
    scale_fill_manual(values = c("Raw count sum" = "#FC8D59", "Parent-depth matched (primary)" = "#F03B20")) +
    labs(
      title = "Artificial OL–MG depth audit",
      subtitle = "Raw count sum = parent OL counts + parent MG counts. Primary depth-matched profile = multinomial downsample of that sum to the mean UMI depth of its two parents.",
      x = NULL, y = "Artificial profile total UMI"
    ) +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 7))
  depth_stem <- file.path(figure_dir, "six_study_artificial_OL_MG_depth_audit")
  save_plot_both(p_depth, depth_stem, 15, 8); print(p_depth)
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = "six-study", figure = "artificial_OL_MG_depth_audit", stem = depth_stem, stringsAsFactors = FALSE)

  historical_states <- historical[historical$CellType1 %in% paste0("MG", 0:5), , drop = FALSE]
  historical_df <- aggregate(as.integer(historical_states$DoubletFinder_call == "Doublet"), by = list(CellType1 = historical_states$CellType1), FUN = sum)
  names(historical_df)[2] <- "n_doublet"
  historical_n <- aggregate(Cell ~ CellType1, data = historical_states, FUN = length); names(historical_n)[2] <- "n_cells"
  historical_df <- merge(historical_df, historical_n, by = "CellType1")
  historical_df$rate <- historical_df$n_doublet / historical_df$n_cells
  historical_df$CellType1 <- factor(historical_df$CellType1, levels = paste0("MG", 0:5))
  write.csv(historical_df, file.path(table_dir, "original_MG0_MG5_DoubletFinder_rates.csv"), row.names = FALSE)
  p_historical <- ggplot(historical_df, aes(CellType1, rate, fill = CellType1)) +
    geom_col(width = 0.72, colour = "white") +
    geom_text(aes(label = sprintf("%.1f%%\\n(n=%d)", 100 * rate, n_cells)), vjust = -0.18, size = 3) +
    scale_fill_manual(values = c(MG0 = "#8DD3C7", MG1 = "#FFFFB3", MG2 = "#BEBADA", MG3 = "#FB8072", MG4 = "#80B1D3", MG5 = "#FDB462")) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Original pipeline: DoubletFinder proportion across MG0–MG5", subtitle = "Original 20,958-cell MG0–MG5 diagnostic universe; MG0 is 514/3,174 = 16.2%. This is not the post-selection reanalysis denominator.", x = "MG state", y = "Doublet proportion") +
    theme_classic(base_size = 10) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 8))
  historical_stem <- file.path(figure_dir, "original_MG0_MG5_DoubletFinder_rates")
  save_plot_both(p_historical, historical_stem, 8, 5); print(p_historical)
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = "six-study", figure = "original_MG0_MG5_DoubletFinder_rates", stem = historical_stem, stringsAsFactors = FALSE)

  p_external <- ggplot(external, aes(group, doublet_rate, fill = method)) +
    geom_col(position = position_dodge(width = 0.78), na.rm = TRUE) +
    geom_text(aes(label = ifelse(is.na(doublet_rate), "NA", sprintf("n=%d", n_available))), position = position_dodge(width = 0.78), vjust = -0.18, size = 2.2, na.rm = TRUE) +
    facet_wrap(~Study, scales = "free_x") +
    scale_fill_manual(values = c(scDblFinder_call = "#00A651", Scrublet_call = "#5B8FF9")) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Non-circular external doublet evidence: scDblFinder and Scrublet", subtitle = "MG00 was defined after reanalysis DoubletFinder removal, so DoubletFinder is intentionally excluded. MG other = MG after excluding MG0 and MG00; MG0 and MG00 each retain shared cells; n is the available-call denominator.", x = NULL, y = "Positive rate") +
    theme_classic(base_size = 8.5) + theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 7))
  external_stem <- file.path(figure_dir, "external_scDblFinder_Scrublet_non_circular_rates")
  save_plot_both(p_external, external_stem, 16, 8); print(p_external)
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = "six-study", figure = "external_scDblFinder_Scrublet_non_circular_rates", stem = external_stem, stringsAsFactors = FALSE)

  p_conditioning <- ggplot(current_df, aes(group, current_reanalysis_DoubletFinder_rate, fill = group)) +
    geom_col(na.rm = TRUE) +
    geom_text(aes(label = ifelse(is.na(current_reanalysis_DoubletFinder_rate), "NA", sprintf("%.1f%%\\n(n=%d)", 100 * current_reanalysis_DoubletFinder_rate, n_available))), vjust = -0.18, size = 2.2, na.rm = TRUE) +
    facet_wrap(~Study, scales = "free_x") +
    scale_fill_manual(values = group_colours[c("MG", "MG other", "MG0", "MG00")]) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = expansion(mult = c(0, 0.20))) +
    labs(title = "Conditioning audit: reanalysis DoubletFinder is not independent MG00 evidence", subtitle = "MG00 was assigned after this exact DoubletFinder removal, hence its zero rate is structural. MG other = MG after excluding MG0 and MG00. This panel explains, rather than validates, post-selection rates.", x = NULL, y = "Reanalysis DoubletFinder positive rate") +
    theme_classic(base_size = 8.5) + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 7))
  conditioning_stem <- file.path(figure_dir, "current_reanalysis_DoubletFinder_conditioning_audit")
  save_plot_both(p_conditioning, conditioning_stem, 16, 8); print(p_conditioning)
  figure_rows[[length(figure_rows) + 1L]] <- data.frame(Study = "six-study", figure = "current_reanalysis_DoubletFinder_conditioning_audit", stem = conditioning_stem, stringsAsFactors = FALSE)

  figures <- do.call(rbind, figure_rows)
  figures$png_exists <- file.exists(paste0(figures$stem, ".png")); figures$pdf_exists <- file.exists(paste0(figures$stem, ".pdf"))
  if (any(!figures$png_exists | !figures$pdf_exists)) stop("At least one revised figure is missing PNG or PDF.")
  write.csv(figures, file.path(table_dir, "revised_figure_manifest.csv"), row.names = FALSE)
  write.csv(data.frame(
    key = c("script", "terminology", "comparison_groups", "UMAP_rule", "primary_artificial_reference", "distance_definition", "independent_score_definition", "raw_count_sum_definition", "depth_matched_definition", "MG0_16pct_definition", "invalid_prior_comparison"),
    value = c(
      "26.OL_MG_revised_UMAP_and_comparison_24server.r", "MG_parent replaces canonical_MG_parent in all six revised .qs objects",
      "OL, MG, MG other (confirmed MG excluding MG0 and MG00), MG0, MG00, and artificial OL–MG in raw-count-sum and parent-depth-matched forms. MG0/MG00 retain shared cells in their own label-specific summaries; their intersection is audit-only, not a primary comparison group.",
      "fixed real-cell reference UMAP; one highlighted label per UMAP; identical coordinate limits within study",
      "same-study/same-batch OL–MG counts summed then multinomial-downsampled to mean parent UMI depth",
      "nearest same-batch parent-depth-matched artificial OL–MG in first 20 fixed reference PCs; lower is closer",
      "for each marker g and cell j: log1p(10000 * count_gj / total_UMI_j), then mean over the independent OL or MG marker panel; dual score = min(OL score, MG score)",
      "integer counts(OL parent) + integer counts(MG parent), total UMI equals parent1 + parent2",
      "multinomial downsample of raw artificial profile to round(mean(parent1 UMI, parent2 UMI)); primary comparison",
      "original-pipeline DoubletFinder across all MG0 cells: 514/3174 = 16.2%",
      "post-selection MG0-only reanalysis DoubletFinder rates are not comparable to 16.2% and are retained only as a conditioning audit"
    ), stringsAsFactors = FALSE
  ), file.path(runtime_dir, "revision_manifest.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "sessionInfo.txt"))
  message("\nRevision complete: ", result_root)
})
