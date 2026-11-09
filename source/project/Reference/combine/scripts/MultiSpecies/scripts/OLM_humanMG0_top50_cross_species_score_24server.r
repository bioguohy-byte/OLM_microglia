# Human OLM (MG0) top50 signature scoring across the original
# Human-Monkey-Mouse integrated microglia object.
# Reference adapted from spatial_module.r:
# AddModuleScore(features = list(top50), ctrl = 10).
# SCP is not installed on this server, so the two-group display below is a
# self-contained ggplot2 violin + Wilcoxon annotation.

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
})

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies")
input_qs <- file.path(
  project_dir,
  "data/processed_data/02.cross_species_integrated_noref_Human_Mouse_Monkey.qs"
)
fig_dir <- file.path(project_dir, "results/figure/OLM_humanMG0_top50_score")
table_dir <- file.path(project_dir, "results/table/OLM_humanMG0_top50_score")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_qs)) stop("Input object is absent: ", input_qs)
obj <- qread(input_qs)
if (!all(c("RNA", "integrated") %in% names(obj@assays))) {
  stop("Expected RNA and integrated assays.")
}
if (!all(c("Human", "Monkey", "Mouse") %in% unique(as.character(obj$Species)))) {
  stop("Missing one or more required species.")
}
DefaultAssay(obj) <- "RNA"

# 1. Human OLM top50: Human MG0 versus Human MG1-MG5 in RNA.
human_obj <- subset(obj, subset = Species == "Human")
DefaultAssay(human_obj) <- "RNA"
Idents(human_obj) <- "CellType1"
other_human_mg <- c("MG1", "MG2", "MG3", "MG4", "MG5")
human_mg0_de <- FindMarkers(
  human_obj,
  ident.1 = "MG0",
  ident.2 = other_human_mg,
  only.pos = TRUE,
  logfc.threshold = 0,
  min.pct = 0.10,
  test.use = "wilcox"
)
# Rank adjusted P value first; exact ties are ranked by decreasing effect size.
human_mg0_de <- human_mg0_de[
  order(human_mg0_de$p_val_adj, -human_mg0_de$avg_log2FC),
  , drop = FALSE
]
human_mg0_top50 <- head(human_mg0_de, 50)
human_mg0_top50 <- cbind(gene = rownames(human_mg0_top50), human_mg0_top50)
rownames(human_mg0_top50) <- NULL
if (nrow(human_mg0_top50) != 50) {
  stop("FindMarkers returned fewer than 50 positive features.")
}
score_genes <- human_mg0_top50$gene
missing_score_genes <- setdiff(score_genes, rownames(obj[["RNA"]]))
if (length(missing_score_genes) > 0L) {
  stop("Top50 genes absent from unified RNA matrix: ",
       paste(missing_score_genes, collapse = ", "))
}
utils::write.csv(
  human_mg0_de,
  file.path(table_dir, "human_MG0_vs_MG1to5_FindMarkers_all_positive.csv"),
  row.names = TRUE
)
utils::write.csv(
  human_mg0_top50,
  file.path(table_dir, "human_MG0_FindMarkers_top50.csv"),
  row.names = FALSE
)

# 2. Score ALL 41,733 cells with the same AddModuleScore design as the reference.
score_prefix <- "Human_OLM_MG0_top50_score"
score_col <- paste0(score_prefix, "1")
obj <- AddModuleScore(
  object = obj,
  features = list(score_genes),
  ctrl = 10,
  name = score_prefix,
  assay = "RNA",
  seed = 123
)

# 3. Reproduce the original script's mouse/monkey OLM-like definitions.
obj$new_celltype <- as.character(obj$CellType1)
obj$new_celltype[
  as.character(obj$integrated_snn_res.0.1) == "1" & obj$Species == "Mouse"
] <- "Mouse_MG0_like"
obj$new_celltype[
  as.character(obj$integrated_snn_res.0.1) == "1" & obj$Species == "Monkey"
] <- "Monkey_MG0_like"
obj$new_celltype <- ifelse(
  obj$new_celltype %in% c("MG0", "Mouse_MG0_like", "Monkey_MG0_like"),
  obj$new_celltype,
  "Other_MG"
)

# 4. Whole-object score FeaturePlot: original UMAP coordinates, all cells.
p_score_all <- FeaturePlot(
  object = obj,
  features = score_col,
  reduction = "umap",
  cols = c("grey85", "#D73027"),
  order = TRUE,
  pt.size = 0.25,
  raster = FALSE
) +
  ggtitle("Human OLM (MG0) FindMarkers top50 score: Human, Monkey and Mouse") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_score_all)
ggsave(
  file.path(fig_dir, "OLM_humanMG0_top50_score_feature_UMAP_all_cells.png"),
  p_score_all, width = 8, height = 6.5, dpi = 300
)
ggsave(
  file.path(fig_dir, "OLM_humanMG0_top50_score_feature_UMAP_all_cells.pdf"),
  p_score_all, width = 8, height = 6.5
)

# 5. Per-species OLM/Other-MG score violin plots with Wilcoxon annotation.
make_score_violin <- function(seu, species_name, positive_label,
                              positive_colour, file_stub) {
  keep_cells <- colnames(seu)[as.character(seu$Species) == species_name]
  sub <- subset(seu, cells = keep_cells)
  sub$score_group <- ifelse(
    as.character(sub$new_celltype) == positive_label,
    positive_label,
    "Other MG"
  )
  sub$score_group <- factor(
    sub$score_group,
    levels = c(positive_label, "Other MG")
  )
  dat <- FetchData(sub, vars = c(score_col, "score_group"))
  colnames(dat) <- c("score", "score_group")
  dat <- dat[is.finite(dat$score) & !is.na(dat$score_group), , drop = FALSE]
  if (length(unique(dat$score_group)) != 2L || any(table(dat$score_group) == 0L)) {
    stop("Both groups are required for ", species_name)
  }
  wt <- stats::wilcox.test(score ~ score_group, data = dat, exact = FALSE)
  max_y <- max(dat$score)
  min_y <- min(dat$score)
  pad <- max((max_y - min_y) * 0.12, 0.10)
  bracket_y <- max_y + pad
  p_label <- paste0(
    "Wilcoxon p ",
    format.pval(wt$p.value, digits = 3, eps = 1e-300)
  )
  fill_cols <- c(positive_colour, "#D9D9D9")
  names(fill_cols) <- c(positive_label, "Other MG")
  p <- ggplot(dat, aes(x = score_group, y = score, fill = score_group)) +
    geom_violin(trim = FALSE, colour = "black", linewidth = 0.35) +
    geom_boxplot(
      width = 0.13,
      outlier.shape = NA,
      colour = "black",
      linewidth = 0.35
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      colour = "white",
      size = 1.3
    ) +
    annotate(
      "segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
      linewidth = 0.4
    ) +
    annotate(
      "segment", x = 1, xend = 1,
      y = bracket_y - pad * 0.12, yend = bracket_y,
      linewidth = 0.4
    ) +
    annotate(
      "segment", x = 2, xend = 2,
      y = bracket_y - pad * 0.12, yend = bracket_y,
      linewidth = 0.4
    ) +
    annotate(
      "text", x = 1.5, y = bracket_y + pad * 0.10,
      label = p_label, size = 3.6
    ) +
    scale_fill_manual(values = fill_cols) +
    labs(
      title = paste0(species_name, ": Human OLM (MG0) top50 module score"),
      subtitle = paste0(positive_label, " versus Other MG"),
      x = NULL,
      y = "AddModuleScore (top50; ctrl = 10)"
    ) +
    coord_cartesian(ylim = c(min_y - pad * 0.05, bracket_y + pad * 0.28)) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 15, hjust = 1)
    )
  group_counts <- table(dat$score_group)
  group_median <- tapply(dat$score, dat$score_group, median)
  group_mean <- tapply(dat$score, dat$score_group, mean)
  smry <- data.frame(
    species = species_name,
    group = levels(dat$score_group),
    n = as.integer(group_counts[levels(dat$score_group)]),
    median_score = as.numeric(group_median[levels(dat$score_group)]),
    mean_score = as.numeric(group_mean[levels(dat$score_group)]),
    wilcoxon_p_value = wt$p.value,
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    smry,
    file.path(table_dir, paste0(file_stub, "_summary.csv")),
    row.names = FALSE
  )
  utils::write.csv(
    dat,
    file.path(table_dir, paste0(file_stub, "_cell_scores.csv")),
    row.names = TRUE
  )
  ggsave(
    file.path(fig_dir, paste0(file_stub, ".png")),
    p, width = 5.4, height = 5, dpi = 300
  )
  ggsave(
    file.path(fig_dir, paste0(file_stub, ".pdf")),
    p, width = 5.4, height = 5
  )
  print(p)
  list(plot = p, summary = smry)
}

human_violin <- make_score_violin(
  obj, "Human", "MG0", "#4DAF4A",
  "OLM_humanMG0_top50_score_violin_Human"
)
monkey_violin <- make_score_violin(
  obj, "Monkey", "Monkey_MG0_like", "#B2DF8A",
  "OLM_humanMG0_top50_score_violin_Monkey"
)
mouse_violin <- make_score_violin(
  obj, "Mouse", "Mouse_MG0_like", "#33A02C",
  "OLM_humanMG0_top50_score_violin_Mouse"
)

# 6. Auditable score/UMAP mapping table. The original .qs remains untouched.
umap_xy <- Embeddings(obj, reduction = "umap")
score_metadata <- data.frame(
  cell_barcode = colnames(obj),
  Species = as.character(obj$Species),
  CellType1 = as.character(obj$CellType1),
  integrated_snn_res_0_1 = as.character(obj$integrated_snn_res.0.1),
  original_cross_species_group = as.character(obj$new_celltype),
  OLM_humanMG0_top50_score = as.numeric(obj[[score_col]][, 1]),
  UMAP_1 = umap_xy[, 1],
  UMAP_2 = umap_xy[, 2],
  row.names = NULL,
  check.names = FALSE
)
utils::write.csv(
  score_metadata,
  file.path(table_dir, "OLM_humanMG0_top50_score_all_cells_UMAP_metadata.csv"),
  row.names = FALSE
)
run_manifest <- data.frame(
  item = c(
    "input_object", "assay_for_DE_and_score", "DE_comparison", "DE_test",
    "only_positive", "min_pct", "logfc_threshold", "top_n",
    "top50_present_in_RNA", "score_method", "score_control_genes",
    "total_cells", "human_cells", "monkey_cells", "mouse_cells"
  ),
  value = c(
    input_qs, "RNA", "Human MG0 vs Human MG1-MG5", "Wilcoxon",
    "TRUE", "0.10", "0", "50", length(score_genes),
    "Seurat::AddModuleScore", "10", ncol(obj),
    sum(obj$Species == "Human"), sum(obj$Species == "Monkey"),
    sum(obj$Species == "Mouse")
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  run_manifest,
  file.path(table_dir, "OLM_humanMG0_top50_score_manifest.csv"),
  row.names = FALSE
)

message("Completed. Figures: ", fig_dir)
message("Completed. Tables: ", table_dir)
message("Top50 present in unified RNA matrix: ", length(score_genes), "/50")