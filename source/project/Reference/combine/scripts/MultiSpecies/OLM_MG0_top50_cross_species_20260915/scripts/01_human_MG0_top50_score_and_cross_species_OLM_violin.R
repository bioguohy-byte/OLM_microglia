#!/usr/bin/env Rscript

# Cross-species OLM analysis using the integrated Human/Monkey/Mouse object.
#
# Biological definitions, retained exactly from the source combine_3.r script:
#   Human OLM       = Species == "Human"  & CellType1 == "MG0"
#   Monkey OLM-like = Species == "Monkey" & integrated_snn_res.0.1 == "1"
#   Mouse OLM-like  = Species == "Mouse"  & integrated_snn_res.0.1 == "1"
#
# The Human MG0 top-50 signature is calculated in the RNA assay from
# Human MG0 versus Human MG1--MG5. It replaces the earlier externally supplied
# reference gene set. Every cell is scored using AddModuleScore(ctrl = 10), as
# in the supplied spatial_module.r reference. UMAP coordinates are read only.

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/OLM_MG0_top50_cross_species_20260915")
input_qs <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/02.cross_species_integrated_noref_Human_Mouse_Monkey.qs")
result_dir <- file.path(project_dir, "results", "01_human_MG0_top50_score_cross_species_OLM")
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

# All temporary/cache locations remain in the project rather than the home directory.
Sys.setenv(
  TMPDIR = tmp_dir,
  TMP = tmp_dir,
  TEMP = tmp_dir,
  XDG_CACHE_HOME = file.path(project_dir, "cache"),
  R_USER_CACHE_DIR = file.path(project_dir, "cache", "R")
)
dir.create(Sys.getenv("XDG_CACHE_HOME"), recursive = TRUE, showWarnings = FALSE)
dir.create(Sys.getenv("R_USER_CACHE_DIR"), recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library"), file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(SCP)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260915)

if (!file.exists(input_qs)) stop("Input object does not exist: ", input_qs)

obj <- qread(input_qs)
required_metadata <- c("Species", "CellType1", "integrated_snn_res.0.1")
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata) > 0L) {
  stop("Input object lacks required metadata: ", paste(missing_metadata, collapse = ", "))
}
if (!"RNA" %in% Assays(obj)) stop("Input object lacks an RNA assay.")
if (!"umap" %in% Reductions(obj)) stop("Input object lacks the stored UMAP reduction.")
if (!all(c("Human", "Monkey", "Mouse") %in% unique(as.character(obj$Species)))) {
  stop("All Human, Monkey and Mouse cells must be present.")
}
DefaultAssay(obj) <- "RNA"

# 1. Human MG0 FindMarkers: only Human cells enter the reference contrast.
human_cells <- colnames(obj)[as.character(obj$Species) == "Human"]
human_obj <- subset(obj, cells = human_cells)
DefaultAssay(human_obj) <- "RNA"
Idents(human_obj) <- "CellType1"

human_levels <- levels(Idents(human_obj))
required_human_types <- c("MG0", "MG1", "MG2", "MG3", "MG4", "MG5")
if (!all(required_human_types %in% human_levels)) {
  stop("Human CellType1 is missing one or more MG0--MG5 labels.")
}

human_mg0_markers <- FindMarkers(
  object = human_obj,
  ident.1 = "MG0",
  ident.2 = c("MG1", "MG2", "MG3", "MG4", "MG5"),
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0
)

human_mg0_markers$gene <- rownames(human_mg0_markers)
human_mg0_markers <- human_mg0_markers[
  order(human_mg0_markers$p_val_adj, -human_mg0_markers$avg_log2FC, human_mg0_markers$gene),
  c("gene", setdiff(colnames(human_mg0_markers), "gene")),
  drop = FALSE
]
if (nrow(human_mg0_markers) < 50L) {
  stop("Fewer than 50 positive Human MG0 markers were returned.")
}

human_mg0_top50 <- human_mg0_markers[seq_len(50L), , drop = FALSE]
top50_genes <- human_mg0_top50$gene
if (length(unique(top50_genes)) != 50L) stop("The top-50 signature contains duplicate genes.")
if (length(setdiff(top50_genes, rownames(obj[["RNA"]]))) > 0L) {
  stop("At least one Human MG0 top-50 gene is absent from the integrated RNA assay.")
}

write.table(
  human_mg0_markers,
  file = file.path(table_dir, "Human_MG0_vs_Human_MG1toMG5_FindMarkers_all_positive.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.csv(
  human_mg0_markers,
  file = file.path(table_dir, "Human_MG0_vs_Human_MG1toMG5_FindMarkers_all_positive.csv"),
  row.names = FALSE
)
write.table(
  human_mg0_top50,
  file = file.path(table_dir, "Human_MG0_FindMarkers_top50_reference_signature.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.csv(
  human_mg0_top50,
  file = file.path(table_dir, "Human_MG0_FindMarkers_top50_reference_signature.csv"),
  row.names = FALSE
)

# 2. Score every cell from the identical Human MG0 top-50 signature.
# This is intentionally the same AddModuleScore control-gene design as spatial_module.r.
obj <- AddModuleScore(
  object = obj,
  assay = "RNA",
  features = list(top50_genes),
  ctrl = 10,
  name = "Human_MG0_top50_score",
  seed = 20260915
)
score_column_from_AddModuleScore <- "Human_MG0_top50_score1"
if (!score_column_from_AddModuleScore %in% colnames(obj[[]])) stop("AddModuleScore did not produce ", score_column_from_AddModuleScore)
# A descriptive alias prevents the automatic trailing "1" from appearing in the final figures.
score_column <- "Human_MG0_top50_module_score"
obj[[score_column]] <- obj[[score_column_from_AddModuleScore]][, 1]

# 3. Apply the author-defined cross-species OLM/MG0-like assignments.
obj$OLM_cross_species_definition <- "Other MG"
obj$OLM_cross_species_definition[
  as.character(obj$Species) == "Human" & as.character(obj$CellType1) == "MG0"
] <- "OLM"
obj$OLM_cross_species_definition[
  as.character(obj$Species) == "Monkey" & as.character(obj$integrated_snn_res.0.1) == "1"
] <- "OLM"
obj$OLM_cross_species_definition[
  as.character(obj$Species) == "Mouse" & as.character(obj$integrated_snn_res.0.1) == "1"
] <- "OLM"
obj$OLM_cross_species_definition <- factor(obj$OLM_cross_species_definition, levels = c("OLM", "Other MG"))

obj$OLM_definition_rule <- "Other MG"
obj$OLM_definition_rule[
  as.character(obj$Species) == "Human" & as.character(obj$CellType1) == "MG0"
] <- "Human: CellType1 == MG0"
obj$OLM_definition_rule[
  as.character(obj$Species) == "Monkey" & as.character(obj$integrated_snn_res.0.1) == "1"
] <- "Monkey: integrated_snn_res.0.1 == 1"
obj$OLM_definition_rule[
  as.character(obj$Species) == "Mouse" & as.character(obj$integrated_snn_res.0.1) == "1"
] <- "Mouse: integrated_snn_res.0.1 == 1"

# 4. Stored integrated UMAP.  This is the exact positive-expression drawing
# scheme used for the prior LPC/CPZ panels: all cells are grey first, then only
# score-positive cells are drawn in ascending score order from yellow to dark red.
score_umap_coordinates <- Embeddings(obj, reduction = "umap")
score_plot_frame <- data.frame(
  UMAP_1 = score_umap_coordinates[, 1],
  UMAP_2 = score_umap_coordinates[, 2],
  module_score = as.numeric(obj[[score_column]][, 1]),
  stringsAsFactors = FALSE
)
score_positive_frame <- score_plot_frame[score_plot_frame$module_score > 0, , drop = FALSE]
score_positive_frame <- score_positive_frame[order(score_positive_frame$module_score), , drop = FALSE]
score_upper_limit <- if (nrow(score_positive_frame) == 0L) {
  1
} else {
  as.numeric(quantile(score_positive_frame$module_score, 0.99, names = FALSE))
}
score_umap_all <- ggplot(score_plot_frame, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    color = "#B8B8B8",
    size = 0.20,
    alpha = 0.60,
    shape = 16,
    stroke = 0
  ) +
  geom_point(
    data = score_positive_frame,
    aes(color = module_score),
    size = 0.20 * 1.55,
    alpha = 1,
    shape = 16,
    stroke = 0
  ) +
  scale_color_gradient(
    low = "#FFD166",
    high = "#B2182B",
    limits = c(0, score_upper_limit),
    oob = scales::squish
  ) +
  labs(
    title = "Human MG0 top-50 module score",
    subtitle = "Human MG0 signature scored across all Human, Monkey and Mouse cells",
    color = "Module\nscore",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  coord_equal() +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )
ggsave(file.path(figure_dir, "01_Human_MG0_top50_module_score_all_species_UMAP.png"), score_umap_all, width = 7.2, height = 6.4, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "01_Human_MG0_top50_module_score_all_species_UMAP.pdf"), score_umap_all, width = 7.2, height = 6.4, device = grDevices::cairo_pdf, bg = "white")

# 5. Per-species OLM versus Other MG module-score violins.
# The FeatureStatPlot settings and red/blue palette duplicate the prior M0-vs-Other reference plots.
olm_other_colours <- c("OLM" = "#E64B35", "Other MG" = "#4DBBD5")

human_score_object <- subset(obj, cells = colnames(obj)[as.character(obj$Species) == "Human"])
human_score_object$OLM_vs_Other_MG <- factor(as.character(human_score_object$OLM_cross_species_definition), levels = c("OLM", "Other MG"))
human_violin <- SCP::FeatureStatPlot(
  srt = human_score_object,
  assay = "RNA",
  slot = "data",
  group.by = "OLM_vs_Other_MG",
  stat.by = score_column,
  add_box = FALSE,
  comparisons = list(c("Other MG", "OLM")),
  bg.by = "OLM_vs_Other_MG",
  palcolor = olm_other_colours,
  bg_palcolor = olm_other_colours,
  legend.position = "none",
  xlab = "",
  theme_args = list(base_size = 10),
  ncol = 1
)
human_violin <- human_violin & labs(title = "Human") & theme(
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 11),
  plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
  plot.margin = ggplot2::margin(t = 24, r = 4, b = 4, l = 4, unit = "pt")
)
ggsave(file.path(figure_dir, "02_Human_OLM_vs_Other_MG_top50_score_violin.png"), human_violin, width = 4.6, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "02_Human_OLM_vs_Other_MG_top50_score_violin.pdf"), human_violin, width = 4.6, height = 6.2, device = grDevices::cairo_pdf, bg = "white")

monkey_score_object <- subset(obj, cells = colnames(obj)[as.character(obj$Species) == "Monkey"])
monkey_score_object$OLM_vs_Other_MG <- factor(as.character(monkey_score_object$OLM_cross_species_definition), levels = c("OLM", "Other MG"))
monkey_violin <- SCP::FeatureStatPlot(
  srt = monkey_score_object,
  assay = "RNA",
  slot = "data",
  group.by = "OLM_vs_Other_MG",
  stat.by = score_column,
  add_box = FALSE,
  comparisons = list(c("Other MG", "OLM")),
  bg.by = "OLM_vs_Other_MG",
  palcolor = olm_other_colours,
  bg_palcolor = olm_other_colours,
  legend.position = "none",
  xlab = "",
  theme_args = list(base_size = 10),
  ncol = 1
)
monkey_violin <- monkey_violin & labs(title = "Monkey") & theme(
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 11),
  plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
  plot.margin = ggplot2::margin(t = 24, r = 4, b = 4, l = 4, unit = "pt")
)
ggsave(file.path(figure_dir, "03_Monkey_OLM_vs_Other_MG_top50_score_violin.png"), monkey_violin, width = 4.6, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "03_Monkey_OLM_vs_Other_MG_top50_score_violin.pdf"), monkey_violin, width = 4.6, height = 6.2, device = grDevices::cairo_pdf, bg = "white")

mouse_score_object <- subset(obj, cells = colnames(obj)[as.character(obj$Species) == "Mouse"])
mouse_score_object$OLM_vs_Other_MG <- factor(as.character(mouse_score_object$OLM_cross_species_definition), levels = c("OLM", "Other MG"))
mouse_violin <- SCP::FeatureStatPlot(
  srt = mouse_score_object,
  assay = "RNA",
  slot = "data",
  group.by = "OLM_vs_Other_MG",
  stat.by = score_column,
  add_box = FALSE,
  comparisons = list(c("Other MG", "OLM")),
  bg.by = "OLM_vs_Other_MG",
  palcolor = olm_other_colours,
  bg_palcolor = olm_other_colours,
  legend.position = "none",
  xlab = "",
  theme_args = list(base_size = 10),
  ncol = 1
)
mouse_violin <- mouse_violin & labs(title = "Mouse") & theme(
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 11),
  plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
  plot.margin = ggplot2::margin(t = 24, r = 4, b = 4, l = 4, unit = "pt")
)
ggsave(file.path(figure_dir, "04_Mouse_OLM_vs_Other_MG_top50_score_violin.png"), mouse_violin, width = 4.6, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "04_Mouse_OLM_vs_Other_MG_top50_score_violin.pdf"), mouse_violin, width = 4.6, height = 6.2, device = grDevices::cairo_pdf, bg = "white")

# 6. Per-species statistics and auditable cell/group assignments.
human_score_values <- as.numeric(human_score_object[[score_column]][, 1])
human_groups <- as.character(human_score_object$OLM_vs_Other_MG)
monkey_score_values <- as.numeric(monkey_score_object[[score_column]][, 1])
monkey_groups <- as.character(monkey_score_object$OLM_vs_Other_MG)
mouse_score_values <- as.numeric(mouse_score_object[[score_column]][, 1])
mouse_groups <- as.character(mouse_score_object$OLM_vs_Other_MG)

human_test <- wilcox.test(human_score_values[human_groups == "OLM"], human_score_values[human_groups == "Other MG"], alternative = "two.sided", exact = FALSE)
monkey_test <- wilcox.test(monkey_score_values[monkey_groups == "OLM"], monkey_score_values[monkey_groups == "Other MG"], alternative = "two.sided", exact = FALSE)
mouse_test <- wilcox.test(mouse_score_values[mouse_groups == "OLM"], mouse_score_values[mouse_groups == "Other MG"], alternative = "two.sided", exact = FALSE)

species_statistics <- data.frame(
  Species = c("Human", "Monkey", "Mouse"),
  OLM_definition = c("CellType1 == MG0", "integrated_snn_res.0.1 == 1", "integrated_snn_res.0.1 == 1"),
  OLM_cells = c(sum(human_groups == "OLM"), sum(monkey_groups == "OLM"), sum(mouse_groups == "OLM")),
  Other_MG_cells = c(sum(human_groups == "Other MG"), sum(monkey_groups == "Other MG"), sum(mouse_groups == "Other MG")),
  OLM_mean_module_score = c(mean(human_score_values[human_groups == "OLM"]), mean(monkey_score_values[monkey_groups == "OLM"]), mean(mouse_score_values[mouse_groups == "OLM"])),
  OLM_median_module_score = c(median(human_score_values[human_groups == "OLM"]), median(monkey_score_values[monkey_groups == "OLM"]), median(mouse_score_values[mouse_groups == "OLM"])),
  Other_MG_mean_module_score = c(mean(human_score_values[human_groups == "Other MG"]), mean(monkey_score_values[monkey_groups == "Other MG"]), mean(mouse_score_values[mouse_groups == "Other MG"])),
  Other_MG_median_module_score = c(median(human_score_values[human_groups == "Other MG"]), median(monkey_score_values[monkey_groups == "Other MG"]), median(mouse_score_values[mouse_groups == "Other MG"])),
  Wilcoxon_W = c(unname(human_test$statistic), unname(monkey_test$statistic), unname(mouse_test$statistic)),
  p_value_unadjusted = c(human_test$p.value, monkey_test$p.value, mouse_test$p.value),
  stringsAsFactors = FALSE
)
species_statistics$p_value_BH_three_species <- p.adjust(species_statistics$p_value_unadjusted, method = "BH")
species_statistics$plot_significance <- ifelse(species_statistics$p_value_unadjusted <= 1e-4, "****", ifelse(species_statistics$p_value_unadjusted <= 1e-3, "***", ifelse(species_statistics$p_value_unadjusted <= 1e-2, "**", ifelse(species_statistics$p_value_unadjusted <= 0.05, "*", "ns"))))
write.table(species_statistics, file.path(table_dir, "Human_Monkey_Mouse_OLM_vs_Other_MG_top50_score_Wilcoxon.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

umap_coordinates <- Embeddings(obj, reduction = "umap")
cell_audit <- data.frame(
  cell_barcode = colnames(obj),
  Species = as.character(obj$Species),
  CellType1 = as.character(obj$CellType1),
  integrated_snn_res_0.1 = as.character(obj$integrated_snn_res.0.1),
  OLM_definition_rule = as.character(obj$OLM_definition_rule),
  OLM_vs_Other_MG = as.character(obj$OLM_cross_species_definition),
  Human_MG0_top50_score = as.numeric(obj[[score_column]][, 1]),
  UMAP_1 = umap_coordinates[, 1],
  UMAP_2 = umap_coordinates[, 2],
  stringsAsFactors = FALSE
)
write.table(cell_audit, file.path(table_dir, "per_cell_Human_MG0_top50_score_OLM_definition_UMAP.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

qsave(obj, file.path(object_dir, "02_cross_species_integrated_with_Human_MG0_top50_score.qs"), preset = "high")

parameter_record <- data.frame(
  item = c(
    "input_object", "Human_marker_contrast", "FindMarkers_assay_layer", "FindMarkers_test", "FindMarkers_only_positive", "FindMarkers_min_pct", "FindMarkers_logfc_threshold", "signature_size", "module_scoring", "module_score_assay", "module_score_ctrl", "module_score_seed", "Human_OLM_definition", "Monkey_OLM_definition", "Mouse_OLM_definition", "UMAP_source", "UMAP_colours", "violin_plotter", "violin_colours", "violin_test", "plot_star_rule"
  ),
  value = c(
    input_qs, "MG0 versus MG1, MG2, MG3, MG4, MG5 among Human cells only", "RNA/data", "two-sided Wilcoxon rank-sum", "TRUE", "0.10", "0", "top 50 positive markers ordered by adjusted p value, then avg_log2FC", "Seurat::AddModuleScore", "RNA", "10", "20260915", "Species == Human and CellType1 == MG0", "Species == Monkey and integrated_snn_res.0.1 == 1", "Species == Mouse and integrated_snn_res.0.1 == 1", "stored umap reduction; not recomputed", "LPC/CPZ convention: all cells #B8B8B8 (size 0.20, alpha 0.60); score >0 cells yellow #FFD166 to dark red #B2182B (size 0.31, alpha 1); positive-score q99 upper limit", "SCP::FeatureStatPlot (same reference-style parameters)", "OLM #E64B35; Other MG #4DBBD5", "two-sided Wilcoxon rank-sum", "**** <=1e-4; *** <=1e-3; ** <=1e-2; * <=0.05; ns otherwise; plotted p is unadjusted"
  ),
  stringsAsFactors = FALSE
)
write.table(parameter_record, file.path(result_dir, "01_parameter_record.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "01_human_MG0_top50_score_cross_species_sessionInfo.txt"))
writeLines(as.character(Sys.time()), con = file.path(log_dir, "01_human_MG0_top50_score_cross_species_completed_at.txt"))
message("Completed: ", result_dir)
