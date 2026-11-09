#!/usr/bin/env Rscript
# Figure-ordered code from the reviewed manuscript analysis package.
# Each panel contains its original quoted R code.
# Source this chapter, then use mg_list_panels() or mg_show_panel("PANEL_ID").
# Source paths in panel records identify the original analysis scripts.
.mg_panels <- new.env(parent=emptyenv())
.mg_order <- character()
mg_register_panel <- function(id, code) {
  if (!id %in% .mg_order) .mg_order <<- c(.mg_order,id)
  assign(id, code, envir=.mg_panels)
}
mg_list_panels <- function() .mg_order
mg_panel_code <- function(id) {
  if (!exists(id,envir=.mg_panels,inherits=FALSE)) stop("Unknown panel: ",id)
  get(id,envir=.mg_panels,inherits=FALSE)
}
mg_show_panel <- function(id) {
  blocks <- mg_panel_code(id)
  if (is.null(blocks)) {cat(id, ": No R block is included for this panel.\n");return(invisible(NULL))}
  for (block in blocks) {
    cat("# Source: ",block$source," [actual ",block$ranges,"; requested ",block$requested_ranges,"]\n# Role: ",block$role,"\n",sep="")
    cat(paste(deparse(block$code,width.cutoff=120L),collapse="\n"),"\n\n")
  }
  invisible(blocks)
}


# Supplementary Figure 4

# S4A  Human OLM top-50 module score over all Human/Monkey/Mouse microglia, positive-only yellow/red UMAP
# Kind: computational
mg_register_panel("S4A", list(
  list(source="source/project/Reference/combine/scripts/MultiSpecies/OLM_MG0_top50_cross_species_20260915/scripts/01_human_MG0_top50_score_and_cross_species_OLM_violin.R", ranges="1-227", requested_ranges="1-227", role="complete historical setup and score/plot context; outputs remain historical", code=quote({
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
  }))
))

# S4B  Within each species: OLM top-50 module-score violin, OLM versus Other MG
# Kind: computational
mg_register_panel("S4B", list(
  list(source="source/project/Reference/combine/scripts/MultiSpecies/OLM_MG0_top50_cross_species_20260915/scripts/01_human_MG0_top50_score_and_cross_species_OLM_violin.R", ranges="1-344", requested_ranges="1-344", role="complete historical setup and score/plot context; outputs remain historical", code=quote({
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
  }))
))

# S4C  Monkey_MG0_like top10 GO Biological Process barplot
# Kind: computational
mg_register_panel("S4C", list(
  list(source="source/project/Reference/combine/scripts/MultiSpecies/scripts/Monkey_MG0_like_GO.r", ranges="1-142", requested_ranges="1-142", role="GO enrichment and actual GO barplot", code=quote({
library(qs)
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(stringr)
library(ggplot2)

# MODULE : GO Enrichment and Barplot for a Single Cluster -------------------

## 1 Parameters --------------------------------------------------------------
# -- 数据输入参数 --
input_filepath  <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/monkey_mg_markers.qs")
target_cluster  <- "Monkey_MG0_like"       # 目标提取的 cluster 名称

# -- 过滤与分析阈值 --
p_val_adj_cut   <- 0.05                   # 显著性阈值
log2FC_cut      <- 0.5                    # Log2FC 阈值
db_org          <- "org.Mmu.eg.db"         # 物种数据库: "org.Mm.eg.db" (鼠) / "org.Hs.eg.db" (人) / "org.Mmu.eg.db" (恒河猴)
ontology_type   <- "BP"                   # GO 类型: "BP", "CC", "MF" 或 "ALL"
top_n_go        <- 10                     # 图中展示的 Top GO 通路数量

# -- 绘图与美化参数 --
cluster_color   <- "#B2DF8A"              # 柱状图的填充颜色
wrap_width      <- 60                     # 文本自动换行的字符数限制

# -- 输出参数 --
out_dir         <- "results/figure/cluster_go/"  # 图片保存的文件夹路径
out_width       <- 8                      # 输出图片宽度
out_height      <- 5                      # 输出图片高度


## 2 Process Data ------------------------------------------------------------
# 动态加载物种数据库
require(db_org, character.only = TRUE)

# 1. 读取 Marker 数据框 (修复：使用 qread 真正载入数据)
marker_df <- qread(input_filepath)

# 2. 过滤目标群的 Marker
cluster_markers <- marker_df %>%
  dplyr::filter(cluster == target_cluster) %>%
  dplyr::filter(p_val_adj < p_val_adj_cut) %>% 
  dplyr::filter(avg_log2FC > log2FC_cut) %>%
  dplyr::arrange(desc(avg_log2FC))

# 3. ID 转换 (SYMBOL -> ENTREZID)
gene_conv <- bitr(cluster_markers$gene, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = db_org)

# 4. 运行 GO 富集分析
ego_cluster <- enrichGO(
  gene          = gene_conv$ENTREZID,
  OrgDb         = db_org,
  ont           = ontology_type,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# 5. 提取富集结果并准备作图数据
plot_data <- as.data.frame(ego_cluster)
dir.create(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/table"), recursive = TRUE)
write.xlsx(plot_data,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/table/monkey_go.xlsx"))

if(nrow(plot_data) > 0) {
  plot_data <- plot_data %>%
    dplyr::arrange(p.adjust) %>%
    head(top_n_go) %>%
    dplyr::mutate(
      logP = -log10(p.adjust),
      # 对过长的 Description 自动换行
      Description = stringr::str_wrap(Description, width = wrap_width) 
    )
  
  # 锁定因子顺序，确保显著性最高的排在最上面
  plot_data$Description <- factor(plot_data$Description, levels = rev(plot_data$Description))
} else {
  warning(paste("Warning: No significant GO terms found for", target_cluster))
}


## 3 Generate Plots ----------------------------------------------------------
# 初始化画图对象为空
p_barplot <- NULL

if(nrow(plot_data) > 0) {
  
  # 计算文字的 X 轴起始坐标 (放在坐标轴起点偏右一点点，确保整齐对齐)
  text_x_pos <- max(plot_data$logP) * 0.02 
  
  p_barplot <- ggplot(plot_data, aes(x = logP, y = Description)) +
    geom_col(width = 0.8, alpha = 0.8, fill = cluster_color, color = "black", linewidth = 0.3) +
    # 文字标签：固定在图表左侧对齐
    geom_text(aes(x = text_x_pos, label = Description), 
              size = 4.5, hjust = 0, color = "black") + 
    theme_minimal() +
    labs(
      title = paste0("Cluster ", target_cluster, " Top ", top_n_go, " GO Terms (", ontology_type, ")"), 
      x = "-log10(p.adjust)", 
      y = NULL
    ) +
    theme(
      axis.text.y  = element_blank(), # 隐藏原生 Y 轴文字
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(), # 隐藏网格线，画面更干净
      axis.line    = element_line(color = "black"),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    # 向右延展 X 轴空间，防止特别长的通路名字越界被裁切
    scale_x_continuous(expand = expansion(mult = c(0, 0.4))) 
}


## 4 Save Results ------------------------------------------------------------
if(!is.null(p_barplot)) {
  # 如果文件夹不存在，自动创建
  if(!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # 替换 cluster 名字里的非法字符（如果有），防止保存报错
  safe_name <- gsub("[^a-zA-Z0-9_]", "_", target_cluster) 
  file_base <- file.path(out_dir, paste0("Cluster_", safe_name, "_GO_", ontology_type))
  
  # 保存 PDF
  ggsave(filename = paste0(file_base, ".pdf"), 
         plot     = p_barplot, 
         width    = out_width, 
         height   = out_height)
  
  # 保存 PNG
  ggsave(filename = paste0(file_base, ".png"), 
         plot     = p_barplot, 
         width    = out_width, 
         height   = out_height, 
         dpi      = 300)
  
  print(paste("✅ GO barplot successfully generated and saved for:", target_cluster))
}
  })),
  list(source="source/project/Reference/combine/scripts/MultiSpecies/scripts/finallmarkers.r", ranges="15-23", requested_ranges="15-23", role="upstream marker-generation provenance; not a subsequent plot step", code=quote({
# Monkey findallmarkers----------------
monkey_wm_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/monkey_hip_mg_v2.qs"))
DefaultAssay(monkey_wm_mg) <- "RNA"
Idents(monkey_wm_mg) <- "new_celltype"
monkey_mg_markers <- FindAllMarkers(monkey_wm_mg, 
                                   only.pos = TRUE, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.25)
qsave(monkey_mg_markers,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/monkey_mg_markers.qs"))
  }))
))

# S4D  Mouse_MG0_like top10 GO Biological Process barplot
# Kind: computational
mg_register_panel("S4D", list(
  list(source="source/project/Reference/combine/scripts/MultiSpecies/scripts/Mouse_M0_like_GO.r", ranges="1-142", requested_ranges="1-142", role="GO enrichment and actual GO barplot", code=quote({
library(qs)
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(stringr)
library(ggplot2)

# MODULE : GO Enrichment and Barplot for a Single Cluster -------------------

## 1 Parameters --------------------------------------------------------------
# -- 数据输入参数 --
input_filepath  <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/mouse_mg_markers.qs")
target_cluster  <- "Mouse_MG0_like"       # 目标提取的 cluster 名称

# -- 过滤与分析阈值 --
p_val_adj_cut   <- 0.05                   # 显著性阈值
log2FC_cut      <- 0.5                    # Log2FC 阈值
db_org          <- "org.Mm.eg.db"         # 物种数据库: "org.Mm.eg.db" (鼠) / "org.Hs.eg.db" (人)
ontology_type   <- "BP"                   # GO 类型: "BP", "CC", "MF" 或 "ALL"
top_n_go        <- 10                     # 图中展示的 Top GO 通路数量

# -- 绘图与美化参数 --
cluster_color   <- "#33A02C"              # 柱状图的填充颜色
wrap_width      <- 60                     # 文本自动换行的字符数限制

# -- 输出参数 --
out_dir         <- "results/figure/cluster_go/"  # 图片保存的文件夹路径
out_width       <- 8                      # 输出图片宽度
out_height      <- 5                      # 输出图片高度


## 2 Process Data ------------------------------------------------------------
# 动态加载物种数据库
require(db_org, character.only = TRUE)

# 1. 读取 Marker 数据框 (修复：使用 qread 真正载入数据)
marker_df <- qread(input_filepath)

# 2. 过滤目标群的 Marker
cluster_markers <- marker_df %>%
  dplyr::filter(cluster == target_cluster) %>%
  dplyr::filter(p_val_adj < p_val_adj_cut) %>% 
  dplyr::filter(avg_log2FC > log2FC_cut) %>%
  dplyr::arrange(desc(avg_log2FC))

# 3. ID 转换 (SYMBOL -> ENTREZID)
gene_conv <- bitr(cluster_markers$gene, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = db_org)

# 4. 运行 GO 富集分析
ego_cluster <- enrichGO(
  gene          = gene_conv$ENTREZID,
  OrgDb         = db_org,
  ont           = ontology_type,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# 5. 提取富集结果并准备作图数据
plot_data <- as.data.frame(ego_cluster)
dir.create(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/table"), recursive = TRUE)
write.xlsx(plot_data,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/table/mouse_go.xlsx"))

if(nrow(plot_data) > 0) {
  plot_data <- plot_data %>%
    dplyr::arrange(p.adjust) %>%
    head(top_n_go) %>%
    dplyr::mutate(
      logP = -log10(p.adjust),
      # 对过长的 Description 自动换行
      Description = stringr::str_wrap(Description, width = wrap_width) 
    )
  
  # 锁定因子顺序，确保显著性最高的排在最上面
  plot_data$Description <- factor(plot_data$Description, levels = rev(plot_data$Description))
} else {
  warning(paste("Warning: No significant GO terms found for", target_cluster))
}


## 3 Generate Plots ----------------------------------------------------------
# 初始化画图对象为空
p_barplot <- NULL

if(nrow(plot_data) > 0) {
  
  # 计算文字的 X 轴起始坐标 (放在坐标轴起点偏右一点点，确保整齐对齐)
  text_x_pos <- max(plot_data$logP) * 0.02 
  
  p_barplot <- ggplot(plot_data, aes(x = logP, y = Description)) +
    geom_col(width = 0.8, alpha = 0.8, fill = cluster_color, color = "black", linewidth = 0.3) +
    # 文字标签：固定在图表左侧对齐
    geom_text(aes(x = text_x_pos, label = Description), 
              size = 4.5, hjust = 0, color = "black") + 
    theme_minimal() +
    labs(
      title = paste0("Cluster ", target_cluster, " Top ", top_n_go, " GO Terms (", ontology_type, ")"), 
      x = "-log10(p.adjust)", 
      y = NULL
    ) +
    theme(
      axis.text.y  = element_blank(), # 隐藏原生 Y 轴文字
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(), # 隐藏网格线，画面更干净
      axis.line    = element_line(color = "black"),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    # 向右延展 X 轴空间，防止特别长的通路名字越界被裁切
    scale_x_continuous(expand = expansion(mult = c(0, 0.4))) 
}


## 4 Save Results ------------------------------------------------------------
if(!is.null(p_barplot)) {
  # 如果文件夹不存在，自动创建
  if(!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # 替换 cluster 名字里的非法字符（如果有），防止保存报错
  safe_name <- gsub("[^a-zA-Z0-9_]", "_", target_cluster) 
  file_base <- file.path(out_dir, paste0("Cluster_", safe_name, "_GO_", ontology_type))
  
  # 保存 PDF
  ggsave(filename = paste0(file_base, ".pdf"), 
         plot     = p_barplot, 
         width    = out_width, 
         height   = out_height)
  
  # 保存 PNG
  ggsave(filename = paste0(file_base, ".png"), 
         plot     = p_barplot, 
         width    = out_width, 
         height   = out_height, 
         dpi      = 300)
  
  print(paste("✅ GO barplot successfully generated and saved for:", target_cluster))
}
  })),
  list(source="source/project/Reference/combine/scripts/MultiSpecies/scripts/finallmarkers.r", ranges="4-12", requested_ranges="4-12", role="upstream marker-generation provenance; not a subsequent plot step", code=quote({
# Mouse findallmarkers----------------
mouse_wm_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/mouse_wm_mg_v2.qs"))
DefaultAssay(mouse_wm_mg) <- "RNA"
Idents(mouse_wm_mg) <- "new_celltype"
mouse_mg_markers <- FindAllMarkers(mouse_wm_mg, 
                                   only.pos = TRUE, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.25)
qsave(mouse_mg_markers,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/mouse_mg_markers.qs"))
  }))
))

# S4E  LPC GSE148676 stored UMAP with merged cluster0-cluster10 labels
# Kind: computational
mg_register_panel("S4E", list(
  list(source="source/external/LPC/scripts/34_gse148676_reference_spatial_module.R", ranges="1-130;158-164", requested_ranges="1-130;158-164", role="complete LPC input/mapping/score context and actual UMAP export", code=quote({
# MODULE : Cell Type Annotation and Module Scoring -------------------
#
# Direct adaptation of spatial_module.r for merged GSE148676 cells.
# The original package calls and plotting structure are retained. Only the
# input/output paths and the GSE148676 cluster mapping differ. The saveRDS()
# directly after AddModuleScore makes the original reload-for-plot step work
# when this script is sourced from a fresh R session.

## 1 Parameters --------------------------------------------------------------
# Load dependencies
library(Seurat)
library(qs)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SCP)

# File paths
project <- Sys.getenv("MG_LPC_ROOT")
file_in_merfish <- file.path(
  project,
  "objects/simple_microglia/mouse_GSE148676_merged_groups_seurat_v5.rds"
)
file_in_markers <- file.path(
  project,
  "raw_data/reference/mouse_mg_markers.qs"
)
file_out_merfish <- file.path(
  project,
  "objects/simple_microglia/mouse_GSE148676_merged_groups_reference_module_score.rds"
)
file_out_umap_pdf <- file.path(
  project,
  "figures/simple_microglia/mouse_GSE148676_reference_module_score_UMAP.pdf"
)
file_out_umap_png <- sub("\\.pdf$", ".png", file_out_umap_pdf)
file_out_violin_pdf <- file.path(
  project,
  "figures/simple_microglia/mouse_GSE148676_reference_module_score_violin.pdf"
)
file_out_violin_png <- sub("\\.pdf$", ".png", file_out_violin_pdf)
dir.create(dirname(file_out_merfish), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(file_out_umap_pdf), recursive = TRUE, showWarnings = FALSE)

# Cell-type mapping dictionary (old merged-cluster label -> plot label)
col_old_celltype <- "merged_cluster"
col_new_celltype <- "plot_celltype"
celltype_dict <- c(
  "0" = "cluster0",
  "1" = "cluster1",
  "2" = "cluster2",
  "3" = "cluster3",
  "4" = "cluster4",
  "5" = "cluster5",
  "6" = "cluster6",
  "7" = "cluster7",
  "8" = "cluster8",
  "9" = "cluster9",
  "10" = "cluster10"
)

# Marker extraction and module-score parameters
marker_cluster_col <- "cluster"
target_cluster <- "cluster2"
marker_gene_col <- "gene"

score_module_name <- "myelin"
score_feature_name <- paste0(score_module_name, "1")
score_ctrl_size <- 10

# Plotting parameter retained from the reference script
plot_pt_size <- 0


## 2 Process Data ------------------------------------------------------------
# 2.1 Read data
merfish <- readRDS(file_in_merfish)
mouse_mg_markers <- qread(file_in_markers)

# 2.2 Replace cell-type names
merfish@meta.data[[col_new_celltype]] <- as.character(
  merfish@meta.data[[col_old_celltype]]
)
for (old_name in names(celltype_dict)) {
  new_name <- celltype_dict[[old_name]]
  merfish@meta.data[[col_new_celltype]][
    merfish@meta.data[[col_new_celltype]] %in% old_name
  ] <- new_name
}

# 2.3 Extract target-cluster genes and intersect with GSE148676 genes
filtered_data <- mouse_mg_markers %>%
  filter(.data[[marker_cluster_col]] == target_cluster)
int_gene <- intersect(filtered_data[[marker_gene_col]], rownames(merfish))

# 2.4 Calculate module score
merfish <- AddModuleScore(
  object = merfish,
  features = list(int_gene),
  ctrl = score_ctrl_size,
  name = score_module_name
)

# The reference SCP function accesses the Seurat v4 RNA@data slot. This
# conversion preserves the existing normalized RNA layer for plotting only.
merfish[["RNA"]] <- as(merfish[["RNA"]], "Assay")

# The reference script re-opens a previously saved scored object for plotting.
# Save here so the same operation is valid in a new session.
saveRDS(merfish, file = file_out_merfish)


## 3 Generate Plots ----------------------------------------------------------
merfish <- readRDS(file_out_merfish)

st_mg_colors_npg <- c(
  "cluster1" = "#E64B35",
  "cluster0" = "#4DBBD5",
  "cluster2" = "#00A087",
  "cluster3" = "#3C5488",
  "cluster4" = "#F39B7F",
  "cluster5" = "#7E6148",
  "cluster6" = "#8E6BBE",
  "cluster7" = "#7F7F7F",
  "cluster8" = "#66C2A5",
  "cluster9" = "#E6AB02",
  "cluster10" = "#A65628"
)

(p1 <- DimPlot(merfish, group.by = col_new_celltype, cols = st_mg_colors_npg))

## 4 Save Results ------------------------------------------------------------
# Save updated Seurat object (with the new cell-type and module score)
saveRDS(merfish, file = file_out_merfish)

# Save standalone UMAP and violin panels
ggsave(filename = file_out_umap_pdf, plot = p1, width = 5, height = 5)
ggsave(filename = file_out_umap_png, plot = p1, width = 5, height = 5, dpi = 300)
  })),
  list(source="source/external/LPC/scripts/31_merge_gse148676_clusters_export.R", ranges="1-68", requested_ranges="1-68", role="upstream cluster merge provenance; not a subsequent plot step", code=quote({
#!/usr/bin/env Rscript

# GSE148676 microglia: merge selected Seurat clusters without changing the
# original seurat_clusters column, then export group-level expression and UMAP.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(Seurat)
})

project <- Sys.getenv("MG_LPC_ROOT")
input_object <- file.path(
  project,
  "objects",
  "simple_microglia",
  "mouse_GSE148676_simple_microglia_seurat_v5.rds"
)
output_object <- file.path(
  project,
  "objects",
  "simple_microglia",
  "mouse_GSE148676_merged_groups_seurat_v5.rds"
)
report_dir <- file.path(project, "reports", "simple_microglia")
figure_dir <- file.path(project, "figures", "simple_microglia")

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

merged_3_9_10 <- c("3", "9", "10")
merged_0_1_2_4_5_8 <- c("0", "1", "2", "4", "5", "8")
all_expected_clusters <- as.character(0:17)
remaining_clusters <- setdiff(
  all_expected_clusters,
  c(merged_3_9_10, merged_0_1_2_4_5_8)
)
new_cluster_levels <- as.character(0:10)
old_to_new_cluster <- c(
  # Final IDs are ordered by merged-group cell count, from largest to smallest.
  setNames(rep("1", length(merged_3_9_10)), merged_3_9_10),
  setNames(rep("0", length(merged_0_1_2_4_5_8)), merged_0_1_2_4_5_8),
  setNames(as.character(2:10), remaining_clusters)
)

message("Reading ", input_object)
microglia <- readRDS(input_object)
DefaultAssay(microglia) <- "RNA"

original_cluster <- as.character(microglia$seurat_clusters)
unexpected_clusters <- setdiff(unique(original_cluster), all_expected_clusters)
if (length(unexpected_clusters) > 0) {
  stop(
    "Unexpected original cluster labels: ",
    paste(unexpected_clusters, collapse = ", ")
  )
}

merged_cluster_source <- original_cluster
merged_cluster_source[original_cluster %in% merged_3_9_10] <- "3_9_10"
merged_cluster_source[original_cluster %in% merged_0_1_2_4_5_8] <- "0_1_2_4_5_8"
microglia$merged_cluster <- factor(
  unname(old_to_new_cluster[original_cluster]),
  levels = new_cluster_levels
)
microglia$merged_cluster_source <- merged_cluster_source
Idents(microglia) <- "merged_cluster"
  }))
))

# S4F  LPC myelin1 signature score by merged clusters; cluster1 compared to each other cluster
# Kind: computational
mg_register_panel("S4F", list(
  list(source="source/external/LPC/scripts/34_gse148676_reference_spatial_module.R", ranges="1-128;132-150;165-166", requested_ranges="1-128;132-150;165-166", role="complete input, signature, palette and violin context", code=quote({
# MODULE : Cell Type Annotation and Module Scoring -------------------
#
# Direct adaptation of spatial_module.r for merged GSE148676 cells.
# The original package calls and plotting structure are retained. Only the
# input/output paths and the GSE148676 cluster mapping differ. The saveRDS()
# directly after AddModuleScore makes the original reload-for-plot step work
# when this script is sourced from a fresh R session.

## 1 Parameters --------------------------------------------------------------
# Load dependencies
library(Seurat)
library(qs)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SCP)

# File paths
project <- Sys.getenv("MG_LPC_ROOT")
file_in_merfish <- file.path(
  project,
  "objects/simple_microglia/mouse_GSE148676_merged_groups_seurat_v5.rds"
)
file_in_markers <- file.path(
  project,
  "raw_data/reference/mouse_mg_markers.qs"
)
file_out_merfish <- file.path(
  project,
  "objects/simple_microglia/mouse_GSE148676_merged_groups_reference_module_score.rds"
)
file_out_umap_pdf <- file.path(
  project,
  "figures/simple_microglia/mouse_GSE148676_reference_module_score_UMAP.pdf"
)
file_out_umap_png <- sub("\\.pdf$", ".png", file_out_umap_pdf)
file_out_violin_pdf <- file.path(
  project,
  "figures/simple_microglia/mouse_GSE148676_reference_module_score_violin.pdf"
)
file_out_violin_png <- sub("\\.pdf$", ".png", file_out_violin_pdf)
dir.create(dirname(file_out_merfish), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(file_out_umap_pdf), recursive = TRUE, showWarnings = FALSE)

# Cell-type mapping dictionary (old merged-cluster label -> plot label)
col_old_celltype <- "merged_cluster"
col_new_celltype <- "plot_celltype"
celltype_dict <- c(
  "0" = "cluster0",
  "1" = "cluster1",
  "2" = "cluster2",
  "3" = "cluster3",
  "4" = "cluster4",
  "5" = "cluster5",
  "6" = "cluster6",
  "7" = "cluster7",
  "8" = "cluster8",
  "9" = "cluster9",
  "10" = "cluster10"
)

# Marker extraction and module-score parameters
marker_cluster_col <- "cluster"
target_cluster <- "cluster2"
marker_gene_col <- "gene"

score_module_name <- "myelin"
score_feature_name <- paste0(score_module_name, "1")
score_ctrl_size <- 10

# Plotting parameter retained from the reference script
plot_pt_size <- 0


## 2 Process Data ------------------------------------------------------------
# 2.1 Read data
merfish <- readRDS(file_in_merfish)
mouse_mg_markers <- qread(file_in_markers)

# 2.2 Replace cell-type names
merfish@meta.data[[col_new_celltype]] <- as.character(
  merfish@meta.data[[col_old_celltype]]
)
for (old_name in names(celltype_dict)) {
  new_name <- celltype_dict[[old_name]]
  merfish@meta.data[[col_new_celltype]][
    merfish@meta.data[[col_new_celltype]] %in% old_name
  ] <- new_name
}

# 2.3 Extract target-cluster genes and intersect with GSE148676 genes
filtered_data <- mouse_mg_markers %>%
  filter(.data[[marker_cluster_col]] == target_cluster)
int_gene <- intersect(filtered_data[[marker_gene_col]], rownames(merfish))

# 2.4 Calculate module score
merfish <- AddModuleScore(
  object = merfish,
  features = list(int_gene),
  ctrl = score_ctrl_size,
  name = score_module_name
)

# The reference SCP function accesses the Seurat v4 RNA@data slot. This
# conversion preserves the existing normalized RNA layer for plotting only.
merfish[["RNA"]] <- as(merfish[["RNA"]], "Assay")

# The reference script re-opens a previously saved scored object for plotting.
# Save here so the same operation is valid in a new session.
saveRDS(merfish, file = file_out_merfish)


## 3 Generate Plots ----------------------------------------------------------
merfish <- readRDS(file_out_merfish)

st_mg_colors_npg <- c(
  "cluster1" = "#E64B35",
  "cluster0" = "#4DBBD5",
  "cluster2" = "#00A087",
  "cluster3" = "#3C5488",
  "cluster4" = "#F39B7F",
  "cluster5" = "#7E6148",
  "cluster6" = "#8E6BBE",
  "cluster7" = "#7F7F7F",
  "cluster8" = "#66C2A5",
  "cluster9" = "#E6AB02",
  "cluster10" = "#A65628"
)

# The first group is the merged candidate group (old clusters 3, 9, and 10),
# preserving the reference comparison pattern of first group versus all others.
target_order <- c(
  "cluster1", "cluster0", "cluster2", "cluster3", "cluster4", "cluster5",
  "cluster6", "cluster7", "cluster8", "cluster9", "cluster10"
)
merfish$plot_celltype <- factor(merfish$plot_celltype, levels = target_order)

celltype <- levels(merfish$plot_celltype)
my_comparisons <- lapply(celltype[-1], function(x) c(x, celltype[1]))

p2 <- FeatureStatPlot(
  srt = merfish, group.by = "plot_celltype",
  stat.by = c(score_feature_name), add_box = FALSE,
  comparisons = my_comparisons,
  bg.by = "plot_celltype",
  palcolor = st_mg_colors_npg,
  bg_palcolor = st_mg_colors_npg
)

ggsave(filename = file_out_violin_pdf, plot = p2, width = 7, height = 5)
ggsave(filename = file_out_violin_png, plot = p2, width = 7, height = 5, dpi = 300)
  }))
))

# S4G  LPC Trim2/Rnf220/Mbp/Mog/Olig2/St18 normalized-expression UMAP panels
# Kind: computational
mg_register_panel("S4G", list(
  list(source="source/external/LPC/scripts/28_gse148676_extended_marker_featureplots.R", ranges="12-57;63-130;197-228", requested_ranges="12-57;63-130;197-228", role="input, case-insensitive matching, stored UMAP and RNA/data; grey background plus positive-expression yellow/red q99-clipped coord_equal panels; separate marker PNG/PDF output", code=quote({
project <- Sys.getenv("MG_LPC_ROOT")
object_path <- file.path(
  project,
  "objects",
  "simple_microglia",
  "mouse_GSE148676_simple_microglia_seurat_v5.rds"
)
figure_dir <- file.path(project, "figures", "simple_microglia")
individual_figure_dir <- file.path(
  figure_dir,
  "mouse_GSE148676_individual_featureplots"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_figure_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading ", object_path)
microglia <- readRDS(object_path)
DefaultAssay(microglia) <- "RNA"

requested_genes <- c("Trim2", "Rnf220", "St18", "Olig2", "Mbp", "Mog")
available_genes <- rownames(microglia)
matched_genes <- vapply(
  requested_genes,
  function(gene) {
    hits <- available_genes[toupper(available_genes) == toupper(gene)]
    if (length(hits) == 0) {
      stop("Gene is absent from the RNA assay: ", gene)
    }
    hits[[1]]
  },
  character(1)
)
names(matched_genes) <- requested_genes
message(
  "Matched features: ",
  paste(names(matched_genes), matched_genes, sep = "=", collapse = ", ")
)

umap_frame <- as.data.frame(Embeddings(microglia, "umap"))
colnames(umap_frame)[1:2] <- c("UMAP_1", "UMAP_2")
cell_order <- rownames(umap_frame)

normalized_matrix <- LayerData(
  microglia[["RNA"]],
  layer = "data"
)[matched_genes, cell_order, drop = FALSE]

for (gene_label in requested_genes) {
  feature_name <- matched_genes[[gene_label]]
  expression <- as.numeric(normalized_matrix[feature_name, ])
  expression_values[[gene_label]] <- expression

  plot_frame <- umap_frame
  plot_frame$expression <- expression
  positive_frame <- plot_frame[plot_frame$expression > 0, , drop = FALSE]
  positive_frame <- positive_frame[
    order(positive_frame$expression),
    ,
    drop = FALSE
  ]

  if (nrow(positive_frame) == 0) {
    upper_limit <- 1
  } else {
    upper_limit <- as.numeric(
      quantile(
        positive_frame$expression,
        0.99,
        names = FALSE
      )
    )
  }

  gene_plots[[gene_label]] <- ggplot(
    plot_frame,
    aes(x = UMAP_1, y = UMAP_2)
  ) +
    geom_point(
      color = "#B8B8B8",
      size = point_size,
      alpha = 0.60,
      shape = 16,
      stroke = 0
    ) +
    geom_point(
      data = positive_frame,
      aes(color = expression),
      size = point_size * 1.55,
      alpha = 1,
      shape = 16,
      stroke = 0
    ) +
    scale_color_gradient(
      low = "#FFD166",
      high = "#B2182B",
      limits = c(0, upper_limit),
      oob = scales::squish
    ) +
    labs(
      title = gene_label,
      color = "Log-normalized\nexpression",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    coord_equal() +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
}

# Export every marker panel as a standalone figure. These files retain the
# same UMAP coordinates and normalized RNA data used in the combined panel.
for (gene_label in requested_genes) {
  individual_path <- file.path(
    individual_figure_dir,
    paste0(
      "mouse_GSE148676_",
      gene_label,
      "_log_normalized_FeaturePlot.png"
    )
  )
  ggsave(
    individual_path,
    gene_plots[[gene_label]],
    width = 7.2,
    height = 6.4,
    dpi = 300,
    bg = "white"
  )
  message("Wrote ", individual_path)

  individual_pdf_path <- sub("\\.png$", ".pdf", individual_path)
  ggsave(
    individual_pdf_path,
    gene_plots[[gene_label]],
    width = 7.2,
    height = 6.4,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  message("Wrote ", individual_pdf_path)
}
  }))
))

# S4H  100% retained-MG composition, LPC OLM and Monkey/Mouse MG0-like
# Kind: computational
mg_register_panel("S4H", list(
  list(source="source/project/Reference/combine/scripts/MultiSpecies/results/MG0_like_proportion_LPC_GSE148676_Monkey_Mouse_20260921/100pct_composition_and_umap_source_audit/plot_LPC_MG0_Monkey_Mouse_100pct_composition_umap_audit_20260921.R", ranges="37-109;113-138;175-225;259-268", requested_ranges="37-109;113-138;175-225;259-268", role="LPC metadata and cross-species input/target definition; counts, all-retained-MG denominators and percentages; 100% stacked bars and annotations; PDF/PNG exports", code=quote({
lpc_metadata_path <- file.path(
  lpc_project_root,
  "reports/simple_microglia/mouse_GSE148676_merged_groups_cell_metadata.csv"
)
multi_species_object_path <- file.path(
  project_root,
  "combine/scripts/MultiSpecies/data/processed_data/02.cross_species_integrated_noref_Human_Mouse_Monkey.qs"
)
if (!file.exists(lpc_metadata_path)) stop("LPC metadata/UMAP export is absent: ", lpc_metadata_path)
if (!file.exists(multi_species_object_path)) stop("Multi-species object is absent: ", multi_species_object_path)

# LPC (GSE148676): use the existing per-cell export from the merged-cluster
# microglia object. Its UMAP is the stored UMAP and cluster 1 is the project's
# designated M0/OLM group.
lpc_metadata <- read.csv(lpc_metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
required_lpc_columns <- c("cell", "merged_cluster", "UMAP_1", "UMAP_2")
if (!all(required_lpc_columns %in% colnames(lpc_metadata))) {
  stop("LPC metadata is missing: ", paste(setdiff(required_lpc_columns, colnames(lpc_metadata)), collapse = ", "))
}
if (anyDuplicated(lpc_metadata$cell) > 0L) stop("LPC metadata contains duplicated cell IDs.")
lpc_frame <- data.frame(
  Dataset = "LPC (GSE148676, mouse)",
  UMAP_1 = as.numeric(lpc_metadata$UMAP_1),
  UMAP_2 = as.numeric(lpc_metadata$UMAP_2),
  Target_group = ifelse(as.character(lpc_metadata$merged_cluster) == "1", "OLM", "Other MG"),
  stringsAsFactors = FALSE
)
lpc_total <- nrow(lpc_frame)
lpc_target <- sum(lpc_frame$Target_group == "OLM")
lpc_cluster_values <- paste(sort(unique(as.character(lpc_metadata$merged_cluster))), collapse = ", ")
if (lpc_total == 0L || lpc_target == 0L) stop("LPC has no cells and/or no merged cluster 1 cells.")
rm(lpc_metadata)
gc(verbose = FALSE)

# Monkey and Mouse: retain the stored integrated UMAP and use the original
# author-defined integrated cluster 1 as MG0-like for each species.
multi_species <- qread(multi_species_object_path)
multi_metadata <- multi_species[[]]
required_multi_columns <- c("Species", "integrated_snn_res.0.1")
if (!all(required_multi_columns %in% colnames(multi_metadata))) {
  stop("Multi-species metadata is missing: ", paste(setdiff(required_multi_columns, colnames(multi_metadata)), collapse = ", "))
}
if (!"umap" %in% Reductions(multi_species)) stop("The multi-species object has no stored UMAP reduction.")
multi_umap <- Embeddings(multi_species, reduction = "umap")
if (!identical(rownames(multi_metadata), rownames(multi_umap))) {
  stop("Multi-species metadata and stored UMAP have non-identical cell order.")
}

make_multispecies_frame <- function(species_name) {
  keep <- as.character(multi_metadata$Species) == species_name
  frame <- data.frame(
    Dataset = species_name,
    UMAP_1 = as.numeric(multi_umap[keep, 1]),
    UMAP_2 = as.numeric(multi_umap[keep, 2]),
    Target_group = ifelse(
      as.character(multi_metadata$integrated_snn_res.0.1[keep]) == "1",
      "MG0-like", "Other MG"
    ),
    stringsAsFactors = FALSE
  )
  if (nrow(frame) == 0L || sum(frame$Target_group == "MG0-like") == 0L) {
    stop(species_name, " has no cells and/or no integrated cluster 1 cells.")
  }
  frame
}
monkey_frame <- make_multispecies_frame("Monkey")
mouse_frame <- make_multispecies_frame("Mouse")
monkey_cluster_values <- paste(sort(unique(as.character(
  multi_metadata$integrated_snn_res.0.1[as.character(multi_metadata$Species) == "Monkey"]
))), collapse = ", ")
mouse_cluster_values <- paste(sort(unique(as.character(
  multi_metadata$integrated_snn_res.0.1[as.character(multi_metadata$Species) == "Mouse"]
))), collapse = ", ")

composition_table <- data.frame(
  Dataset = c("LPC (GSE148676, mouse)", "Monkey", "Mouse"),
  Target_label = c("OLM", "MG0-like", "MG0-like"),
  Target_definition = c(
    "merged_cluster == 1",
    "integrated_snn_res.0.1 == 1",
    "integrated_snn_res.0.1 == 1"
  ),
  Target_cells = c(
    lpc_target,
    sum(monkey_frame$Target_group == "MG0-like"),
    sum(mouse_frame$Target_group == "MG0-like")
  ),
  Total_microglia_cells = c(nrow(lpc_frame), nrow(monkey_frame), nrow(mouse_frame)),
  stringsAsFactors = FALSE
)
composition_table$Target_percent <- 100 * composition_table$Target_cells / composition_table$Total_microglia_cells
composition_table$Other_MG_cells <- composition_table$Total_microglia_cells - composition_table$Target_cells
composition_table$Other_MG_percent <- 100 - composition_table$Target_percent
composition_table$Dataset <- factor(
  composition_table$Dataset,
  # coord_flip displays the last factor level at the top. This ordering gives
  # the requested LPC, Monkey, Mouse top-to-bottom reading sequence.
  levels = c("Mouse", "Monkey", "LPC (GSE148676, mouse)")
)
write.csv(composition_table, file.path(result_dir, "MG0_like_100pct_composition_counts.csv"), row.names = FALSE)

composition_plot <- ggplot(composition_long, aes(x = Dataset, y = Percent, fill = Segment)) +
  # The reference composition plot places the target subtype at the right end
  # of the 100% bar, so Other MG is stacked first and the target follows it.
  geom_col(
    width = 0.72,
    color = "#1A1A1A",
    linewidth = 0.35,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    data = composition_table,
    aes(x = Dataset, y = Other_MG_percent / 2, label = sprintf("Other MG %.1f%%", Other_MG_percent)),
    inherit.aes = FALSE,
    size = 4.0,
    color = "#303030"
  ) +
  geom_text(
    data = composition_table,
    aes(x = Dataset, y = 103.0, label = target_annotation),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.9,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("Other MG" = "#D0D0D0", "OLM" = "#4DAF4A", "MG0-like" = "#8FCB6A"),
    breaks = c("OLM", "MG0-like", "Other MG"),
    name = "Microglial group"
  ) +
  scale_y_continuous(
    limits = c(0, 115),
    breaks = seq(0, 100, by = 25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_flip(clip = "off") +
  labs(
    title = "MG0/MG0-like composition among all microglia",
    subtitle = "Every bar is normalized to 100% of the retained microglial cells in its input dataset",
    x = NULL,
    y = "Proportion"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.y = element_text(face = "bold", size = 11),
    axis.title.x = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    plot.margin = margin(t = 10, r = 165, b = 8, l = 8, unit = "pt")
  )

ggsave(
  file.path(result_dir, "MG0_like_100pct_composition.pdf"),
  composition_plot,
  width = 10.8, height = 4.9, device = grDevices::cairo_pdf, bg = "white"
)
ggsave(
  file.path(result_dir, "MG0_like_100pct_composition.png"),
  composition_plot,
  width = 10.8, height = 4.9, dpi = 300, bg = "white"
)
  }))
))
