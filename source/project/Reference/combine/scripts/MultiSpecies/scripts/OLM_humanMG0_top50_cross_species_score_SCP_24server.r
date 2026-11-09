# Human OLM (MG0) top-50 signature score in the integrated human/monkey/mouse object
#
# The module-score and violin-plot parameters below deliberately follow the user's
# reference script `spatial_module.r`:
#   AddModuleScore(features = list(...), ctrl = 10, name = ...)
#   SCP::FeatureStatPlot(..., add_box = FALSE, bg.by = group.by,
#                        palcolor = ..., bg_palcolor = ...)
# This script is a new analysis file.  It does not modify `combine_3.r` or its object.

scp_lib <- file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library")
.libPaths(c(scp_lib, .libPaths()))

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(SCP)
  library(ggplot2)
})

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies")
input_qs <- file.path(project_dir, "data/processed_data/02.cross_species_integrated_noref_Human_Mouse_Monkey.qs")
figure_dir <- file.path(project_dir, "results/figure/OLM_humanMG0_top50_score_SCP")
table_dir <- file.path(project_dir, "results/table/OLM_humanMG0_top50_score_SCP")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, stem, width = 7, height = 6) {
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 300, bg = "white"
  )
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, bg = "white"
  )
}

obj <- qread(input_qs)
DefaultAssay(obj) <- "RNA"
stopifnot(all(c("Species", "CellType1", "integrated_snn_res.0.1") %in% colnames(obj@meta.data)))

# 1. Find the human OLM/MG0 top-50 genes.  The contrast is human MG0 vs human MG1--MG5.
human_obj <- subset(obj, subset = Species == "Human")
DefaultAssay(human_obj) <- "RNA"
Idents(human_obj) <- "CellType1"
human_mg0_de <- FindMarkers(
  human_obj,
  ident.1 = "MG0",
  ident.2 = c("MG1", "MG2", "MG3", "MG4", "MG5"),
  only.pos = TRUE,
  logfc.threshold = 0,
  min.pct = 0.10,
  test.use = "wilcox"
)
human_mg0_de$gene <- rownames(human_mg0_de)
human_mg0_de <- human_mg0_de[order(human_mg0_de$p_val_adj, -human_mg0_de$avg_log2FC), , drop = FALSE]
human_mg0_top50 <- head(human_mg0_de, 50)
score_genes <- human_mg0_top50$gene
missing_score_genes <- setdiff(score_genes, rownames(obj))
stopifnot(length(missing_score_genes) == 0L, length(score_genes) == 50L)
write.csv(human_mg0_top50, file.path(table_dir, "human_MG0_vs_MG1toMG5_FindMarkers_top50.csv"), row.names = FALSE)

# 2. Score every cell in the same integrated UMAP object.
# `ctrl = 10` is the reference-script setting.  The seed makes control-gene sampling reproducible.
score_module_name <- "Human_OLM_MG0_top50_score"
score_feature_name <- paste0(score_module_name, "1")
score_ctrl_size <- 10
set.seed(20260915)
obj <- AddModuleScore(
  object = obj,
  features = list(score_genes),
  ctrl = score_ctrl_size,
  name = score_module_name
)

# Recreate the historical cross-species OLM labels, without changing the source object.
obj$OLM_cross_species_label <- as.character(obj$CellType1)
obj$OLM_cross_species_label[obj$Species == "Mouse" & obj$integrated_snn_res.0.1 == "1"] <- "Mouse_MG0_like"
obj$OLM_cross_species_label[obj$Species == "Monkey" & obj$integrated_snn_res.0.1 == "1"] <- "Monkey_MG0_like"
obj$OLM_cross_species_label[!obj$OLM_cross_species_label %in% c("MG0", "Mouse_MG0_like", "Monkey_MG0_like")] <- "Other_MG"

# 3. Feature plots: identical integrated UMAP coordinates; score colour is grey to red.
score_cols <- c("grey85", "#D73027")
plot_score_umap <- function(cells, label, stem) {
  p <- FeaturePlot(
    obj,
    features = score_feature_name,
    reduction = "umap",
    cells = cells,
    cols = score_cols,
    pt.size = 0.15,
    order = TRUE
  ) + ggtitle(label) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(p, stem, width = 7.5, height = 6)
}

plot_score_umap(colnames(obj), "Human OLM/MG0 top-50 module score: all species", "01_score_UMAP_all_species")
plot_score_umap(colnames(obj)[obj$Species == "Human"], "Human OLM/MG0 top-50 module score: Human", "02_score_UMAP_human")
plot_score_umap(colnames(obj)[obj$Species == "Monkey"], "Human OLM/MG0 top-50 module score: Monkey", "03_score_UMAP_monkey")
plot_score_umap(colnames(obj)[obj$Species == "Mouse"], "Human OLM/MG0 top-50 module score: Mouse", "04_score_UMAP_mouse")

# 4. Reference-style SCP FeatureStatPlot comparisons.
# Every species is a distinct complete plot: its MG0-like group versus all other MG cells of that species.
make_species_scp_plot <- function(species, target_label, target_group, stem, target_colour) {
  species_cells <- colnames(obj)[obj$Species == species]
  srt <- subset(obj, cells = species_cells)
  group_col <- paste0("score_group_", tolower(species))
  srt[[group_col]] <- factor(
    ifelse(srt$OLM_cross_species_label == target_label, target_group, "Other MG"),
    levels = c("Other MG", target_group)
  )
  cols <- c("Other MG" = "#D9D9D9")
  cols[target_group] <- target_colour
  p <- SCP::FeatureStatPlot(
    srt = srt,
    group.by = group_col,
    stat.by = c(score_feature_name),
    add_box = FALSE,
    comparisons = list(c(target_group, "Other MG")),
    bg.by = group_col,
    palcolor = cols,
    bg_palcolor = cols
  )
  save_plot(p, stem, width = 5.5, height = 6)
}

make_species_scp_plot("Human", "MG0", "Human MG0", "05_SCP_score_human_MG0_vs_other_MG", "#4DAF4A")
make_species_scp_plot("Monkey", "Monkey_MG0_like", "Monkey MG0-like", "06_SCP_score_monkey_MG0like_vs_other_MG", "#33A02C")
make_species_scp_plot("Mouse", "Mouse_MG0_like", "Mouse MG0-like", "07_SCP_score_mouse_MG0like_vs_other_MG", "#B2DF8A")

# 5. One horizontal, descriptive cross-species comparison requested by the user.
# This is not a species-adjusted test: it visualizes score distributions only.
obj$score_group_cross_species <- "Other MG"
obj$score_group_cross_species[obj$Species == "Human" & obj$OLM_cross_species_label == "MG0"] <- "Human MG0"
obj$score_group_cross_species[obj$Species == "Monkey" & obj$OLM_cross_species_label == "Monkey_MG0_like"] <- "Monkey MG0-like"
obj$score_group_cross_species[obj$Species == "Mouse" & obj$OLM_cross_species_label == "Mouse_MG0_like"] <- "Mouse MG0-like"
obj$score_group_cross_species <- factor(
  obj$score_group_cross_species,
  levels = c("Other MG", "Human MG0", "Monkey MG0-like", "Mouse MG0-like")
)
cross_groups <- levels(obj$score_group_cross_species)
cross_cols <- c(
  "Other MG" = "#D9D9D9",
  "Human MG0" = "#4DAF4A",
  "Monkey MG0-like" = "#33A02C",
  "Mouse MG0-like" = "#B2DF8A"
)
cross_comparisons <- lapply(cross_groups[-1], function(group) c(group, "Other MG"))
p_cross_species <- SCP::FeatureStatPlot(
  srt = obj,
  group.by = "score_group_cross_species",
  stat.by = c(score_feature_name),
  add_box = FALSE,
  comparisons = cross_comparisons,
  bg.by = "score_group_cross_species",
  palcolor = cross_cols,
  bg_palcolor = cross_cols
)
save_plot(p_cross_species, "08_SCP_score_cross_species_MG0like_vs_other_MG", width = 8.5, height = 6)

# 6. Audit tables for the score and group definitions.
score_audit <- obj@meta.data[, c("Species", "CellType1", "integrated_snn_res.0.1", "OLM_cross_species_label", "score_group_cross_species", score_feature_name), drop = FALSE]
score_audit$cell_barcode <- rownames(score_audit)
score_audit <- score_audit[, c("cell_barcode", setdiff(colnames(score_audit), "cell_barcode")), drop = FALSE]
write.csv(score_audit, file.path(table_dir, "per_cell_OLM_top50_score_and_group_audit.csv"), row.names = FALSE)
score_summary <- aggregate(
  score_audit[[score_feature_name]],
  by = list(Species = score_audit$Species, Group = score_audit$score_group_cross_species),
  FUN = function(x) c(n = length(x), mean = mean(x), median = median(x), sd = sd(x))
)
score_summary <- data.frame(
  Species = score_summary$Species,
  Group = score_summary$Group,
  do.call(rbind, score_summary$x),
  row.names = NULL,
  check.names = FALSE
)
write.csv(score_summary, file.path(table_dir, "OLM_top50_score_group_summary.csv"), row.names = FALSE)

qsave(obj, file.path(table_dir, "02_cross_species_integrated_with_OLM_top50_score.qs"), preset = "high")
writeLines(c(
  "Human OLM/MG0 top-50 score analysis completed.",
  paste("Top-50 genes:", paste(score_genes, collapse = ", ")),
  "AddModuleScore ctrl = 10; seed = 20260915.",
  "SCP FeatureStatPlot settings follow spatial_module.r."
), file.path(table_dir, "README.txt"))

message("Completed. Figures: ", figure_dir)
message("Audit tables and scored object: ", table_dir)
