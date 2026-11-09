#!/usr/bin/env Rscript

# 100% MG composition and source-audit UMAPs for the requested LPC, Monkey and
# Mouse inputs. All labels and UMAP coordinates are read from existing outputs
# or objects; no source data, clustering, labels or embeddings are changed.

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
lpc_project_root <- Sys.getenv("MG_LPC_ROOT")
result_dir <- file.path(
  project_root,
  "combine/scripts/MultiSpecies/results/MG0_like_proportion_LPC_GSE148676_Monkey_Mouse_20260921",
  "100pct_composition_and_umap_source_audit"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# Prevent any R cache or temporary data from being written to a home directory.
tmp_dir <- file.path(result_dir, "tmp")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(
  TMPDIR = tmp_dir,
  TMP = tmp_dir,
  TEMP = tmp_dir,
  XDG_CACHE_HOME = file.path(result_dir, "cache"),
  R_USER_CACHE_DIR = file.path(result_dir, "cache", "R")
)
dir.create(Sys.getenv("XDG_CACHE_HOME"), recursive = TRUE, showWarnings = FALSE)
dir.create(Sys.getenv("R_USER_CACHE_DIR"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library"), file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

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
rm(multi_species, multi_metadata, multi_umap)
gc(verbose = FALSE)

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

source_validation_table <- data.frame(
  Dataset = c("LPC (GSE148676, mouse)", "Monkey", "Mouse"),
  Input_artifact = c(lpc_metadata_path, multi_species_object_path, multi_species_object_path),
  Input_contents = c(
    "52,253 retained LPC microglial cells with merged_cluster and stored UMAP coordinates",
    "Monkey subset of the stored Human-Monkey-Mouse integrated microglia object",
    "Mouse subset of the stored Human-Monkey-Mouse integrated microglia object"
  ),
  UMAP_used = c("Stored LPC UMAP exported from the merged-cluster object", "Stored integrated UMAP", "Stored integrated UMAP"),
  Available_cluster_values = c(lpc_cluster_values, monkey_cluster_values, mouse_cluster_values),
  Target_label = as.character(composition_table$Target_label),
  Target_definition = composition_table$Target_definition,
  Target_cells = composition_table$Target_cells,
  Total_microglia_cells = composition_table$Total_microglia_cells,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write.csv(source_validation_table, file.path(result_dir, "source_dataset_validation.csv"), row.names = FALSE)

# Each bar sums to 100%; the right-side annotation supplies the exact target
# count and target percentage, making both the numerator and denominator auditable.
composition_long <- rbind(
  data.frame(Dataset = composition_table$Dataset, Segment = "Other MG", Percent = composition_table$Other_MG_percent),
  data.frame(Dataset = composition_table$Dataset, Segment = as.character(composition_table$Target_label), Percent = composition_table$Target_percent)
)
composition_long$Dataset <- factor(composition_long$Dataset, levels = levels(composition_table$Dataset))
composition_long$Segment <- factor(composition_long$Segment, levels = c("Other MG", "OLM", "MG0-like"))
composition_table$target_annotation <- sprintf(
  "%s: %.1f%% (%s/%s)",
  composition_table$Target_label,
  composition_table$Target_percent,
  format(composition_table$Target_cells, big.mark = ","),
  format(composition_table$Total_microglia_cells, big.mark = ",")
)

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

make_audit_umap <- function(frame, title, target_label) {
  other <- frame[frame$Target_group == "Other MG", , drop = FALSE]
  target <- frame[frame$Target_group == target_label, , drop = FALSE]
  ggplot() +
    geom_point(
      data = other,
      aes(x = UMAP_1, y = UMAP_2),
      color = "#C9C9C9", size = 0.22, alpha = 0.65
    ) +
    geom_point(
      data = target,
      aes(x = UMAP_1, y = UMAP_2, color = Target_group),
      size = 0.33, alpha = 0.95
    ) +
    scale_color_manual(values = c("OLM" = "#4DAF4A", "MG0-like" = "#8FCB6A"), name = "Target group") +
    coord_equal() +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      legend.position = "bottom",
      axis.title = element_text(size = 10)
    )
}

lpc_umap_plot <- make_audit_umap(lpc_frame, "LPC (GSE148676): OLM vs Other MG", "OLM")
monkey_umap_plot <- make_audit_umap(monkey_frame, "Monkey: MG0-like vs Other MG", "MG0-like")
mouse_umap_plot <- make_audit_umap(mouse_frame, "Mouse: MG0-like vs Other MG", "MG0-like")
umap_audit_plot <- (lpc_umap_plot + monkey_umap_plot + mouse_umap_plot) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

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
ggsave(
  file.path(result_dir, "source_data_audit_UMAPs.pdf"),
  umap_audit_plot,
  width = 14.4, height = 5.3, device = grDevices::cairo_pdf, bg = "white"
)
ggsave(
  file.path(result_dir, "source_data_audit_UMAPs.png"),
  umap_audit_plot,
  width = 14.4, height = 5.3, dpi = 300, bg = "white"
)

writeLines(
  c(
    "100% composition and UMAP source audit",
    "",
    "LPC/GSE148676: all retained microglia are the denominator; merged_cluster == 1 is labelled OLM.",
    "Monkey: all Monkey cells in the stored integrated Human-Monkey-Mouse microglia object are the denominator; integrated_snn_res.0.1 == 1 is labelled MG0-like.",
    "Mouse: all Mouse cells in the stored integrated Human-Monkey-Mouse microglia object are the denominator; integrated_snn_res.0.1 == 1 is labelled MG0-like.",
    "The UMAPs use stored coordinates. No UMAP, clustering, label or source object was recalculated or modified."
  ),
  file.path(result_dir, "README.txt")
)

print(composition_table[, c("Dataset", "Target_label", "Target_cells", "Total_microglia_cells", "Target_percent")])
message("Wrote 100% composition and UMAP audit results to ", result_dir)
