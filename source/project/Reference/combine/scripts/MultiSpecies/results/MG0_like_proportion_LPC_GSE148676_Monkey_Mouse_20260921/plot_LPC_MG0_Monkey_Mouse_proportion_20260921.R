#!/usr/bin/env Rscript

# Compare the proportion of the pre-defined MG0/MG0-like population among all
# microglia in three independent analysis inputs. The source data are read only.
# LPC/GSE148676: merged_cluster == "1" (the project's designated M0 cluster).
# Multi-species Monkey/Mouse: integrated_snn_res.0.1 == "1" (the original
# project's MG0-like definition for each non-human species).

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
lpc_project_root <- Sys.getenv("MG_LPC_ROOT")
result_dir <- file.path(
  project_root,
  "combine/scripts/MultiSpecies/results/MG0_like_proportion_LPC_GSE148676_Monkey_Mouse_20260921"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# Keep temporary files and R caches outside the home directory.
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
})

lpc_cell_table_path <- file.path(
  lpc_project_root,
  "reports/simple_microglia/mouse_GSE148676_M0_module_score_per_cell.csv"
)
multi_species_object_path <- file.path(
  project_root,
  "combine/scripts/MultiSpecies/data/processed_data/02.cross_species_integrated_noref_Human_Mouse_Monkey.qs"
)

if (!file.exists(lpc_cell_table_path)) stop("LPC cell table does not exist: ", lpc_cell_table_path)
if (!file.exists(multi_species_object_path)) stop("Multi-species object does not exist: ", multi_species_object_path)

# The LPC cell table is the cell-level export made directly from the project's
# mouse_GSE148676_merged_groups_M0_module_score.rds object. It contains one row
# per retained microglial cell and its existing merged_cluster label.
lpc_cells <- read.csv(lpc_cell_table_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("cell", "merged_cluster") %in% colnames(lpc_cells))) {
  stop("The LPC cell table must contain cell and merged_cluster columns.")
}
if (anyDuplicated(lpc_cells$cell) > 0L) stop("The LPC cell table contains duplicated cell identifiers.")
lpc_total <- nrow(lpc_cells)
lpc_target <- sum(as.character(lpc_cells$merged_cluster) == "1")
if (lpc_total == 0L || lpc_target == 0L) stop("LPC has no total cells and/or no merged cluster 1 cells.")
rm(lpc_cells)
gc(verbose = FALSE)

# Use the unmodified integrated object and the original cluster-1 rule rather
# than the later descriptive OLM alias column.
multi_species <- qread(multi_species_object_path)
multi_metadata <- multi_species[[]]
required_multi_columns <- c("Species", "integrated_snn_res.0.1")
if (!all(required_multi_columns %in% colnames(multi_metadata))) {
  stop("The multi-species object lacks: ", paste(setdiff(required_multi_columns, colnames(multi_metadata)), collapse = ", "))
}

species_labels <- as.character(multi_metadata$Species)
cluster_labels <- as.character(multi_metadata$integrated_snn_res.0.1)
count_species <- function(species_name) {
  keep <- species_labels == species_name
  total <- sum(keep)
  target <- sum(keep & cluster_labels == "1")
  if (total == 0L || target == 0L) {
    stop(species_name, " has no total cells and/or no integrated cluster 1 cells.")
  }
  c(total = total, target = target)
}
monkey_counts <- count_species("Monkey")
mouse_counts <- count_species("Mouse")
rm(multi_species, multi_metadata)
gc(verbose = FALSE)

composition_table <- data.frame(
  Dataset = c("LPC (GSE148676, mouse)", "Monkey", "Mouse"),
  Target_label = c("MG0", "MG0-like", "MG0-like"),
  Target_definition = c(
    "LPC: merged_cluster == 1",
    "Monkey: integrated_snn_res.0.1 == 1",
    "Mouse: integrated_snn_res.0.1 == 1"
  ),
  Target_cells = c(lpc_target, monkey_counts[["target"]], mouse_counts[["target"]]),
  Total_microglia_cells = c(lpc_total, monkey_counts[["total"]], mouse_counts[["total"]]),
  Percent_of_microglia = c(
    100 * lpc_target / lpc_total,
    100 * monkey_counts[["target"]] / monkey_counts[["total"]],
    100 * mouse_counts[["target"]] / mouse_counts[["total"]]
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
composition_table$Dataset <- factor(
  composition_table$Dataset,
  levels = c("LPC (GSE148676, mouse)", "Monkey", "Mouse")
)
composition_table$display_label <- sprintf(
  "%.1f%%\n(%s / %s)",
  composition_table$Percent_of_microglia,
  format(composition_table$Target_cells, big.mark = ","),
  format(composition_table$Total_microglia_cells, big.mark = ",")
)

write.csv(
  composition_table[, setdiff(colnames(composition_table), "display_label"), drop = FALSE],
  file = file.path(result_dir, "MG0_like_proportion_counts.csv"),
  row.names = FALSE
)

upper_limit <- max(20, ceiling((max(composition_table$Percent_of_microglia) + 2) / 5) * 5)
bar_plot <- ggplot(composition_table, aes(x = Dataset, y = Percent_of_microglia, fill = Dataset)) +
  geom_col(width = 0.68, color = "#202020", linewidth = 0.35, show.legend = FALSE) +
  geom_text(
    aes(label = display_label),
    vjust = -0.25,
    lineheight = 0.95,
    size = 4.1,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    "LPC (GSE148676, mouse)" = "#D55E00",
    "Monkey" = "#0072B2",
    "Mouse" = "#009E73"
  )) +
  scale_y_continuous(
    limits = c(0, upper_limit),
    breaks = seq(0, upper_limit, by = 5),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "MG0/MG0-like proportion among microglia",
    subtitle = "LPC: designated merged cluster 1; Monkey and Mouse: original integrated cluster 1",
    x = NULL,
    y = "MG0/MG0-like cells (% of all microglia)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(face = "bold", size = 11),
    axis.title.y = element_text(face = "bold", size = 11),
    axis.line = element_line(linewidth = 0.6),
    plot.margin = margin(t = 10, r = 12, b = 10, l = 10, unit = "pt")
  )

ggsave(
  filename = file.path(result_dir, "MG0_like_proportion_among_microglia.pdf"),
  plot = bar_plot,
  width = 8.4,
  height = 5.6,
  device = grDevices::cairo_pdf,
  bg = "white"
)
ggsave(
  filename = file.path(result_dir, "MG0_like_proportion_among_microglia.png"),
  plot = bar_plot,
  width = 8.4,
  height = 5.6,
  dpi = 300,
  bg = "white"
)

writeLines(
  c(
    "MG0/MG0-like composition provenance",
    "",
    paste("LPC cell-level input:", lpc_cell_table_path),
    "LPC target: merged_cluster == 1; denominator: all 52,253 retained LPC microglial cells.",
    paste("Multi-species input:", multi_species_object_path),
    "Monkey target: Species == Monkey and integrated_snn_res.0.1 == 1; denominator: all Monkey microglia.",
    "Mouse target: Species == Mouse and integrated_snn_res.0.1 == 1; denominator: all Mouse microglia.",
    "No cell labels or source objects were changed."
  ),
  con = file.path(result_dir, "README.txt")
)

print(composition_table[, c("Dataset", "Target_label", "Target_cells", "Total_microglia_cells", "Percent_of_microglia")])
message("Wrote results to ", result_dir)
