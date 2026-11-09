#!/usr/bin/env Rscript
# Descriptive QC plots for the original human cross-study MG0-MG5 object.
# Input object is read only. No QC filtering, relabelling, or resampling occurs.

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
combine_root <- file.path(project_root, "combine")
input_file <- file.path(combine_root, "data/processed_data/mg_cca_meta_v2.qs")
output_dir <- file.path(combine_root, "results/MG0_MG5_QC_cross_study")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Input object not found: ", input_file)
obj <- qread(input_file)
required <- c("CellType1", "nCount_RNA", "nFeature_RNA", "Study")
if (!all(required %in% colnames(obj@meta.data))) {
  stop("Required metadata missing: ", paste(setdiff(required, colnames(obj@meta.data)), collapse = ", "))
}

mg_order <- paste0("MG", 0:5)
mg_palette <- c(
  "MG0" = "#4DAF4A", "MG1" = "#A65628", "MG2" = "#377EB8",
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG5" = "#8B0000"
)

# Analysis universe: all cells that belong to the original cross-study MG0-MG5
# categories. Values are pre-existing cell-level metadata, not recalculated.
qc_data <- obj@meta.data %>%
  transmute(
    CellType = factor(as.character(CellType1), levels = mg_order),
    Study = as.character(Study),
    nCount_RNA = as.numeric(nCount_RNA),
    nFeature_RNA = as.numeric(nFeature_RNA)
  ) %>%
  filter(!is.na(CellType), !is.na(nCount_RNA), !is.na(nFeature_RNA))

if (!nrow(qc_data)) stop("No non-missing MG0-MG5 QC records available.")

# Per-cell pooled summary: answers the requested mean detected-gene statistic.
qc_summary <- qc_data %>%
  group_by(CellType) %>%
  summarise(
    cells = n(),
    mean_nFeature_RNA = mean(nFeature_RNA),
    median_nFeature_RNA = median(nFeature_RNA),
    mean_nCount_RNA = mean(nCount_RNA),
    median_nCount_RNA = median(nCount_RNA),
    .groups = "drop"
  )

# Study-resolved table: keeps cross-study depth heterogeneity auditable.
qc_by_study <- qc_data %>%
  group_by(Study, CellType) %>%
  summarise(
    cells = n(),
    mean_nFeature_RNA = mean(nFeature_RNA),
    median_nFeature_RNA = median(nFeature_RNA),
    q1_nFeature_RNA = quantile(nFeature_RNA, 0.25),
    q3_nFeature_RNA = quantile(nFeature_RNA, 0.75),
    mean_nCount_RNA = mean(nCount_RNA),
    median_nCount_RNA = median(nCount_RNA),
    q1_nCount_RNA = quantile(nCount_RNA, 0.25),
    q3_nCount_RNA = quantile(nCount_RNA, 0.75),
    .groups = "drop"
  )

# Requested bar chart: mean number of detected genes (nFeature_RNA).
p_mean_nfeature <- ggplot(qc_summary, aes(x = CellType, y = mean_nFeature_RNA, fill = CellType)) +
  geom_col(width = 0.75, color = "black") +
  geom_text(aes(label = paste0(round(mean_nFeature_RNA), "\n(n=", comma(cells), ")")),
            vjust = -0.25, size = 3.6) +
  scale_fill_manual(values = mg_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = comma) +
  labs(
    x = NULL, y = "Mean detected genes per cell (nFeature_RNA)",
    title = "Human cross-study MG0-MG5: mean detected genes"
  ) +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(color = "black"))

# Requested QC distributions: raw cell-level nCount and nFeature as separate
# violin plots. No data transformation is applied to the displayed values.
plot_violin <- function(y_col, y_label, title_text) {
  ggplot(qc_data, aes(x = CellType, y = .data[[y_col]], fill = CellType)) +
    geom_violin(scale = "width", trim = FALSE, color = "black", linewidth = 0.35) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.35) +
    scale_fill_manual(values = mg_palette) +
    scale_y_continuous(labels = comma) +
    labs(x = NULL, y = y_label, title = title_text) +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text = element_text(color = "black"))
}

p_ncount <- plot_violin(
  "nCount_RNA", "RNA counts per cell (nCount_RNA)",
  "Human cross-study MG0-MG5: nCount_RNA"
)
p_nfeature <- plot_violin(
  "nFeature_RNA", "Detected genes per cell (nFeature_RNA)",
  "Human cross-study MG0-MG5: nFeature_RNA"
)

save_plot <- function(plot, stem, width = 7, height = 5) {
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot, width = width, height = height, bg = "white")
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300, bg = "white")
}

save_plot(p_mean_nfeature, "MG0_MG5_mean_nFeature_RNA_bar")
save_plot(p_ncount, "MG0_MG5_nCount_RNA_violin")
save_plot(p_nfeature, "MG0_MG5_nFeature_RNA_violin")
write_csv(qc_summary, file.path(output_dir, "MG0_MG5_QC_pooled_summary.csv"))
write_csv(qc_by_study, file.path(output_dir, "MG0_MG5_QC_summary_by_Study.csv"))
writeLines(
  c(
    "Human cross-study MG0-MG5 QC plots.",
    "Input: original combine/data/processed_data/mg_cca_meta_v2.qs; 21,244 cells.",
    "Displayed values are existing nCount_RNA and nFeature_RNA metadata; no cells were filtered or changed.",
    "The plots pool cells across studies and are descriptive. Study-specific library depth and chemistry can influence them; consult MG0_MG5_QC_summary_by_Study.csv for study-resolved values.",
    paste0("Seurat version: ", as.character(packageVersion("Seurat")))
  ),
  file.path(output_dir, "README.txt")
)

message("Completed MG0-MG5 QC plots: ", output_dir)
