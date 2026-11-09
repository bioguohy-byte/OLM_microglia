# QC comparison across MG0--MG5 -----------------------------------------------
# Input: integrated, annotated microglia object on server 24.
# Output: one PNG and one PDF per QC metric, a combined figure in both formats,
#         and tabular descriptive/global-test results for verification.

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine")
input_file <- file.path(project_dir, "data/processed_data/mg_cca_meta_v2.qs")
figure_dir <- file.path(project_dir, "results/figure/MG0_MG5_QC")
table_dir <- file.path(project_dir, "results/table/MG0_MG5_QC")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Input object not found: ", input_file)

# The metric calculations all use the same RNA count layer. This avoids mixing
# source-specific QC metadata that may have been calculated by different rules.
mg_obj <- qread(input_file, nthreads = 16)
if (!"CellType1" %in% colnames(mg_obj@meta.data)) {
  stop("CellType1 is absent from the input object.")
}

rna_counts <- GetAssayData(mg_obj, assay = "RNA", layer = "counts")
if (!identical(colnames(rna_counts), rownames(mg_obj@meta.data))) {
  stop("RNA count columns and metadata rows are not aligned.")
}

mg_levels <- paste0("MG", 0:5)
cell_state <- factor(mg_obj$CellType1, levels = mg_levels)
if (anyNA(cell_state)) stop("Unexpected CellType1 labels detected.")

total_umi <- Matrix::colSums(rna_counts)
detected_genes <- Matrix::colSums(rna_counts != 0)

percentage_from_gene_pattern <- function(count_matrix, total_counts, pattern) {
  matched_genes <- grep(pattern, rownames(count_matrix), value = TRUE)
  if (!length(matched_genes)) stop("No genes matched pattern: ", pattern)
  100 * Matrix::colSums(count_matrix[matched_genes, , drop = FALSE]) / total_counts
}

# MT- captures human mitochondrial genes. RPL/RPS captures cytosolic ribosomal
# genes; mitochondrial ribosomal genes (MRPL/MRPS) are deliberately excluded.
percent_mt <- percentage_from_gene_pattern(rna_counts, total_umi, "^MT-")
percent_ribo <- percentage_from_gene_pattern(rna_counts, total_umi, "^RP[SL]")

qc_data <- data.frame(
  CellType1 = cell_state,
  nCount_RNA = as.numeric(total_umi),
  nFeature_RNA = as.numeric(detected_genes),
  percent_mt = as.numeric(percent_mt),
  percent_ribo = as.numeric(percent_ribo),
  check.names = FALSE
)

metric_specs <- data.frame(
  metric = c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_ribo"),
  label = c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "percent.ribo"
  ),
  log_scale = c(TRUE, TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

qc_long <- do.call(rbind, lapply(metric_specs$metric, function(metric_name) {
  data.frame(
    CellType1 = qc_data$CellType1,
    metric = metric_name,
    value = qc_data[[metric_name]],
    stringsAsFactors = FALSE
  )
}))
qc_long$metric <- factor(qc_long$metric, levels = metric_specs$metric)

qc_summary <- qc_long %>%
  group_by(metric, CellType1) %>%
  summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    q1 = unname(quantile(value, 0.25)),
    q3 = unname(quantile(value, 0.75)),
    .groups = "drop"
  ) %>%
  left_join(metric_specs[, c("metric", "label")], by = "metric")

qc_kruskal <- do.call(rbind, lapply(metric_specs$metric, function(metric_name) {
  test_result <- kruskal.test(qc_data[[metric_name]] ~ qc_data$CellType1)
  group_number <- nlevels(qc_data$CellType1)
  cell_number <- nrow(qc_data)
  epsilon_squared <- max(0, (unname(test_result$statistic) - group_number + 1) /
                            (cell_number - group_number))
  data.frame(
    metric = metric_name,
    statistic = unname(test_result$statistic),
    df = unname(test_result$parameter),
    p_value = test_result$p.value,
    epsilon_squared = epsilon_squared,
    stringsAsFactors = FALSE
  )
})) %>%
  left_join(metric_specs[, c("metric", "label")], by = "metric")

write.csv(qc_summary,
          file.path(table_dir, "MG0_MG5_QC_summary.csv"),
          row.names = FALSE)
write.csv(qc_kruskal,
          file.path(table_dir, "MG0_MG5_QC_kruskal_wallis.csv"),
          row.names = FALSE)

# The final integrated metadata intentionally retains no comparable continuous
# doublet score. Record this explicitly instead of mixing binary source-specific
# DoubletFinder calls with a numeric score.
doublet_candidates <- grep(
  "doublet|pANN|scDblFinder|doublet_score",
  colnames(mg_obj@meta.data),
  value = TRUE,
  ignore.case = TRUE
)
doublet_status <- data.frame(
  status = if (length(doublet_candidates)) "candidate columns found" else "no comparable doublet score available",
  candidate_columns = if (length(doublet_candidates)) paste(doublet_candidates, collapse = "; ") else NA_character_,
  next_step = if (length(doublet_candidates)) {
    "Validate whether candidate columns are continuous and available for every study before plotting."
  } else {
    "Compute one harmonized doublet score per Study/Sample from RNA counts before adding this QC panel."
  },
  stringsAsFactors = FALSE
)
write.csv(doublet_status,
          file.path(table_dir, "MG0_MG5_doublet_score_status.csv"),
          row.names = FALSE)

mg_colors <- c(
  MG0 = "#4DAF4A", MG1 = "#A65628", MG2 = "#377EB8",
  MG3 = "#984EA3", MG4 = "#E6AB02", MG5 = "#8B0000"
)

make_qc_plot <- function(metric_name, metric_label, use_log_scale) {
  metric_test <- qc_kruskal[qc_kruskal$metric == metric_name, , drop = FALSE]
  caption_text <- sprintf(
    "Kruskal-Wallis: H = %.1f, p = %.3g, epsilon-squared = %.3f",
    metric_test$statistic, metric_test$p_value, metric_test$epsilon_squared
  )

  plot_object <- ggplot(qc_data, aes(x = CellType1, y = .data[[metric_name]], fill = CellType1)) +
    geom_violin(trim = TRUE, scale = "width", linewidth = 0.25, color = "grey25") +
    geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.3, fill = "white") +
    scale_fill_manual(values = mg_colors, drop = FALSE) +
    labs(
      title = paste0("MG0-MG5 QC: ", metric_label),
      x = NULL,
      y = metric_label,
      caption = caption_text
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0)
    )

  if (use_log_scale) {
    plot_object <- plot_object + scale_y_continuous(trans = "log10")
  }
  plot_object
}

qc_plots <- list()
for (row_index in seq_len(nrow(metric_specs))) {
  metric_name <- metric_specs$metric[row_index]
  plot_object <- make_qc_plot(
    metric_name = metric_name,
    metric_label = metric_specs$label[row_index],
    use_log_scale = metric_specs$log_scale[row_index]
  )
  qc_plots[[metric_name]] <- plot_object

  figure_base <- file.path(figure_dir, paste0("MG0_MG5_QC_", metric_name))
  ggsave(paste0(figure_base, ".png"), plot_object,
         width = 7, height = 5, dpi = 300, bg = "white")
  ggsave(paste0(figure_base, ".pdf"), plot_object,
         width = 7, height = 5, bg = "white")
}

combined_plot <- wrap_plots(qc_plots, ncol = 2) +
  plot_annotation(title = "QC comparison across MG0-MG5")
combined_base <- file.path(figure_dir, "MG0_MG5_QC_combined")
ggsave(paste0(combined_base, ".png"), combined_plot,
       width = 14, height = 10, dpi = 300, bg = "white")
ggsave(paste0(combined_base, ".pdf"), combined_plot,
       width = 14, height = 10, bg = "white")

message("QC comparison completed.")
message("Figures: ", figure_dir)
message("Tables: ", table_dir)
message("Doublet score status: ", doublet_status$status)
