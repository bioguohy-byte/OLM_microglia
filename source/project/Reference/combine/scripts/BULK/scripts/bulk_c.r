# 1 Parameters ----

options(stringsAsFactors = FALSE)
rm(list = ls())

check_and_load_packages <- function(pkgs, bioc_pkgs = character()) {
  missing_cran <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing_cran) > 0 || length(missing_bioc) > 0) {
    msg <- c()
    
    if (length(missing_cran) > 0) {
      msg <- c(msg, paste0("Missing CRAN packages: ", paste(missing_cran, collapse = ", ")))
    }
    
    if (length(missing_bioc) > 0) {
      msg <- c(msg, paste0("Missing Bioconductor packages: ", paste(missing_bioc, collapse = ", ")))
    }
    
    stop(
      paste(
        c(
          msg,
          "Please install required packages before running this script.",
          "Example:",
          "install.packages(c('ggplot2', 'pheatmap', 'RColorBrewer'))",
          "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')",
          "BiocManager::install(c('DESeq2'))"
        ),
        collapse = "\n"
      )
    )
  }
  
  invisible(lapply(c(pkgs, bioc_pkgs), function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

check_and_load_packages(
  pkgs = c("ggplot2", "pheatmap", "RColorBrewer", "grid"),
  bioc_pkgs = c("DESeq2")
)

params <- list(
  input_file = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/processed_data/All_gene_counts.list"),
  output_dir = "bulk_analysis_output",
  sep = "\t",
  gene_column = 1,
  sample_columns = NULL,
  force_integer_counts = TRUE,
  aggregate_duplicated_genes = TRUE,
  min_count = 10,
  min_samples = 2,
  design_formula = ~ condition,
  contrast = c("condition", "treat", "con"),
  alpha = 0.05,
  top_variable_genes = 50,
  p_adj_cutoff = 0.05,
  log2fc_cutoff = 1
)

sample_info <- data.frame(
  sample = c("con1", "con2", "con3", "con4", "treat1", "treat2", "treat3", "treat4"),
  condition = c("con", "con", "con", "con", "treat", "treat", "treat", "treat"),
  row.names = c("con1", "con2", "con3", "con4", "treat1", "treat2", "treat3", "treat4"),
  stringsAsFactors = FALSE
)

sample_info$condition <- factor(sample_info$condition, levels = c("con", "treat"))

dir.create(params$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(params$output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(params$output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(params$output_dir, "rdata"), showWarnings = FALSE, recursive = TRUE)

message("Parameters loaded successfully.")


# 2 Process Data --------------

read_count_matrix <- function(file, sep = "\t", gene_column = 1, sample_columns = NULL) {
  raw_df <- read.delim(
    file = file,
    sep = sep,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  
  if (ncol(raw_df) < 2) {
    stop("Input file must contain at least one gene column and one sample column.")
  }
  
  if (is.numeric(gene_column)) {
    if (gene_column < 1 || gene_column > ncol(raw_df)) {
      stop("gene_column is out of range.")
    }
    gene_col_name <- colnames(raw_df)[gene_column]
  } else {
    if (!gene_column %in% colnames(raw_df)) {
      stop("Specified gene_column not found in input file.")
    }
    gene_col_name <- gene_column
  }
  
  if (is.null(sample_columns)) {
    sample_columns <- setdiff(colnames(raw_df), gene_col_name)
  } else {
    if (!all(sample_columns %in% colnames(raw_df))) {
      stop("Some sample_columns are not present in the input file.")
    }
  }
  
  gene_ids <- raw_df[[gene_col_name]]
  if (any(is.na(gene_ids)) || any(gene_ids == "")) {
    stop("Gene ID column contains NA or empty values.")
  }
  
  count_df <- raw_df[, sample_columns, drop = FALSE]
  count_df[] <- lapply(count_df, function(x) suppressWarnings(as.numeric(x)))
  
  if (anyNA(count_df)) {
    stop("Count matrix contains NA after numeric conversion. Please check the input file format.")
  }
  
  count_mat <- as.matrix(count_df)
  rownames(count_mat) <- gene_ids
  
  count_mat
}

validate_and_clean_counts <- function(count_mat, force_integer_counts = TRUE, aggregate_duplicated_genes = TRUE) {
  if (!is.matrix(count_mat)) {
    stop("count_mat must be a matrix.")
  }
  
  if (!all(is.finite(count_mat))) {
    stop("Count matrix contains non-finite values.")
  }
  
  if (any(count_mat < 0)) {
    stop("Count matrix contains negative values.")
  }
  
  if (aggregate_duplicated_genes && anyDuplicated(rownames(count_mat))) {
    message("Duplicated gene IDs detected. Aggregating duplicated rows by sum.")
    count_mat <- rowsum(count_mat, group = rownames(count_mat), reorder = FALSE)
  }
  
  if (force_integer_counts) {
    if (any(abs(count_mat - round(count_mat)) > 1e-8)) {
      message("Non-integer values detected. Rounding to nearest integer.")
    }
    count_mat <- round(count_mat)
    storage.mode(count_mat) <- "integer"
  }
  
  count_mat
}

filter_low_expression <- function(count_mat, min_count = 10, min_samples = 2) {
  keep <- rowSums(count_mat >= min_count) >= min_samples
  list(
    filtered_matrix = count_mat[keep, , drop = FALSE],
    keep = keep
  )
}

count_mat <- read_count_matrix(
  file = params$input_file,
  sep = params$sep,
  gene_column = params$gene_column,
  sample_columns = params$sample_columns
)

count_mat <- validate_and_clean_counts(
  count_mat = count_mat,
  force_integer_counts = params$force_integer_counts,
  aggregate_duplicated_genes = params$aggregate_duplicated_genes
)

if (!all(sample_info$sample %in% colnames(count_mat))) {
  stop("Some samples in sample_info are not present in the count matrix.")
}

count_mat <- count_mat[, sample_info$sample, drop = FALSE]
sample_info <- sample_info[colnames(count_mat), , drop = FALSE]

if (!identical(rownames(sample_info), colnames(count_mat))) {
  stop("Sample names/order mismatch between sample_info and count matrix.")
}

filtered_res <- filter_low_expression(
  count_mat = count_mat,
  min_count = params$min_count,
  min_samples = params$min_samples
)

count_mat_filtered <- filtered_res$filtered_matrix

if (nrow(count_mat_filtered) == 0) {
  stop("No genes passed filtering. Please relax min_count/min_samples.")
}

message("Original genes: ", nrow(count_mat))
message("Filtered genes: ", nrow(count_mat_filtered))

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_mat_filtered,
  colData = sample_info,
  design = params$design_formula
)

dds <- DESeq2::DESeq(dds)

res <- DESeq2::results(
  dds,
  contrast = params$contrast,
  alpha = params$alpha
)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]
res_df <- res_df[order(res_df$padj, res_df$pvalue), , drop = FALSE]

norm_counts <- DESeq2::counts(dds, normalized = TRUE)
vst_obj <- DESeq2::vst(dds, blind = FALSE)
vst_mat <- SummarizedExperiment::assay(vst_obj)

res_df$regulation <- "NotSig"
res_df$regulation[
  !is.na(res_df$padj) &
    !is.na(res_df$log2FoldChange) &
    res_df$padj < params$p_adj_cutoff &
    res_df$log2FoldChange >= params$log2fc_cutoff
] <- "Up"

res_df$regulation[
  !is.na(res_df$padj) &
    !is.na(res_df$log2FoldChange) &
    res_df$padj < params$p_adj_cutoff &
    res_df$log2FoldChange <= -params$log2fc_cutoff
] <- "Down"

deg_up <- res_df[res_df$regulation == "Up", , drop = FALSE]
deg_down <- res_df[res_df$regulation == "Down", , drop = FALSE]
deg_all <- res_df[res_df$regulation %in% c("Up", "Down"), , drop = FALSE]

summary_stats <- data.frame(
  metric = c("input_genes", "filtered_genes", "significant_up", "significant_down", "significant_total"),
  value = c(nrow(count_mat), nrow(count_mat_filtered), nrow(deg_up), nrow(deg_down), nrow(deg_all)),
  stringsAsFactors = FALSE
)

message("Data processing finished.")


# 3 Generate Plots ----

plot_dir <- file.path(params$output_dir, "plots")

libsize_df <- data.frame(
  sample = colnames(count_mat),
  library_size = colSums(count_mat),
  condition = sample_info$condition,
  stringsAsFactors = FALSE
)

p_libsize <- ggplot(libsize_df, aes(x = sample, y = library_size, fill = condition)) +
  geom_col(width = 0.7) +
  theme_bw(base_size = 12) +
  labs(title = "Library Size", x = "Sample", y = "Raw Read Counts") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(plot_dir, "01_library_size_barplot.pdf"),
  plot = p_libsize,
  width = 8,
  height = 5
)

pca_data <- DESeq2::plotPCA(vst_obj, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3.5) +
  theme_bw(base_size = 12) +
  labs(
    title = "PCA of Samples",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(plot_dir, "02_PCA_plot.pdf"),
  plot = p_pca,
  width = 7,
  height = 6
)

sample_dists <- dist(t(SummarizedExperiment::assay(vst_obj)))
sample_dist_mat <- as.matrix(sample_dists)

pdf(file.path(plot_dir, "03_sample_distance_heatmap.pdf"), width = 7, height = 6)
pheatmap::pheatmap(
  sample_dist_mat,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255),
  main = "Sample Distance Heatmap"
)
dev.off()

gene_var <- apply(vst_mat, 1, var, na.rm = TRUE)
valid_var <- !is.na(gene_var) & gene_var > 0
vst_mat_heatmap <- vst_mat[valid_var, , drop = FALSE]
gene_var_heatmap <- gene_var[valid_var]

if (nrow(vst_mat_heatmap) == 0) {
  warning("No genes with variance > 0 were found. Skip top variable genes heatmap.")
} else {
  top_n <- min(params$top_variable_genes, nrow(vst_mat_heatmap))
  top_genes <- names(sort(gene_var_heatmap, decreasing = TRUE))[seq_len(top_n)]
  
  annotation_col <- data.frame(condition = sample_info$condition)
  rownames(annotation_col) <- rownames(sample_info)
  
  ph <- pheatmap::pheatmap(
    vst_mat_heatmap[top_genes, , drop = FALSE],
    scale = "row",
    annotation_col = annotation_col,
    show_rownames = TRUE,
    fontsize_row = 7,
    main = paste0("Top ", top_n, " Variable Genes"),
    silent = TRUE
  )
  
  pdf(file.path(plot_dir, "04_top_variable_genes_heatmap.pdf"), width = 8, height = 10)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()
}

pdf(file.path(plot_dir, "05_MA_plot.pdf"), width = 7, height = 6)
DESeq2::plotMA(res, alpha = params$alpha, main = "MA Plot")
dev.off()

volcano_df <- res_df
volcano_df$neg_log10_padj <- -log10(volcano_df$padj)
volcano_df$neg_log10_padj[is.infinite(volcano_df$neg_log10_padj)] <- NA

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj, color = regulation)) +
  geom_point(alpha = 0.7, size = 1.8, na.rm = TRUE) +
  scale_color_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "NotSig" = "grey70")) +
  geom_vline(xintercept = c(-params$log2fc_cutoff, params$log2fc_cutoff), linetype = 2, color = "black") +
  geom_hline(yintercept = -log10(params$p_adj_cutoff), linetype = 2, color = "black") +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Volcano Plot: ", params$contrast[2], " vs ", params$contrast[3]),
    x = "log2 Fold Change",
    y = "-log10 adjusted P value"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(plot_dir, "06_volcano_plot.pdf"),
  plot = p_volcano,
  width = 7,
  height = 6
)

message("Plots generated successfully.")


# 4 Save Results ----

table_dir <- file.path(params$output_dir, "tables")
rdata_dir <- file.path(params$output_dir, "rdata")

write.csv(
  sample_info,
  file = file.path(table_dir, "sample_info.csv"),
  row.names = TRUE,
  quote = FALSE
)

write.csv(
  data.frame(gene_id = rownames(count_mat), count_mat, check.names = FALSE),
  file = file.path(table_dir, "raw_counts.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(gene_id = rownames(count_mat_filtered), count_mat_filtered, check.names = FALSE),
  file = file.path(table_dir, "filtered_counts.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(gene_id = rownames(norm_counts), norm_counts, check.names = FALSE),
  file = file.path(table_dir, "normalized_counts_DESeq2.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(gene_id = rownames(vst_mat), vst_mat, check.names = FALSE),
  file = file.path(table_dir, "vst_transformed_counts.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  res_df,
  file = file.path(table_dir, "DEG_all_results.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  deg_up,
  file = file.path(table_dir, "DEG_up.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  deg_down,
  file = file.path(table_dir, "DEG_down.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  deg_all,
  file = file.path(table_dir, "DEG_significant.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  summary_stats,
  file = file.path(table_dir, "summary_statistics.csv"),
  row.names = FALSE,
  quote = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(params$output_dir, "sessionInfo.txt")
)

save(
  params,
  sample_info,
  count_mat,
  count_mat_filtered,
  dds,
  res,
  res_df,
  norm_counts,
  vst_mat,
  file = file.path(rdata_dir, "bulk_analysis_workspace.RData")
)

message("All results saved to: ", normalizePath(params$output_dir))
message("Bulk RNA-seq analysis completed successfully.")
