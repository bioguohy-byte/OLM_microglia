#!/usr/bin/env Rscript

# QC violin plots for the six original full-cell study objects.
# Historical MG cells are identified by exact barcode matching to
# combine/data/processed_data/mg_cca_meta_v2.qs.  Their CellType labels are
# replaced for plotting only: MG0 is displayed as OLM and MG1-MG5 are retained.

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
out_dir <- file.path(project_root, "combine/results/qc_by_source_celltype_with_MG0_MG5_20260920")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

final_mg_path <- file.path(project_root, "combine/data/processed_data/mg_cca_meta_v2.qs")
study_manifest <- data.frame(
  Study = c("Su", "Chia", "Fullard", "Jakel", "Grubman", "Yao"),
  historical_Study = c(
    "Su et al., 2022",
    "Chia et al., 2025",
    "Fullard et al., 2021",
    "Jakel et al., 2019",
    "Grubman et al., 2019",
    "Yao et al., 2024"
  ),
  source_path = c(
    file.path(project_root, "Su_et_al./data/processed_data/Su_CSC_alldata_v2.qs"),
    file.path(project_root, "Chia_et_al./data/processed_data/Chia_v2.qs"),
    file.path(project_root, "Fullard_et_al./data/processed_data/Fullard_v2.qs"),
    file.path(project_root, "Jakel_et_al./data/processed_data/Jakel_v2.qs"),
    file.path(project_root, "Grubman_et_al./data/processed_data/Grubman_v2.qs"),
    file.path(project_root, "Yao_et_al./data/processed_data/Yao_v2.qs")
  ),
  stringsAsFactors = FALSE
)

required_meta <- c("CellType", "nCount_RNA", "nFeature_RNA")
mg_levels <- c("OLM", "MG1", "MG2", "MG3", "MG4", "MG5")
mg_palette <- c(
  OLM = "#4DAF4A",
  MG1 = "#A65628",
  MG2 = "#377EB8",
  MG3 = "#984EA3",
  MG4 = "#E6AB02",
  MG5 = "#8B0000"
)

final_mg <- qread(final_mg_path)
final_meta <- final_mg[[]]
if (!all(c("Study", "CellType1") %in% colnames(final_meta))) {
  stop("The historical MG object is missing Study and/or CellType1 metadata.")
}
if (!all(final_meta$CellType1 %in% paste0("MG", 0:5))) {
  stop("Historical MG object contains CellType1 values outside MG0-MG5.")
}

make_violin <- function(plot_data, metric, y_label, palette, n_categories) {
  ggplot(plot_data, aes(x = plot_CellType, y = .data[[metric]], fill = plot_CellType)) +
    geom_violin(scale = "width", trim = FALSE, colour = "#333333", linewidth = 0.22) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", colour = "#333333", linewidth = 0.22) +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.5, fill = "black", colour = "black") +
    scale_fill_manual(values = palette, drop = FALSE) +
    scale_y_continuous(labels = label_number(big.mark = ",")) +
    labs(x = NULL, y = y_label) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "none",
      plot.margin = margin(6, 10, 10, 6)
    )
}

summary_list <- list()
mapping_audit_list <- list()
pdf_path <- file.path(out_dir, "QC_violin_by_source_CellType_with_historical_MG0_MG5.pdf")
grDevices::pdf(pdf_path, width = 14, height = 12, onefile = TRUE, useDingbats = FALSE)

for (i in seq_len(nrow(study_manifest))) {
  study <- study_manifest$Study[[i]]
  historical_study <- study_manifest$historical_Study[[i]]
  source_path <- study_manifest$source_path[[i]]
  if (!file.exists(source_path)) {
    stop(sprintf("Source object does not exist for %s: %s", study, source_path))
  }

  source_obj <- qread(source_path)
  source_meta <- source_obj[[]]
  missing_cols <- setdiff(required_meta, colnames(source_meta))
  if (length(missing_cols) > 0L) {
    stop(sprintf("%s is missing required metadata: %s", study, paste(missing_cols, collapse = ", ")))
  }

  final_idx <- which(final_meta$Study == historical_study)
  final_cells <- rownames(final_meta)[final_idx]
  if (length(final_cells) == 0L) {
    stop(sprintf("No historical MG cells were found for %s in mg_cca_meta_v2.qs.", study))
  }
  if (!all(final_cells %in% rownames(source_meta))) {
    absent <- setdiff(final_cells, rownames(source_meta))
    stop(sprintf("%s historical MG cell barcodes are absent from the full-cell object for %s (first: %s).",
      length(absent), study, absent[[1]]))
  }

  source_meta$plot_CellType <- as.character(source_meta$CellType)
  source_meta$plot_CellType[is.na(source_meta$plot_CellType) | source_meta$plot_CellType == ""] <- "Unannotated"
  mg_state <- rep(NA_character_, nrow(source_meta))
  names(mg_state) <- rownames(source_meta)
  mg_state[final_cells] <- as.character(final_meta[final_cells, "CellType1", drop = TRUE])
  source_meta$plot_CellType[!is.na(mg_state)] <- ifelse(mg_state[!is.na(mg_state)] == "MG0", "OLM", mg_state[!is.na(mg_state)])

  non_mg_levels <- sort(setdiff(unique(source_meta$plot_CellType), mg_levels))
  plot_levels <- c(mg_levels[mg_levels %in% source_meta$plot_CellType], non_mg_levels)
  source_meta$plot_CellType <- factor(source_meta$plot_CellType, levels = plot_levels)

  other_levels <- setdiff(plot_levels, names(mg_palette))
  other_palette <- setNames(hue_pal(l = 60, c = 100)(length(other_levels)), other_levels)
  plot_palette <- c(mg_palette, other_palette)[plot_levels]

  plot_data <- data.frame(
    plot_CellType = source_meta$plot_CellType,
    nCount_RNA = as.numeric(source_meta$nCount_RNA),
    nFeature_RNA = as.numeric(source_meta$nFeature_RNA),
    stringsAsFactors = FALSE
  )
  if (anyNA(plot_data$nCount_RNA) || anyNA(plot_data$nFeature_RNA)) {
    stop(sprintf("%s contains missing nCount_RNA or nFeature_RNA values.", study))
  }

  group_sizes <- table(plot_data$plot_CellType)
  p_count <- make_violin(plot_data, "nCount_RNA", "nCount_RNA (total RNA counts)", plot_palette, length(plot_levels))
  p_feature <- make_violin(plot_data, "nFeature_RNA", "nFeature_RNA (detected genes)", plot_palette, length(plot_levels))
  single_plot_theme <- theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 10)
  )
  p_count_single <- p_count +
    labs(
      title = sprintf("%s: nCount_RNA by source CellType", study),
      subtitle = sprintf("Historical microglia are subdivided as OLM (MG0) and MG1-MG5; total cells = %s", comma(nrow(plot_data)))
    ) +
    single_plot_theme
  p_feature_single <- p_feature +
    labs(
      title = sprintf("%s: nFeature_RNA by source CellType", study),
      subtitle = sprintf("Historical microglia are subdivided as OLM (MG0) and MG1-MG5; total cells = %s", comma(nrow(plot_data)))
    ) +
    single_plot_theme
  combined_plot <- p_count / p_feature +
    plot_annotation(
      title = sprintf("%s: QC distributions by source CellType", study),
      subtitle = sprintf("Historical microglia are subdivided as OLM (MG0) and MG1-MG5; total cells = %s", comma(nrow(plot_data))),
      theme = theme(plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 10))
    )
  print(combined_plot)
  ggsave(
    filename = file.path(out_dir, sprintf("%s_QC_violin_by_source_CellType_with_historical_MG0_MG5.png", study)),
    plot = combined_plot, width = 14, height = 12, units = "in", dpi = 220, bg = "white"
  )
  single_width <- max(12, length(plot_levels) * 0.65)
  ggsave(
    filename = file.path(out_dir, sprintf("%s_nCount_RNA_by_source_CellType_with_historical_MG0_MG5.pdf", study)),
    plot = p_count_single, width = single_width, height = 7.5, units = "in", bg = "white"
  )
  ggsave(
    filename = file.path(out_dir, sprintf("%s_nFeature_RNA_by_source_CellType_with_historical_MG0_MG5.pdf", study)),
    plot = p_feature_single, width = single_width, height = 7.5, units = "in", bg = "white"
  )
  ggsave(
    filename = file.path(out_dir, sprintf("%s_nCount_RNA_by_source_CellType_with_historical_MG0_MG5.png", study)),
    plot = p_count_single, width = single_width, height = 7.5, units = "in", dpi = 220, bg = "white"
  )
  ggsave(
    filename = file.path(out_dir, sprintf("%s_nFeature_RNA_by_source_CellType_with_historical_MG0_MG5.png", study)),
    plot = p_feature_single, width = single_width, height = 7.5, units = "in", dpi = 220, bg = "white"
  )

  group_split <- split(plot_data, plot_data$plot_CellType, drop = TRUE)
  study_summary <- do.call(rbind, lapply(names(group_split), function(group_name) {
    group_data <- group_split[[group_name]]
    data.frame(
      Study = study,
      plot_CellType = group_name,
      n_cells = nrow(group_data),
      nCount_RNA_median = median(group_data$nCount_RNA),
      nCount_RNA_mean = mean(group_data$nCount_RNA),
      nFeature_RNA_median = median(group_data$nFeature_RNA),
      nFeature_RNA_mean = mean(group_data$nFeature_RNA),
      stringsAsFactors = FALSE
    )
  }))
  summary_list[[study]] <- study_summary

  mapping_audit_list[[study]] <- data.frame(
    Study = study,
    historical_CellType1 = paste0("MG", 0:5),
    plotted_label = mg_levels,
    n_cells = as.integer(table(factor(mg_state[!is.na(mg_state)], levels = paste0("MG", 0:5)))),
    stringsAsFactors = FALSE
  )
}
grDevices::dev.off()

summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, file.path(out_dir, "QC_summary_by_study_and_plot_CellType.csv"), row.names = FALSE)
mapping_audit_df <- do.call(rbind, mapping_audit_list)
write.csv(mapping_audit_df, file.path(out_dir, "historical_MG0_MG5_barcode_mapping_audit.csv"), row.names = FALSE)

writeLines(c(
  "QC violin plot provenance",
  "",
  "Input metrics: nCount_RNA and nFeature_RNA stored in each original full-cell Seurat object.",
  "Non-microglial labels: original CellType values from each full-cell object.",
  "Microglial labels: exact cell barcodes in combine/data/processed_data/mg_cca_meta_v2.qs.",
  "Display mapping: MG0 -> OLM; MG1-MG5 retain their historical labels.",
  "No source Seurat object is modified or saved by this script."
), file.path(out_dir, "README.txt"))

message(sprintf("Wrote PDF: %s", pdf_path))
message(sprintf("Wrote %d study-specific PNG files and two CSV audit tables to %s", nrow(study_manifest), out_dir))
