#!/usr/bin/env Rscript

# Rebuild the current six-study diagnostic figures after the reanalysis manual
# cluster-2 label was renamed from MG0 to MG00.  Doublet calls are read from the
# previously completed three-method diagnostic workflow and are not recomputed.

rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
invisible(gc())

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(qs)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
old_stage3_root <- file.path(project_root, "combine/results/full_count_reanalysis/03_current_MG0_doublet_diagnosis")
rename_root <- file.path(project_root, "combine/results/full_count_reanalysis/02_microglia_integration_MG00")
result_root <- file.path(project_root, "combine/results/full_count_reanalysis/03_current_MG00_doublet_diagnosis")
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
object_dir <- file.path(result_root, "object")
runtime_dir <- file.path(result_root, "runtime")
scrublet_dir <- file.path(old_stage3_root, "scrublet_output")
for (directory in c(figure_dir, table_dir, object_dir, runtime_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

object_file <- file.path(rename_root, "object/current_MG00_QC_before_doublet_removal_three_method_partial.qs")
if (!file.exists(object_file)) stop("MG00 diagnostic object is missing: ", object_file)
if (!dir.exists(scrublet_dir)) stop("Previously completed Scrublet calls are missing: ", scrublet_dir)

study_order <- c(
  "Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
  "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"
)
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width = 7.5, height = 4.5) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
    dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
    device = grDevices::pdf, bg = "white", limitsize = FALSE)
}

diagnostic <- qs::qread(object_file)
meta <- diagnostic[[]]
meta$Cell <- rownames(meta)
required <- c(
  "Study", "diagnostic_cluster", "MG00_status", "MG00_transfer_score",
  "DoubletFinder_call", "scDblFinder_call", "DF_partition"
)
if (length(setdiff(required, colnames(meta)))) {
  stop("MG00 diagnostic metadata lacks: ", paste(setdiff(required, colnames(meta)), collapse = ", "))
}
if (any(c("MG0_status", "MG0_transfer_score") %in% colnames(meta))) {
  stop("An obsolete reanalysis MG0 metadata field remains in the MG00 object.")
}

scrublet_files <- list.files(scrublet_dir, pattern = "^Scrublet_.*\\.csv$", full.names = TRUE)
scrublet_files <- scrublet_files[basename(scrublet_files) != "Scrublet_run_summary.csv"]
if (!length(scrublet_files)) stop("No Scrublet per-cell files are available.")
scrublet <- do.call(rbind, lapply(scrublet_files, read.csv, stringsAsFactors = FALSE, check.names = FALSE))
if (anyDuplicated(scrublet$Cell) || !setequal(scrublet$Cell, meta$Cell)) {
  stop("Scrublet calls do not exactly match the MG00 QC-passed diagnostic cells.")
}
scrublet <- scrublet[match(meta$Cell, scrublet$Cell), , drop = FALSE]
meta$Scrublet_score <- scrublet$Scrublet_score
meta$Scrublet_call <- factor(scrublet$Scrublet_call, levels = c("Singlet", "Doublet"))
meta$Scrublet_threshold <- scrublet$Scrublet_threshold
meta$DoubletFinder_call <- factor(meta$DoubletFinder_call, levels = c("Singlet", "Doublet"))
meta$scDblFinder_call <- factor(meta$scDblFinder_call, levels = c("Singlet", "Doublet"))
meta$MG00_status <- factor(meta$MG00_status,
  levels = c("Other predicted cluster", "MG00 (manual cluster 2 transfer)"))

umap <- Seurat::Embeddings(diagnostic, "umap_diagnostic")
meta$UMAP_1 <- umap[meta$Cell, 1]
meta$UMAP_2 <- umap[meta$Cell, 2]
if (!setequal(study_order, unique(as.character(meta$Study)))) {
  stop("Expected exactly the six selected studies in the MG00 diagnostic object.")
}

umap_base <- function(dat, title, subtitle = NULL) {
  ggplot(dat, aes(UMAP_1, UMAP_2)) +
    coord_equal() +
    labs(title = title, subtitle = subtitle, x = "Diagnostic UMAP 1", y = "Diagnostic UMAP 2") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      plot.margin = margin(t = 12, r = 14, b = 10, l = 14),
      legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 8)
    )
}
make_cluster_plot <- function(dat, study) {
  umap_base(dat, "Diagnostic clusters", study) +
    geom_point(aes(colour = factor(diagnostic_cluster)), size = 0.48, alpha = 0.84) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2)))
}
make_MG00_plot <- function(dat, study) {
  dat <- dat[order(dat$MG00_status), , drop = FALSE]
  umap_base(dat, "MG00", study) +
    geom_point(aes(colour = MG00_status), size = 0.48, alpha = 0.88) +
    scale_colour_manual(values = c("Other predicted cluster" = "#D9D9D9", "MG00 (manual cluster 2 transfer)" = "#C51B7D")) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
}
make_doublet_plot <- function(dat, column, method, study) {
  dat <- dat[order(dat[[column]]), , drop = FALSE]
  umap_base(dat, method, study) +
    geom_point(aes(colour = .data[[column]]), size = 0.48, alpha = 0.88) +
    scale_colour_manual(values = c(Singlet = "#BDBDBD", Doublet = "#D7301F")) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
}

plot_types <- c("diagnostic_cluster", "MG00", "DoubletFinder", "scDblFinder", "Scrublet")
plot_store <- setNames(lapply(plot_types, function(x) list()), plot_types)
for (study in study_order) {
  dat <- meta[meta$Study == study, , drop = FALSE]
  plots <- list(
    diagnostic_cluster = make_cluster_plot(dat, study),
    MG00 = make_MG00_plot(dat, study),
    DoubletFinder = make_doublet_plot(dat, "DoubletFinder_call", "DoubletFinder", study),
    scDblFinder = make_doublet_plot(dat, "scDblFinder_call", "scDblFinder", study),
    Scrublet = make_doublet_plot(dat, "Scrublet_call", "Scrublet", study)
  )
  for (plot_type in names(plots)) {
    plot_store[[plot_type]][[study]] <- plots[[plot_type]] +
      labs(title = study, subtitle = NULL) +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5), legend.position = "right")
    save_plot_both(
      plots[[plot_type]],
      file.path(figure_dir, paste0("UMAP_", safe_name(study), "_", plot_type))
    )
  }
}
for (plot_type in plot_types) {
  composite <- wrap_plots(plot_store[[plot_type]][study_order], ncol = 3) +
    plot_annotation(title = paste0("Six-study diagnostic UMAPs: ", plot_type))
  save_plot_both(
    composite,
    file.path(figure_dir, paste0("Composite_UMAP_", plot_type, "_six_studies")),
    15, 9
  )
  print(composite)
}

long_calls <- rbind(
  data.frame(Cell = meta$Cell, Study = meta$Study, MG00_status = meta$MG00_status, method = "DoubletFinder", call = meta$DoubletFinder_call),
  data.frame(Cell = meta$Cell, Study = meta$Study, MG00_status = meta$MG00_status, method = "scDblFinder", call = meta$scDblFinder_call),
  data.frame(Cell = meta$Cell, Study = meta$Study, MG00_status = meta$MG00_status, method = "Scrublet", call = meta$Scrublet_call)
)
long_calls$method <- factor(long_calls$method, levels = c("DoubletFinder", "scDblFinder", "Scrublet"))
long_calls$call <- factor(long_calls$call, levels = c("Singlet", "Doublet"))
proportion_table <- aggregate(as.integer(long_calls$call == "Doublet"),
  by = list(Study = long_calls$Study, MG00_status = long_calls$MG00_status, method = long_calls$method), FUN = sum)
names(proportion_table)[4] <- "n_doublet"
denominator_table <- aggregate(Cell ~ Study + MG00_status + method, data = long_calls, FUN = length)
names(denominator_table)[4] <- "n_cells"
proportion_table <- merge(proportion_table, denominator_table,
  by = c("Study", "MG00_status", "method"), sort = FALSE)
proportion_table$doublet_proportion <- proportion_table$n_doublet / proportion_table$n_cells
pooled_table <- aggregate(cbind(n_doublet, n_cells) ~ MG00_status + method, data = proportion_table, FUN = sum)
pooled_table$doublet_proportion <- pooled_table$n_doublet / pooled_table$n_cells

p_current_bar <- ggplot(pooled_table, aes(MG00_status, doublet_proportion, fill = MG00_status)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%\\n(n=%d)", 100 * doublet_proportion, n_cells)), vjust = -0.18, size = 3) +
  facet_wrap(~method, nrow = 1) +
  scale_fill_manual(values = c("Other predicted cluster" = "#BDBDBD", "MG00 (manual cluster 2 transfer)" = "#C51B7D")) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = expansion(mult = c(0, 0.17))) +
  labs(
    title = "Predicted doublet proportion: MG00 versus other reanalysis clusters",
    subtitle = "MG00 is the user-selected final integrated cluster 2, transferred to QC-passed pre-removal cells.",
    x = NULL, y = "Predicted doublet proportion"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none")
save_plot_both(p_current_bar, file.path(figure_dir, "Doublet_proportion_MG00_vs_other_three_methods"), 15, 5)
print(p_current_bar)

diagnostic$Scrublet_score <- meta$Scrublet_score
diagnostic$Scrublet_call <- meta$Scrublet_call
diagnostic$Scrublet_threshold <- meta$Scrublet_threshold
qs::qsave(diagnostic, file.path(object_dir, "current_MG00_QC_before_doublet_removal_all_three_calls.qs"), preset = "high")
write.csv(meta, file.path(table_dir, "six_study_cell_metadata_with_all_three_doublet_calls_MG00.csv"), row.names = FALSE)
write.csv(proportion_table, file.path(table_dir, "doublet_proportions_by_study_MG00_vs_other.csv"), row.names = FALSE)
write.csv(pooled_table, file.path(table_dir, "doublet_proportions_pooled_MG00_vs_other.csv"), row.names = FALSE)
write.csv(long_calls, file.path(table_dir, "per_cell_long_three_method_calls_MG00.csv"), row.names = FALSE)
agreement <- meta[, c("Cell", "Study", "MG00_status", "DoubletFinder_call", "scDblFinder_call", "Scrublet_call"), drop = FALSE]
agreement$n_methods_calling_doublet <- rowSums(cbind(
  agreement$DoubletFinder_call == "Doublet", agreement$scDblFinder_call == "Doublet", agreement$Scrublet_call == "Doublet"
), na.rm = TRUE)
write.csv(agreement, file.path(table_dir, "per_cell_three_method_agreement_MG00.csv"), row.names = FALSE)

writeLines(c(
  "Input: MG00-labelled QC-passed pre-removal diagnostic object (16,931 cells).",
  "MG00 status: user-selected final integrated cluster 2 transferred to pre-removal cells; it replaces the reanalysis MG0 status only.",
  "DoubletFinder, scDblFinder and Scrublet calls are the existing calls for these exact cells and technical partitions.",
  "No counts, QC flag, UMAP coordinate, cluster or doublet call was changed by the label rename.",
  "Original CellType1 MG0--MG5 labels remain in their separate historical-pipeline outputs."
), file.path(runtime_dir, "MG00_stage3_method_manifest.txt"))
writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "MG00_stage3_R_sessionInfo.txt"))
message("MG00 stage-3 figures and the all-three-calls MG00 object are ready: ", result_root)
