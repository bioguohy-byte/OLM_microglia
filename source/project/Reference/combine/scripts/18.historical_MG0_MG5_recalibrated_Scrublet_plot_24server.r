#!/usr/bin/env Rscript

# Historical MG0--MG5 doublet diagnosis, recalibrated Scrublet plotting branch.
# Reads a new, non-destructive Scrublet call table and preserves historical
# MG0--MG5 labels plus the previously computed DoubletFinder/scDblFinder calls.
# Every single figure and every composite is saved as both PNG and PDF, and
# composites are printed for direct inspection in RStudio's Plots history.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
source_root <- file.path(project_root, "combine/results/MG0_MG5_doublet_diagnosis")
result_root <- file.path(project_root, "combine/results/MG0_MG5_doublet_diagnosis_recalibrated_scrublet")
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
runtime_dir <- file.path(result_root, "runtime")
for (directory in c(figure_dir, table_dir, runtime_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)

study_order <- c(
  "Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
  "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"
)
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = grDevices::pdf, bg = "white", limitsize = FALSE)
}

source_file <- file.path(source_root, "table/stage1_cell_metadata_before_scrublet.csv")
scrublet_file <- file.path(table_dir, "recalibrated_Scrublet_per_cell.csv")
parameter_file <- file.path(table_dir, "recalibrated_Scrublet_parameters_by_partition.csv")
if (!file.exists(source_file) || !file.exists(scrublet_file) || !file.exists(parameter_file)) {
  stop("Missing source metadata or recalibrated Scrublet outputs. Run script 17 first.")
}

historical <- read.csv(source_file, check.names = FALSE, stringsAsFactors = FALSE)
scrublet <- read.csv(scrublet_file, check.names = FALSE, stringsAsFactors = FALSE)
parameters <- read.csv(parameter_file, check.names = FALSE, stringsAsFactors = FALSE)
required_historical <- c("Cell", "Study", "Sample", "CellType1", "DoubletFinder_call", "scDblFinder_call", "UMAP_1", "UMAP_2")
required_scrublet <- c("Cell", "Study", "Sample", "Scrublet_partition", "Scrublet_score", "Scrublet_recalibrated_call", "Scrublet_recalibrated_threshold", "expected_doublet_rate")
if (length(setdiff(required_historical, names(historical)))) stop("Historical table lacks: ", paste(setdiff(required_historical, names(historical)), collapse = ", "))
if (length(setdiff(required_scrublet, names(scrublet)))) stop("Scrublet table lacks: ", paste(setdiff(required_scrublet, names(scrublet)), collapse = ", "))
if (anyDuplicated(historical$Cell) || anyDuplicated(scrublet$Cell)) stop("Duplicate cells found in input tables.")

idx <- match(historical$Cell, scrublet$Cell)
if (anyNA(idx)) stop("Recalibrated Scrublet calls do not cover all historical cells.")
if (any(historical$Study != scrublet$Study[idx])) stop("Study mapping differs between historical metadata and recalibrated Scrublet calls.")
# The historical stage-1 table has its own study-level `expected_doublet_rate`.
# Keep it unchanged for provenance and retain the new partition-aware quantities
# under Scrublet-specific names; this avoids duplicate names after cbind().
scrublet_for_join <- scrublet[idx, setdiff(names(scrublet), c("Cell", "Study", "Sample")), drop = FALSE]
names(scrublet_for_join)[names(scrublet_for_join) == "expected_doublet_rate"] <- "Scrublet_expected_doublet_rate"
names(scrublet_for_join)[names(scrublet_for_join) == "expected_doublet_n"] <- "Scrublet_expected_doublet_n"
dat <- cbind(historical, scrublet_for_join)
dat$Study <- factor(dat$Study, levels = study_order)
dat$CellType1 <- factor(dat$CellType1, levels = paste0("MG", 0:5))
dat$MG0_status <- factor(ifelse(dat$CellType1 == "MG0", "MG0", "Other MG state"), levels = c("Other MG state", "MG0"))
dat$DoubletFinder_call <- factor(dat$DoubletFinder_call, levels = c("Singlet", "Doublet"))
dat$scDblFinder_call <- factor(dat$scDblFinder_call, levels = c("Singlet", "Doublet"))
dat$Scrublet_recalibrated_call <- factor(dat$Scrublet_recalibrated_call, levels = c("Singlet", "Doublet"))

# Audit statement for the *historical source-pipeline* processing. This does not
# describe this project's three-method diagnostic calculation.
historical_doublet_status <- data.frame(
  Study = study_order,
  historical_doublet_status = c(
    "Historical source pipeline: no doublet-calling/removal step documented",
    "Historical source pipeline: no doublet-calling/removal step documented",
    "Historical source pipeline: no doublet-calling/removal step documented",
    "Historical source pipeline: DoubletFinder result and Singlet filter present; saved call block is commented",
    "Historical source pipeline: no doublet-calling/removal step documented",
    "Historical source pipeline: DoubletFinder per orig.ident; Singlets retained"
  ),
  stringsAsFactors = FALSE
)
dat <- merge(dat, historical_doublet_status, by = "Study", all.x = TRUE, sort = FALSE)
dat$Study <- factor(dat$Study, levels = study_order)
write.csv(historical_doublet_status, file.path(table_dir, "historical_source_pipeline_doublet_status.csv"), row.names = FALSE)
write.csv(dat, file.path(table_dir, "six_study_cell_metadata_with_recalibrated_Scrublet.csv"), row.names = FALSE)

umap_base <- function(subdat, title, subtitle = NULL) {
  ggplot(subdat, aes(UMAP_1, UMAP_2)) +
    coord_equal() +
    labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 7.5, hjust = 0.5),
      legend.title = element_blank()
    )
}
make_mg0_plot <- function(subdat, title) {
  subdat <- subdat[order(subdat$MG0_status), , drop = FALSE]
  umap_base(subdat, title) +
    geom_point(aes(colour = MG0_status), size = 0.45, alpha = 0.85) +
    scale_colour_manual(values = c("Other MG state" = "#d9d9d9", MG0 = "#c51b7d"))
}
make_doublet_plot <- function(subdat, column, title) {
  subdat <- subdat[order(subdat[[column]]), , drop = FALSE]
  umap_base(subdat, title) +
    geom_point(aes(colour = .data[[column]]), size = 0.45, alpha = 0.85) +
    scale_colour_manual(values = c(Singlet = "#bdbdbd", Doublet = "#d7301f"))
}
make_score_plot <- function(subdat, study, study_parameters) {
  ggplot(subdat, aes(x = Scrublet_score)) +
    geom_histogram(bins = 80, fill = "#9ecae1", colour = "white", linewidth = 0.08) +
    geom_vline(
      data = study_parameters[is.finite(study_parameters$recalibrated_threshold), , drop = FALSE],
      aes(xintercept = recalibrated_threshold), colour = "#c51b7d", linewidth = 0.45
    ) +
    labs(
      title = paste0(study, ": recalibrated Scrublet score distribution"),
      subtitle = "Purple lines: sample-proxy partition thresholds; pooled small samples are documented in the parameter table",
      x = "Scrublet score", y = "Cell count"
    ) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 7.5, hjust = 0.5))
}

plot_store <- list(MG0 = list(), DoubletFinder = list(), scDblFinder = list(), Scrublet_recalibrated = list())
score_store <- list()
for (study in study_order) {
  subdat <- dat[as.character(dat$Study) == study, , drop = FALSE]
  if (!nrow(subdat)) stop("No cells found for ", study)
  study_safe <- safe_name(study)
  plots <- list(
    MG0 = make_mg0_plot(subdat, paste0(study, ": historical MG0")),
    DoubletFinder = make_doublet_plot(subdat, "DoubletFinder_call", paste0(study, ": DoubletFinder")),
    scDblFinder = make_doublet_plot(subdat, "scDblFinder_call", paste0(study, ": scDblFinder")),
    Scrublet_recalibrated = make_doublet_plot(subdat, "Scrublet_recalibrated_call", paste0(study, ": recalibrated Scrublet"))
  )
  for (plot_name in names(plots)) {
    plot_store[[plot_name]][[study]] <- plots[[plot_name]]
    save_plot_both(plots[[plot_name]], file.path(figure_dir, paste0("UMAP_", study_safe, "_", plot_name)), 5, 4.5)
  }
  panel <- (plots$MG0 + plots$DoubletFinder) / (plots$scDblFinder + plots$Scrublet_recalibrated) +
    plot_annotation(
      title = paste0(study, ": historical MG0 and all three doublet callers"),
      subtitle = unique(subdat$historical_doublet_status)
    )
  save_plot_both(panel, file.path(figure_dir, paste0("Four_panel_", study_safe, "_historical_MG0_and_doublet_callers")), 10, 9)
  print(panel)
  score_plot <- make_score_plot(subdat, study, parameters[parameters$Study == study, , drop = FALSE])
  score_store[[study]] <- score_plot
  save_plot_both(score_plot, file.path(figure_dir, paste0("Scrublet_score_distribution_", study_safe)), 7, 5)
  print(score_plot)
}

composite_all_methods <- wrap_plots(lapply(study_order, function(study) {
  (plot_store$MG0[[study]] + plot_store$DoubletFinder[[study]]) /
    (plot_store$scDblFinder[[study]] + plot_store$Scrublet_recalibrated[[study]]) +
    plot_annotation(title = study)
}), ncol = 3) +
  plot_annotation(title = "Historical MG0 validation: MG0 and three recalibrated doublet callers")
save_plot_both(composite_all_methods, file.path(figure_dir, "Composite_historical_MG0_all_three_doublet_callers_six_studies"), 18, 18)
print(composite_all_methods)

long_calls <- rbind(
  data.frame(Cell = dat$Cell, Study = dat$Study, CellType1 = dat$CellType1, method = "DoubletFinder", call = dat$DoubletFinder_call),
  data.frame(Cell = dat$Cell, Study = dat$Study, CellType1 = dat$CellType1, method = "scDblFinder", call = dat$scDblFinder_call),
  data.frame(Cell = dat$Cell, Study = dat$Study, CellType1 = dat$CellType1, method = "Scrublet_recalibrated", call = dat$Scrublet_recalibrated_call)
)
long_calls$method <- factor(long_calls$method, levels = c("DoubletFinder", "scDblFinder", "Scrublet_recalibrated"))
long_calls$call <- factor(long_calls$call, levels = c("Singlet", "Doublet"))
proportion_table <- aggregate(as.integer(long_calls$call == "Doublet"), by = list(Study = long_calls$Study, CellType1 = long_calls$CellType1, method = long_calls$method), FUN = sum)
colnames(proportion_table)[4] <- "n_doublet"
denominator_table <- aggregate(Cell ~ Study + CellType1 + method, data = long_calls, FUN = length)
colnames(denominator_table)[4] <- "n_cells"
proportion_table <- merge(proportion_table, denominator_table, by = c("Study", "CellType1", "method"), sort = FALSE)
proportion_table$doublet_proportion <- proportion_table$n_doublet / proportion_table$n_cells
pooled_table <- aggregate(cbind(n_doublet, n_cells) ~ CellType1 + method, data = proportion_table, FUN = sum)
pooled_table$doublet_proportion <- pooled_table$n_doublet / pooled_table$n_cells
pooled_table$CellType1 <- factor(pooled_table$CellType1, levels = paste0("MG", 0:5))
write.csv(proportion_table, file.path(table_dir, "recalibrated_doublet_proportions_by_study_and_MG_state.csv"), row.names = FALSE)
write.csv(pooled_table, file.path(table_dir, "recalibrated_doublet_proportions_pooled_MG0_MG5.csv"), row.names = FALSE)

make_barplot <- function(subdat, title) {
  ggplot(subdat, aes(CellType1, doublet_proportion, fill = CellType1)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", 100 * doublet_proportion, n_cells)), vjust = -0.18, size = 3) +
    scale_fill_manual(values = c(MG0 = "#8dd3c7", MG1 = "#ffffb3", MG2 = "#bebada", MG3 = "#fb8072", MG4 = "#80b1d3", MG5 = "#fdb462"), drop = FALSE) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = expansion(mult = c(0, 0.16))) +
    labs(title = title, x = "Historical MG state", y = "Predicted doublet proportion") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none")
}

barplots <- list()
for (method in levels(pooled_table$method)) {
  method_title <- switch(method,
    DoubletFinder = "DoubletFinder: doublet proportion by MG state",
    scDblFinder = "scDblFinder: doublet proportion by MG state",
    Scrublet_recalibrated = "Scrublet (recalibrated): doublet proportion by MG state"
  )
  plot <- make_barplot(pooled_table[pooled_table$method == method, , drop = FALSE], method_title)
  barplots[[method]] <- plot
  save_plot_both(plot, file.path(figure_dir, paste0("Doublet_proportion_historical_MG0_MG5_", method)), 7, 5)
}
all_barplots <- wrap_plots(barplots, ncol = 3) + plot_annotation(title = "Historical MG0--MG5: three caller-specific doublet proportions")
save_plot_both(all_barplots, file.path(figure_dir, "Composite_doublet_proportion_historical_MG0_MG5_all_three_methods"), 18, 5.5)
print(all_barplots)

writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "recalibrated_Scrublet_plotting_R_sessionInfo.txt"))
writeLines(c(
  "Scope: the same QC-passed historical MG0--MG5 diagnostic cells used by the original three-method diagnostic branch.",
  "DoubletFinder and scDblFinder calls are retained verbatim from that branch.",
  "Scrublet calls come from the non-destructive recalibrated branch; they use Sample as the narrowest available technical proxy, pooling only samples with fewer than 300 cells for scoring.",
  "A recalibrated Scrublet threshold selects the expected number of high-score cells per partition. It is a diagnostic call and does not delete, relabel, or overwrite any cell.",
  "Historical source-pipeline doublet status is documented separately and must not be confused with this project's newly calculated calls."
), file.path(runtime_dir, "recalibrated_Scrublet_plotting_manifest.txt"))
message("Recalibrated historical Scrublet figures are ready: ", figure_dir)
