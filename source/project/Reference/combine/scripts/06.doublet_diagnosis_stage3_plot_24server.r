#!/usr/bin/env Rscript

# MG0--MG5 doublet diagnosis stage 3: join Scrublet and create figures/tables.
# Every individual figure is saved as PNG/PDF; the six-study composite views are
# also printed so they remain inspectable in RStudio's Plots history.

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
result_root <- file.path(project_root, "combine/results/MG0_MG5_doublet_diagnosis")
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
scrublet_dir <- file.path(result_root, "scrublet_output")
runtime_dir <- file.path(result_root, "runtime")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
study_order <- c("Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
                 "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024")
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = grDevices::pdf, bg = "white", limitsize = FALSE)
}

stage1 <- read.csv(file.path(table_dir, "stage1_cell_metadata_before_scrublet.csv"), check.names = FALSE, stringsAsFactors = FALSE)
scrublet_files <- file.path(scrublet_dir, paste0("Scrublet_", safe_name(study_order), ".csv"))
if (any(!file.exists(scrublet_files))) stop("Missing Scrublet result files; finish stage 2 first.")
scrublet <- do.call(rbind, lapply(scrublet_files, read.csv, stringsAsFactors = FALSE))
if (anyDuplicated(scrublet$Cell)) stop("Duplicate cell IDs in Scrublet results.")
idx <- match(stage1$Cell, scrublet$Cell)
if (anyNA(idx)) stop("Some stage-1 cells do not have a Scrublet result.")
if (any(stage1$Study != scrublet$Study[idx])) stop("Study labels disagree between stages.")
stage1$Scrublet_score <- scrublet$Scrublet_score[idx]
stage1$Scrublet_call <- scrublet$Scrublet_call[idx]
stage1$Scrublet_threshold <- scrublet$Scrublet_threshold[idx]
stage1$CellType1 <- factor(stage1$CellType1, levels = paste0("MG", 0:5))
stage1$MG0_status <- factor(ifelse(stage1$CellType1 == "MG0", "MG0", "Other MG state"), levels = c("Other MG state", "MG0"))
stage1$DoubletFinder_call <- factor(stage1$DoubletFinder_call, levels = c("Singlet", "Doublet"))
stage1$scDblFinder_call <- factor(stage1$scDblFinder_call, levels = c("Singlet", "Doublet"))
stage1$Scrublet_call <- factor(stage1$Scrublet_call, levels = c("Singlet", "Doublet"))
stage1$diagnostic_cluster <- factor(stage1$diagnostic_cluster)

umap_base <- function(dat, title) ggplot(dat, aes(UMAP_1, UMAP_2)) + coord_equal() +
  labs(title = title, x = "UMAP 1", y = "UMAP 2") + theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.title = element_blank())
make_cluster_plot <- function(dat, title) umap_base(dat, title) + geom_point(aes(colour = diagnostic_cluster), size = 0.45, alpha = 0.8) + guides(colour = guide_legend(override.aes = list(size = 2)))
make_mg0_plot <- function(dat, title) {
  dat <- dat[order(dat$MG0_status), , drop = FALSE]
  umap_base(dat, title) + geom_point(aes(colour = MG0_status), size = 0.45, alpha = 0.85) + scale_colour_manual(values = c("Other MG state" = "#d9d9d9", MG0 = "#c51b7d"))
}
make_doublet_plot <- function(dat, column, title) {
  dat <- dat[order(dat[[column]]), , drop = FALSE]
  umap_base(dat, title) + geom_point(aes(colour = .data[[column]]), size = 0.45, alpha = 0.85) + scale_colour_manual(values = c(Singlet = "#bdbdbd", Doublet = "#d7301f"))
}

plot_types <- c("diagnostic_cluster", "MG0", "DoubletFinder", "scDblFinder", "Scrublet")
plot_store <- setNames(vector("list", length(plot_types)), plot_types)
for (pt in plot_types) plot_store[[pt]] <- list()
for (study in study_order) {
  dat <- stage1[stage1$Study == study, , drop = FALSE]
  plots <- list(diagnostic_cluster = make_cluster_plot(dat, paste0(study, ": diagnostic clusters")),
    MG0 = make_mg0_plot(dat, paste0(study, ": MG0")),
    DoubletFinder = make_doublet_plot(dat, "DoubletFinder_call", paste0(study, ": DoubletFinder")),
    scDblFinder = make_doublet_plot(dat, "scDblFinder_call", paste0(study, ": scDblFinder")),
    Scrublet = make_doublet_plot(dat, "Scrublet_call", paste0(study, ": Scrublet")))
  for (pt in names(plots)) {
    plot_store[[pt]][[study]] <- plots[[pt]]
    save_plot_both(plots[[pt]], file.path(figure_dir, paste0("UMAP_", safe_name(study), "_", pt)), 5, 4.5)
  }
}
composite_umaps <- setNames(vector("list", length(plot_types)), plot_types)
for (pt in plot_types) {
  composite <- wrap_plots(plot_store[[pt]][study_order], ncol = 3) + plot_annotation(title = paste0("Six-study diagnostic UMAPs: ", pt))
  composite_umaps[[pt]] <- composite
  save_plot_both(composite, file.path(figure_dir, paste0("Composite_UMAP_", pt, "_six_studies")), 15, 9)
  print(composite)
}

long_calls <- rbind(data.frame(Cell = stage1$Cell, Study = stage1$Study, CellType1 = stage1$CellType1, method = "DoubletFinder", call = stage1$DoubletFinder_call),
  data.frame(Cell = stage1$Cell, Study = stage1$Study, CellType1 = stage1$CellType1, method = "scDblFinder", call = stage1$scDblFinder_call),
  data.frame(Cell = stage1$Cell, Study = stage1$Study, CellType1 = stage1$CellType1, method = "Scrublet", call = stage1$Scrublet_call))
long_calls$method <- factor(long_calls$method, levels = c("DoubletFinder", "scDblFinder", "Scrublet"))
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
utils::write.csv(stage1, file.path(table_dir, "six_study_cell_metadata_with_all_doublet_calls.csv"), row.names = FALSE)
utils::write.csv(proportion_table, file.path(table_dir, "doublet_proportions_by_study_and_MG_state.csv"), row.names = FALSE)
utils::write.csv(pooled_table, file.path(table_dir, "doublet_proportions_pooled_MG0_MG5.csv"), row.names = FALSE)

make_barplot <- function(dat, title, facet = FALSE) {
  p <- ggplot(dat, aes(CellType1, doublet_proportion, fill = CellType1)) + geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", 100 * doublet_proportion, n_cells)), vjust = -0.18, size = 3) +
    scale_fill_manual(values = c(MG0 = "#8dd3c7", MG1 = "#ffffb3", MG2 = "#bebada", MG3 = "#fb8072", MG4 = "#80b1d3", MG5 = "#fdb462"), drop = FALSE) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = expansion(mult = c(0, 0.16))) +
    labs(title = title, x = "MG state", y = "Predicted doublet proportion") + theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none")
  if (facet) p <- p + facet_wrap(~method, nrow = 1)
  p
}
barplots_by_method <- setNames(vector("list", length(levels(pooled_table$method))), levels(pooled_table$method))
for (method in levels(pooled_table$method)) {
  p <- make_barplot(pooled_table[pooled_table$method == method, , drop = FALSE], paste0(method, ": doublet proportion across MG0--MG5"))
  barplots_by_method[[method]] <- p
  save_plot_both(p, file.path(figure_dir, paste0("Doublet_proportion_MG0_MG5_", method)), 7, 5)
  print(p)
}
p_bar_composite <- make_barplot(pooled_table, "Predicted doublet proportion across MG0--MG5", facet = TRUE)
save_plot_both(p_bar_composite,
  file.path(figure_dir, "Composite_Doublet_proportion_MG0_MG5_three_methods"), 15, 5)
print(p_bar_composite)

agreement <- stage1[, c("Cell", "Study", "Sample", "CellType1", "DoubletFinder_call", "scDblFinder_call", "Scrublet_call")]
agreement$n_methods_calling_doublet <- rowSums(cbind(agreement$DoubletFinder_call == "Doublet", agreement$scDblFinder_call == "Doublet", agreement$Scrublet_call == "Doublet"), na.rm = TRUE)
utils::write.csv(agreement, file.path(table_dir, "per_cell_doublet_method_agreement.csv"), row.names = FALSE)
writeLines(c(
  "Stage 3 UMAPs originate from the stage-1 study-specific count-based diagnostic objects.",
  "Each study has individual PNG/PDF UMAPs for diagnostic clusters, MG0, DoubletFinder, scDblFinder and Scrublet. Six-study composites are supplementary summaries.",
  "Doublet proportions are pooled across studies per MG0--MG5 state; study-stratified proportions are also saved.",
  "Calls are diagnostic only: no cell is deleted or relabelled by this analysis.",
  "Limitation: after initial microglia selection, heterotypic MG--non-MG doublets may be underestimated."
), file.path(runtime_dir, "stage3_method_manifest.txt"))
message("Stage 3 complete: ", figure_dir)
