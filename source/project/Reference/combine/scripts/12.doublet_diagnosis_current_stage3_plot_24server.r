#!/usr/bin/env Rscript

# Stage 3, part 3: join Scrublet calls, render six-study diagnostics, and
# save all individual and composite figures.  No analytical object is changed.

rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
invisible(gc())

suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(qs) })

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
result_root <- file.path(project_root, "combine/results/full_count_reanalysis/03_current_MG0_doublet_diagnosis")
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
runtime_dir <- file.path(result_root, "runtime")
scrublet_dir <- file.path(result_root, "scrublet_output")
object_file <- file.path(result_root, "object", "current_MG0_QC_before_doublet_removal_three_method_partial.qs")
stage2_object_file <- file.path(project_root, "combine/results/full_count_reanalysis/02_microglia_integration/object/integrated_Microglia_processed.qs")
for (d in c(figure_dir, table_dir, runtime_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
         dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
         device = grDevices::pdf, bg = "white", limitsize = FALSE)
}

if (!file.exists(object_file)) stop("Stage-1 diagnostic object is missing.")
if (!dir.exists(scrublet_dir)) stop("Scrublet output directory is missing; run stage 2 first.")

diagnostic <- qs::qread(object_file)
meta <- diagnostic[[]]
meta$Cell <- rownames(meta)
required <- c("Study", "diagnostic_cluster", "MG0_status", "DoubletFinder_call", "scDblFinder_call", "DF_partition")
if (length(setdiff(required, colnames(meta)))) stop("Stage-1 metadata is incomplete.")
scrublet_files <- list.files(scrublet_dir, pattern = "^Scrublet_.*\\.csv$", full.names = TRUE)
# The run summary uses the same prefix but contains one row per partition rather
# than one row per cell, so it must not enter the per-cell concatenation.
scrublet_files <- scrublet_files[basename(scrublet_files) != "Scrublet_run_summary.csv"]
if (!length(scrublet_files)) stop("No Scrublet result files found.")
scrublet <- do.call(rbind, lapply(scrublet_files, read.csv, stringsAsFactors = FALSE, check.names = FALSE))
if (anyDuplicated(scrublet$Cell) || !setequal(scrublet$Cell, meta$Cell)) stop("Scrublet calls do not exactly match the QC-passed diagnostic cells.")
scrublet <- scrublet[match(meta$Cell, scrublet$Cell), , drop = FALSE]
meta$Scrublet_score <- scrublet$Scrublet_score
meta$Scrublet_call <- factor(scrublet$Scrublet_call, levels = c("Singlet", "Doublet"))
meta$Scrublet_threshold <- scrublet$Scrublet_threshold
meta$DoubletFinder_call <- factor(meta$DoubletFinder_call, levels = c("Singlet", "Doublet"))
meta$scDblFinder_call <- factor(meta$scDblFinder_call, levels = c("Singlet", "Doublet"))
meta$MG0_status <- factor(meta$MG0_status, levels = c("Other predicted cluster", "MG0 (manual cluster 2 transfer)"))
umap <- Seurat::Embeddings(diagnostic, "umap_diagnostic")
meta$UMAP_1 <- umap[meta$Cell, 1]
meta$UMAP_2 <- umap[meta$Cell, 2]

study_order <- c("Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
                 "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024")
study_order <- study_order[study_order %in% unique(as.character(meta$Study))]
if (length(study_order) != 6L) stop("Expected all six studies in the diagnostic object.")

umap_base <- function(dat, title) {
  ggplot(dat, aes(UMAP_1, UMAP_2)) + coord_equal() +
    labs(title = title, x = "Diagnostic UMAP 1", y = "Diagnostic UMAP 2") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.title = element_blank())
}
make_cluster_plot <- function(dat, title) umap_base(dat, title) +
  geom_point(aes(colour = factor(diagnostic_cluster)), size = 0.45, alpha = 0.82) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
make_mg0_plot <- function(dat, title) {
  dat <- dat[order(dat$MG0_status), , drop = FALSE]
  umap_base(dat, title) + geom_point(aes(colour = MG0_status), size = 0.45, alpha = 0.86) +
    scale_colour_manual(values = c("Other predicted cluster" = "#d9d9d9", "MG0 (manual cluster 2 transfer)" = "#c51b7d"))
}
make_doublet_plot <- function(dat, column, title) {
  dat <- dat[order(dat[[column]]), , drop = FALSE]
  umap_base(dat, title) + geom_point(aes(colour = .data[[column]]), size = 0.45, alpha = 0.86) +
    scale_colour_manual(values = c(Singlet = "#bdbdbd", Doublet = "#d7301f"))
}

plot_types <- c("diagnostic_cluster", "MG0", "DoubletFinder", "scDblFinder", "Scrublet")
plot_store <- setNames(vector("list", length(plot_types)), plot_types)
for (plot_type in plot_types) plot_store[[plot_type]] <- list()
for (study in study_order) {
  dat <- meta[meta$Study == study, , drop = FALSE]
  plots <- list(
    diagnostic_cluster = make_cluster_plot(dat, paste0(study, ": diagnostic clusters")),
    MG0 = make_mg0_plot(dat, paste0(study, ": MG0 transferred from manual cluster 2")),
    DoubletFinder = make_doublet_plot(dat, "DoubletFinder_call", paste0(study, ": DoubletFinder")),
    scDblFinder = make_doublet_plot(dat, "scDblFinder_call", paste0(study, ": scDblFinder")),
    Scrublet = make_doublet_plot(dat, "Scrublet_call", paste0(study, ": Scrublet"))
  )
  for (plot_type in names(plots)) {
    # Keep the individual figures fully titled, but use compact study-only
    # titles inside six-study composites so that three-column panels remain
    # legible rather than clipping long method descriptions.
    plot_store[[plot_type]][[study]] <- plots[[plot_type]] +
      labs(title = study) +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
    save_plot_both(plots[[plot_type]], file.path(figure_dir, paste0("UMAP_", safe_name(study), "_", plot_type)), 5, 4.5)
  }
}
for (plot_type in plot_types) {
  composite <- wrap_plots(plot_store[[plot_type]][study_order], ncol = 3) +
    plot_annotation(title = paste0("Six-study QC-passed diagnostic UMAPs: ", plot_type))
  save_plot_both(composite, file.path(figure_dir, paste0("Composite_UMAP_", plot_type, "_six_studies")), 15, 9)
  print(composite)
}

long_calls <- rbind(
  data.frame(Cell = meta$Cell, Study = meta$Study, MG0_status = meta$MG0_status, method = "DoubletFinder", call = meta$DoubletFinder_call),
  data.frame(Cell = meta$Cell, Study = meta$Study, MG0_status = meta$MG0_status, method = "scDblFinder", call = meta$scDblFinder_call),
  data.frame(Cell = meta$Cell, Study = meta$Study, MG0_status = meta$MG0_status, method = "Scrublet", call = meta$Scrublet_call)
)
long_calls$method <- factor(long_calls$method, levels = c("DoubletFinder", "scDblFinder", "Scrublet"))
long_calls$call <- factor(long_calls$call, levels = c("Singlet", "Doublet"))
proportion_table <- aggregate(as.integer(long_calls$call == "Doublet"), by = list(Study = long_calls$Study, MG0_status = long_calls$MG0_status, method = long_calls$method), FUN = sum)
names(proportion_table)[4] <- "n_doublet"
denominator_table <- aggregate(Cell ~ Study + MG0_status + method, data = long_calls, FUN = length)
names(denominator_table)[4] <- "n_cells"
proportion_table <- merge(proportion_table, denominator_table, by = c("Study", "MG0_status", "method"), sort = FALSE)
proportion_table$doublet_proportion <- proportion_table$n_doublet / proportion_table$n_cells
pooled_table <- aggregate(cbind(n_doublet, n_cells) ~ MG0_status + method, data = proportion_table, FUN = sum)
pooled_table$doublet_proportion <- pooled_table$n_doublet / pooled_table$n_cells

make_barplot <- function(dat) {
  ggplot(dat, aes(MG0_status, doublet_proportion, fill = MG0_status)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.1f%%\\n(n=%d)", 100 * doublet_proportion, n_cells)), vjust = -0.18, size = 3) +
    facet_wrap(~method, nrow = 1) +
    scale_fill_manual(values = c("Other predicted cluster" = "#bdbdbd", "MG0 (manual cluster 2 transfer)" = "#c51b7d")) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = expansion(mult = c(0, 0.17))) +
    labs(title = "Predicted doublet proportion: current MG0 versus other reanalysis clusters",
         subtitle = "MG0 is the user-selected final cluster 2, transferred to QC-passed pre-removal cells.",
         x = NULL, y = "Predicted doublet proportion") +
    theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none")
}
p_current_bar <- make_barplot(pooled_table)
save_plot_both(p_current_bar, file.path(figure_dir, "Doublet_proportion_current_MG0_vs_other_three_methods"), 15, 5)
print(p_current_bar)

# The reanalysis has a manual definition only for MG0.  It has no current
# MG1--MG5 assignment, so an MG0--MG5 plot would be falsely labelled.  Write a
# transparent selection template for the later user-defined state assignments.
reference <- qs::qread(stage2_object_file)
state_template <- as.data.frame(table(integrated_cluster = as.character(reference$integrated_cluster)), stringsAsFactors = FALSE)
names(state_template)[2] <- "n_cells"
state_template$assign_MG_state <- ifelse(state_template$integrated_cluster == "2", "MG0", NA_character_)
state_template$user_note <- ifelse(state_template$integrated_cluster == "2", "User-defined MG0 (2026-09-07)", NA_character_)
write.csv(state_template, file.path(table_dir, "current_integrated_MG0_MG5_state_selection_TEMPLATE.csv"), row.names = FALSE)

write.csv(meta, file.path(table_dir, "six_study_cell_metadata_with_all_three_doublet_calls.csv"), row.names = FALSE)
write.csv(proportion_table, file.path(table_dir, "doublet_proportions_by_study_current_MG0_vs_other.csv"), row.names = FALSE)
write.csv(pooled_table, file.path(table_dir, "doublet_proportions_pooled_current_MG0_vs_other.csv"), row.names = FALSE)
write.csv(long_calls, file.path(table_dir, "per_cell_long_three_method_calls.csv"), row.names = FALSE)
agreement <- meta[, c("Cell", "Study", "MG0_status", "DoubletFinder_call", "scDblFinder_call", "Scrublet_call"), drop = FALSE]
agreement$n_methods_calling_doublet <- rowSums(cbind(agreement$DoubletFinder_call == "Doublet", agreement$scDblFinder_call == "Doublet", agreement$Scrublet_call == "Doublet"), na.rm = TRUE)
write.csv(agreement, file.path(table_dir, "per_cell_three_method_agreement.csv"), row.names = FALSE)
writeLines(c(
  "Input: Stage-2 selected-Microglia object after repeated basic QC and before DoubletFinder removal (16,931 cells).",
  "Diagnostic UMAP: existing uncorrected PCA-neighbour UMAP built after basic QC from shared cross-study features; no batch correction was applied to this diagnostic embedding.",
  "DoubletFinder: recorded Stage-2 calls on the exact input/technical partitions. scDblFinder and Scrublet: newly computed independently on the same raw-count technical partitions.",
  "MG0 overlay: manual final integrated cluster 2 transferred to pre-removal diagnostic cells; scores are retained. It is not a source-object CellType1 label.",
  "The requested current MG0--MG5 bar plot is intentionally withheld: MG1--MG5 have not yet been manually defined in the new reanalysis. The current MG0-vs-other bar plot is not a substitute for it.",
  "No doublet call deletes or relabels a cell in this diagnostic analysis."
), file.path(runtime_dir, "stage3_method_manifest.txt"))
message("Stage 3 plotting complete. MG0-vs-other diagnostics are ready; current MG1--MG5 labels are still required for a correctly labelled MG0--MG5 bar plot.")
