#!/usr/bin/env Rscript

# Corrected full-cell diagnosis, stage 3.
# Calls are joined to each complete QC-passed study object, then historical
# MG0--MG5 is used only as a post hoc label for the requested proportions.

rm(list = ls(envir = .GlobalEnv, all.names = TRUE), envir = .GlobalEnv)
invisible(gc())
suppressPackageStartupMessages({ library(Seurat); library(qs); library(ggplot2) })

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
result_root <- file.path(project_root, "combine/results/fullcell_historical_MG0_MG5_doublet_diagnosis_v2")
object_dir <- file.path(result_root, "object")
input_dir <- file.path(result_root, "scrublet_input")
scrublet_dir <- file.path(result_root, "scrublet_output")
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
runtime_dir <- file.path(result_root, "runtime")
for (directory in c(figure_dir, table_dir, runtime_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
study_order <- c("Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
                 "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024")
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width = 7.5, height = 5) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 320, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = grDevices::pdf, bg = "white", limitsize = FALSE)
}
plot_base <- function(dat, title) {
  ggplot(dat, aes(UMAP_1, UMAP_2)) + coord_equal() +
    labs(title = title, subtitle = "All QC-passed real cells in this study; historical MG labels only overlay after calling", x = "Full-cell diagnostic UMAP 1", y = "Full-cell diagnostic UMAP 2") +
    theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = .5), plot.subtitle = element_text(hjust = .5, size = 8.5), legend.title = element_blank(), legend.position = "bottom")
}
make_cluster <- function(dat, title) plot_base(dat, title) + geom_point(aes(colour = factor(diagnostic_cluster)), size = .40, alpha = .80) + guides(colour = guide_legend(override.aes = list(size = 2)))
make_label <- function(dat, field, positive, colour, title) {
  dat$overlay <- ifelse(as.logical(dat[[field]]), positive, "Other QC-passed real cells")
  dat$overlay <- factor(dat$overlay, levels = c("Other QC-passed real cells", positive))
  dat <- dat[order(dat$overlay), , drop = FALSE]
  plot_base(dat, title) + geom_point(aes(colour = overlay), size = .42, alpha = .85) +
    scale_colour_manual(values = c("Other QC-passed real cells" = "#D9D9D9", setNames(colour, positive))) + guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
}
make_call <- function(dat, field, title) {
  dat$overlay <- factor(dat[[field]], levels = c("Singlet", "Doublet")); dat <- dat[order(dat$overlay), , drop = FALSE]
  plot_base(dat, title) + geom_point(aes(colour = overlay), size = .42, alpha = .85) +
    scale_colour_manual(values = c(Singlet = "#BDBDBD", Doublet = "#D7301F")) + guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))
}

manifest <- read.csv(file.path(input_dir, "scrublet_manifest.csv"), stringsAsFactors = FALSE, check.names = FALSE)
metadata_all <- list(); output_stems <- character()
for (study in study_order) {
  safe <- safe_name(study)
  partial_file <- file.path(object_dir, paste0(safe, "_all_QC_passed_fullcell_partial_calls.qs"))
  if (!file.exists(partial_file)) stop("Missing corrected full-cell partial object for ", study)
  object <- qs::qread(partial_file)
  meta <- object[[]]; meta$Cell <- rownames(meta)
  files <- file.path(scrublet_dir, paste0("Scrublet_", manifest$safe_partition[manifest$Study == study], ".csv"))
  if (!length(files) || any(!file.exists(files))) stop("Missing corrected full-cell Scrublet calls for ", study)
  scrublet <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE, check.names = FALSE))
  if (anyDuplicated(scrublet$Cell) || !setequal(scrublet$Cell, meta$Cell)) stop(study, ": Scrublet calls do not exactly match all QC-passed cells.")
  scrublet <- scrublet[match(meta$Cell, scrublet$Cell), , drop = FALSE]
  object$Scrublet_score <- scrublet$Scrublet_score
  object$Scrublet_call <- factor(scrublet$Scrublet_call, levels = c("Singlet", "Doublet"))
  object$Scrublet_threshold <- scrublet$Scrublet_threshold
  object$DoubletFinder_call <- factor(object$DoubletFinder_call, levels = c("Singlet", "Doublet"))
  object$scDblFinder_call <- factor(object$scDblFinder_call, levels = c("Singlet", "Doublet"))
  umap <- Seurat::Embeddings(object, "umap_fullcell_diagnostic")
  meta <- object[[]]; meta$Cell <- rownames(meta); meta$UMAP_1 <- umap[meta$Cell, 1]; meta$UMAP_2 <- umap[meta$Cell, 2]
  meta$historical_MG0 <- as.logical(meta$historical_MG0); meta$MG00 <- as.logical(meta$MG00)
  plots <- list(
    diagnostic_cluster = make_cluster(meta, paste0(study, ": full-cell diagnostic clusters")),
    historical_MG0 = make_label(meta, "historical_MG0", "MG0", "#2C7FB8", paste0(study, ": MG0 (original MG0--MG5 label)")),
    MG00 = make_label(meta, "MG00", "MG00", "#C51B7D", paste0(study, ": MG00")),
    DoubletFinder = make_call(meta, "DoubletFinder_call", paste0(study, ": DoubletFinder")),
    scDblFinder = make_call(meta, "scDblFinder_call", paste0(study, ": scDblFinder")),
    Scrublet = make_call(meta, "Scrublet_call", paste0(study, ": Scrublet"))
  )
  for (name in names(plots)) {
    stem <- file.path(figure_dir, paste0("UMAP_", safe, "_", name)); save_plot_both(plots[[name]], stem); print(plots[[name]]); output_stems <- c(output_stems, stem)
  }
  qs::qsave(object, file.path(object_dir, paste0(safe, "_all_QC_passed_fullcell_all_three_calls.qs")), preset = "high")
  write.csv(meta, file.path(table_dir, paste0(safe, "_all_QC_passed_cells_with_three_calls_and_label_mapping.csv")), row.names = FALSE)
  # Source objects carry study-specific annotation columns.  Only retain the
  # common audit fields required for the cross-study calls and barcode mapping;
  # this prevents a column mismatch from changing any per-study output.
  common_audit_fields <- c("Cell", "Study", "historical_MG_state", "historical_MG0", "MG00",
                           "DoubletFinder_call", "scDblFinder_call", "Scrublet_call", "Scrublet_score",
                           "Scrublet_threshold", "DF_partition", "diagnostic_cluster")
  metadata_all[[study]] <- meta[, common_audit_fields, drop = FALSE]
  rm(object, meta, scrublet, umap); invisible(gc())
}

all_meta <- do.call(rbind, metadata_all)
long_calls <- rbind(
  data.frame(Cell = all_meta$Cell, Study = all_meta$Study, historical_MG_state = all_meta$historical_MG_state, method = "DoubletFinder", call = all_meta$DoubletFinder_call),
  data.frame(Cell = all_meta$Cell, Study = all_meta$Study, historical_MG_state = all_meta$historical_MG_state, method = "scDblFinder", call = all_meta$scDblFinder_call),
  data.frame(Cell = all_meta$Cell, Study = all_meta$Study, historical_MG_state = all_meta$historical_MG_state, method = "Scrublet", call = all_meta$Scrublet_call)
)
long_calls <- long_calls[!is.na(long_calls$historical_MG_state), , drop = FALSE]
long_calls$historical_MG_state <- factor(as.character(long_calls$historical_MG_state), levels = paste0("MG", 0:5))
long_calls$method <- factor(long_calls$method, levels = c("DoubletFinder", "scDblFinder", "Scrublet"))
long_calls$call <- factor(long_calls$call, levels = c("Singlet", "Doublet"))
summary_by_study <- aggregate(as.integer(long_calls$call == "Doublet"), by = list(Study = long_calls$Study, historical_MG_state = long_calls$historical_MG_state, method = long_calls$method), sum)
names(summary_by_study)[4] <- "n_doublet"
denom_by_study <- aggregate(Cell ~ Study + historical_MG_state + method, data = long_calls, length); names(denom_by_study)[4] <- "n_cells"
summary_by_study <- merge(summary_by_study, denom_by_study, by = c("Study", "historical_MG_state", "method"), sort = FALSE)
summary_by_study$doublet_proportion <- summary_by_study$n_doublet / summary_by_study$n_cells
summary_pooled <- aggregate(cbind(n_doublet, n_cells) ~ historical_MG_state + method, data = summary_by_study, sum)
summary_pooled$doublet_proportion <- summary_pooled$n_doublet / summary_pooled$n_cells
for (method in levels(summary_pooled$method)) {
  dat <- summary_pooled[summary_pooled$method == method, , drop = FALSE]
  plot <- ggplot(dat, aes(historical_MG_state, doublet_proportion, fill = historical_MG_state)) +
    geom_col(width = .72) + geom_text(aes(label = sprintf("%.1f%%\\n(n=%d)", 100 * doublet_proportion, n_cells)), vjust = -.18, size = 3.4) +
    scale_fill_manual(values = c(MG0 = "#C51B7D", MG1 = "#FEE08B", MG2 = "#91BFDB", MG3 = "#FDAE61", MG4 = "#ABDDA4", MG5 = "#BDBDBD")) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), expand = expansion(mult = c(0, .17))) +
    labs(title = paste0("MG0--MG5 doublet proportion: ", method), subtitle = "Calls computed on every QC-passed real cell first, then restricted to cells mapped to the original MG0--MG5 labels for this summary", x = NULL, y = "Predicted doublet proportion") +
    theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = .5), plot.subtitle = element_text(hjust = .5, size = 8.5), legend.position = "none")
  stem <- file.path(figure_dir, paste0("Doublet_proportion_historical_MG0_MG5_", method, "_pooled_fullcell_calls")); save_plot_both(plot, stem, 8, 5); print(plot); output_stems <- c(output_stems, stem)
}
mapping_coverage <- aggregate(Cell ~ Study, data = all_meta, length); names(mapping_coverage)[2] <- "n_all_QC_passed_real_cells"
mapped <- aggregate(!is.na(all_meta$historical_MG_state), by = list(Study = all_meta$Study), sum); names(mapped)[2] <- "n_mapped_historical_MG0_MG5"
mapping_coverage <- merge(mapping_coverage, mapped, by = "Study"); mapping_coverage$mapping_fraction <- mapping_coverage$n_mapped_historical_MG0_MG5 / mapping_coverage$n_all_QC_passed_real_cells
mapping_coverage$Study <- factor(mapping_coverage$Study, levels = study_order)
coverage_plot <- ggplot(mapping_coverage, aes(Study, mapping_fraction)) +
  geom_col(fill = "#2C7FB8", width = .70) +
  geom_text(aes(label = sprintf("%.1f%%\\n(%d/%d)", 100 * mapping_fraction, n_mapped_historical_MG0_MG5, n_all_QC_passed_real_cells)), vjust = -.18, size = 3.1) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.05), expand = expansion(mult = c(0, .02))) +
  labs(title = "Original MG0--MG5 barcode mapping coverage", subtitle = "Of all QC-passed real cells in each study; this mapping is applied only after full-cell calls", x = NULL, y = "Mapped fraction") +
  theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = .5), plot.subtitle = element_text(hjust = .5, size = 8.5), axis.text.x = element_text(angle = 25, hjust = 1))
coverage_stem <- file.path(figure_dir, "Original_MG0_MG5_barcode_mapping_coverage_among_all_QC_passed_cells")
save_plot_both(coverage_plot, coverage_stem, 8.5, 5); print(coverage_plot); output_stems <- c(output_stems, coverage_stem)
pair_check <- data.frame(figure_stem = output_stems, PNG_exists = file.exists(paste0(output_stems, ".png")), PDF_exists = file.exists(paste0(output_stems, ".pdf")))
if (any(!pair_check$PNG_exists | !pair_check$PDF_exists)) stop("A corrected full-cell plot is missing PNG or PDF.")
write.csv(all_meta, file.path(table_dir, "six_study_all_QC_passed_cells_with_all_three_calls_and_label_mapping.csv"), row.names = FALSE)
write.csv(long_calls, file.path(table_dir, "historical_MG0_MG5_cells_three_method_long_calls.csv"), row.names = FALSE)
write.csv(summary_by_study, file.path(table_dir, "doublet_proportions_by_study_historical_MG0_MG5_fullcell_calls.csv"), row.names = FALSE)
write.csv(summary_pooled, file.path(table_dir, "doublet_proportions_pooled_historical_MG0_MG5_fullcell_calls.csv"), row.names = FALSE)
write.csv(mapping_coverage, file.path(table_dir, "historical_MG0_MG5_mapping_coverage_among_all_QC_passed_cells.csv"), row.names = FALSE)
write.csv(pair_check, file.path(table_dir, "PNG_PDF_pair_check.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "stage3_R_sessionInfo.txt"))
writeLines(c(
  "Corrected scope: all real-cell doublet calls and diagnostic UMAPs use all QC-passed real cells, separately in each study.",
  "The MG0--MG5 bars use those full-cell calls only after a Study + original_barcode mapping to the original MG0--MG5 labels.",
  "The denominator n shown over each bar is the mapped historical MG-state cell count, not the full-study cell count.",
  "All individual UMAPs and bar figures were rendered in R and saved as matched PNG/PDF files."
), file.path(runtime_dir, "stage3_scope_manifest.txt"))
message("Corrected full-cell stage 3 complete: ", figure_dir)
