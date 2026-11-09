#!/usr/bin/env Rscript

# Artificial-doublet visualisation, revision 3.
#
# This plotting-only revision splits the prior combined overlay into TWO full
# figures per study and per caller, with all coordinates reused unchanged:
#   A. artificial profiles made by that caller over grey real-cell context;
#   B. that caller's labelled real doublets over grey real-cell context.
# Artificial profiles and caller-labelled real doublets are deliberately never
# coloured together in the same figure. No counts, doublet scores/calls, or
# dimensional reductions are recalculated.

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
combine_root <- file.path(project_root, "combine")
result_root <- file.path(combine_root, "results/artificial_doublet_visualization")
object_dir <- file.path(result_root, "object")
figure_dir <- file.path(result_root, "figure")
table_dir <- file.path(result_root, "table")
runtime_dir <- file.path(result_root, "runtime")
for (directory in c(figure_dir, table_dir, runtime_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

study_order <- c(
  "Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
  "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"
)
callers <- data.frame(
  method = c("DoubletFinder", "scDblFinder", "Scrublet"),
  artificial_origin = c("DoubletFinder artificial", "scDblFinder artificial", "Scrublet artificial"),
  call_column = c("DoubletFinder_call", "scDblFinder_call", "Scrublet_call"),
  artificial_colour = c("#E64B35", "#4DBBD5", "#00A087"),
  stringsAsFactors = FALSE
)

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
save_plot_both <- function(plot, stem, width = 8, height = 6) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
    dpi = 320, bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
    device = grDevices::pdf, bg = "white")
}

add_embedding <- function(object, study, call_column) {
  needed <- c("profile_origin", call_column)
  if (!all(needed %in% colnames(object[[]]))) {
    stop(study, ": required metadata absent: ", paste(setdiff(needed, colnames(object[[]])), collapse = ", "))
  }
  embedding <- Seurat::Embeddings(object, "umap_artificial_doublet")
  dat <- object[[]]
  dat$Cell <- rownames(dat)
  dat$Artificial_UMAP_1 <- embedding[dat$Cell, 1]
  dat$Artificial_UMAP_2 <- embedding[dat$Cell, 2]
  dat
}

base_theme <- function() {
  theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 10)
    )
}

make_artificial_only_plot <- function(dat, study, method, artificial_origin, artificial_colour) {
  display <- dat[as.character(dat$profile_origin) %in% c("Real cell", artificial_origin), , drop = FALSE]
  display$display_class <- ifelse(
    as.character(display$profile_origin) == artificial_origin,
    paste0(method, " artificial doublet"),
    "Real cell"
  )
  levels <- c("Real cell", paste0(method, " artificial doublet"))
  display$display_class <- factor(display$display_class, levels = levels)
  display <- display[order(display$display_class), , drop = FALSE]

  ggplot(display, aes(Artificial_UMAP_1, Artificial_UMAP_2)) +
    geom_point(aes(colour = display_class), size = 0.48, alpha = 0.88) +
    scale_colour_manual(values = setNames(c("#BDBDBD", artificial_colour), levels)) +
    coord_equal() +
    labs(
      title = paste0(study, ": ", method),
      subtitle = "Caller-specific artificial doublets",
      x = "Artificial-doublet UMAP 1", y = "Artificial-doublet UMAP 2", colour = NULL
    ) +
    base_theme() +
    guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1)))
}

make_called_real_doublet_only_plot <- function(dat, study, method, call_column) {
  display <- dat[as.character(dat$profile_origin) == "Real cell", , drop = FALSE]
  display$display_class <- ifelse(
    as.character(display[[call_column]]) == "Doublet",
    paste0(method, "-called real doublet"),
    "Other real cell"
  )
  levels <- c("Other real cell", paste0(method, "-called real doublet"))
  display$display_class <- factor(display$display_class, levels = levels)
  display <- display[order(display$display_class), , drop = FALSE]

  ggplot(display, aes(Artificial_UMAP_1, Artificial_UMAP_2)) +
    geom_point(aes(colour = display_class), size = 0.58, alpha = 0.88) +
    scale_colour_manual(values = setNames(c("#BDBDBD", "#7B3294"), levels)) +
    coord_equal() +
    labs(
      title = paste0(study, ": ", method),
      subtitle = "Caller-labelled real doublets only",
      x = "Artificial-doublet UMAP 1", y = "Artificial-doublet UMAP 2", colour = NULL
    ) +
    base_theme() +
    guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1)))
}

pair_rows <- list()
row_id <- 0L
for (study in study_order) {
  object_file <- file.path(object_dir, paste0(safe_name(study), "_real_and_artificial_doublet_UMAP.qs"))
  if (!file.exists(object_file)) stop("Missing study-specific artificial-doublet object: ", object_file)
  object <- qs::qread(object_file)
  for (caller_id in seq_len(nrow(callers))) {
    caller <- callers[caller_id, , drop = FALSE]
    dat <- add_embedding(object, study, caller$call_column)

    artificial_plot <- make_artificial_only_plot(
      dat, study, caller$method, caller$artificial_origin, caller$artificial_colour
    )
    artificial_stem <- file.path(figure_dir, paste0(
      "UMAP_", safe_name(study), "_", caller$method, "_artificial_doublets_only"
    ))
    save_plot_both(artificial_plot, artificial_stem)
    print(artificial_plot)

    real_plot <- make_called_real_doublet_only_plot(dat, study, caller$method, caller$call_column)
    real_stem <- file.path(figure_dir, paste0(
      "UMAP_", safe_name(study), "_", caller$method, "_called_real_doublets_only"
    ))
    save_plot_both(real_plot, real_stem)
    print(real_plot)

    for (item in list(
      list(kind = "artificial_doublets_only", stem = artificial_stem),
      list(kind = "called_real_doublets_only", stem = real_stem)
    )) {
      row_id <- row_id + 1L
      pair_rows[[row_id]] <- data.frame(
        Study = study,
        method = caller$method,
        figure_type = item$kind,
        figure_stem = item$stem,
        PNG_exists = file.exists(paste0(item$stem, ".png")),
        PDF_exists = file.exists(paste0(item$stem, ".pdf")),
        stringsAsFactors = FALSE
      )
    }
  }
  rm(object)
  invisible(gc())
}

pair_check <- do.call(rbind, pair_rows)
write.csv(pair_check,
  file.path(table_dir, "per_method_artificial_and_called_real_doublet_separate_PNG_PDF_pair_check.csv"),
  row.names = FALSE
)
if (any(!pair_check$PNG_exists | !pair_check$PDF_exists)) {
  stop("At least one separate artificial/real-doublet figure lacks a PNG/PDF pair.")
}
writeLines(capture.output(sessionInfo()),
  file.path(runtime_dir, "per_method_artificial_and_called_real_doublet_separate_R_sessionInfo.txt")
)
message("Separate artificial-doublet and caller-labelled real-doublet plots are complete: ", figure_dir)
