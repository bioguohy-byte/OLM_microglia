#!/usr/bin/env Rscript

# Artificial-doublet visualisation, revision 2.
#
# This plotting-only script reads the previously saved, study-specific Seurat
# objects that contain real cells plus all three classes of artificial profiles.
# For each study and caller, it makes ONE UMAP that displays only:
#   1) other real cells;
#   2) artificial profiles constructed by that caller; and
#   3) real cells labelled Doublet by that same caller.
# The common UMAP coordinates are reused unchanged.  No count matrix, doublet
# score, call, or dimensional reduction is recalculated in this revision.

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
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(runtime_dir, recursive = TRUE, showWarnings = FALSE)

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

make_method_overlay <- function(object, study, method, artificial_origin,
                                call_column, artificial_colour) {
  if (!all(c("profile_origin", call_column) %in% colnames(object[[]]))) {
    stop(study, ": required per-cell metadata is absent for ", method)
  }
  embedding <- Seurat::Embeddings(object, "umap_artificial_doublet")
  dat <- object[[]]
  dat$Cell <- rownames(dat)
  dat$Artificial_UMAP_1 <- embedding[dat$Cell, 1]
  dat$Artificial_UMAP_2 <- embedding[dat$Cell, 2]
  dat <- dat[as.character(dat$profile_origin) %in% c("Real cell", artificial_origin), , drop = FALSE]
  real_doublet <- as.character(dat[[call_column]]) == "Doublet"
  dat$display_class <- ifelse(
    as.character(dat$profile_origin) == artificial_origin,
    paste0(method, " artificial doublet"),
    ifelse(real_doublet, paste0(method, "-called real doublet"), "Other real cell")
  )
  levels <- c("Other real cell", paste0(method, " artificial doublet"), paste0(method, "-called real doublet"))
  dat$display_class <- factor(dat$display_class, levels = levels)
  dat <- dat[order(dat$display_class), , drop = FALSE]

  ggplot(dat, aes(Artificial_UMAP_1, Artificial_UMAP_2)) +
    geom_point(aes(colour = display_class), size = 0.48, alpha = 0.88) +
    scale_colour_manual(values = setNames(c("#BDBDBD", artificial_colour, "#7B3294"), levels)) +
    coord_equal() +
    labs(
      title = paste0(study, ": ", method),
      subtitle = "Artificial profiles and the caller-labelled real doublets",
      x = "Artificial-doublet UMAP 1", y = "Artificial-doublet UMAP 2",
      colour = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 10)
    ) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2.4, alpha = 1)))
}

pair_rows <- list()
row_id <- 0L
for (study in study_order) {
  object_file <- file.path(object_dir, paste0(safe_name(study), "_real_and_artificial_doublet_UMAP.qs"))
  if (!file.exists(object_file)) stop("Missing study-specific artificial-doublet object: ", object_file)
  object <- qs::qread(object_file)
  for (caller_id in seq_len(nrow(callers))) {
    caller <- callers[caller_id, , drop = FALSE]
    plot <- make_method_overlay(
      object = object, study = study,
      method = caller$method,
      artificial_origin = caller$artificial_origin,
      call_column = caller$call_column,
      artificial_colour = caller$artificial_colour
    )
    stem <- file.path(figure_dir, paste0(
      "UMAP_", safe_name(study), "_", caller$method,
      "_artificial_and_called_real_doublets"
    ))
    save_plot_both(plot, stem)
    print(plot)
    row_id <- row_id + 1L
    pair_rows[[row_id]] <- data.frame(
      Study = study,
      method = caller$method,
      figure_stem = stem,
      PNG_exists = file.exists(paste0(stem, ".png")),
      PDF_exists = file.exists(paste0(stem, ".pdf")),
      stringsAsFactors = FALSE
    )
  }
  rm(object)
  invisible(gc())
}

pair_check <- do.call(rbind, pair_rows)
write.csv(pair_check, file.path(table_dir, "per_method_artificial_and_called_real_doublet_PNG_PDF_pair_check.csv"), row.names = FALSE)
if (any(!pair_check$PNG_exists | !pair_check$PDF_exists)) {
  stop("At least one per-method artificial-doublet overlay lacks a PNG/PDF pair.")
}
writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "per_method_artificial_and_called_real_doublet_R_sessionInfo.txt"))
message("Per-method artificial-doublet overlays are complete: ", figure_dir)
