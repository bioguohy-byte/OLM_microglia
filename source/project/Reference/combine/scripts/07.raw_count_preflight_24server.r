#!/usr/bin/env Rscript

# Full reanalysis preflight: confirm that each study has a full-cell RNA count
# matrix with integer values and identify an available sample/batch field.
# This script is read-only with respect to source objects.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(qs)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
result_root <- file.path(project_root, "combine/results/full_count_reanalysis")
figure_dir <- file.path(result_root, "00_preflight/figure")
table_dir <- file.path(result_root, "00_preflight/table")
runtime_dir <- file.path(result_root, "00_preflight/runtime")
for (d in c(figure_dir, table_dir, runtime_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

input_manifest <- data.frame(
  Study = c("Su et al., 2022", "Chia et al., 2025", "Fullard et al., 2021",
            "Jakel et al., 2019", "Grubman et al., 2019", "Yao et al., 2024"),
  Source_object = c(
    file.path(project_root, "Su_et_al./data/processed_data/Su_CSC_alldata_v2.qs"),
    file.path(project_root, "Chia_et_al./data/processed_data/Chia.qs"),
    file.path(project_root, "Fullard_et_al./data/processed_data/Fullard.qs"),
    file.path(project_root, "Jakel_et_al./data/processed_data/Jakel.qs"),
    file.path(project_root, "Grubman_et_al./data/processed_data/Grubman.qs"),
    file.path(project_root, "Yao_et_al./data/processed_data/Yao_stroke.qs")
  ),
  Batch_field = c("Sample", "sample", "subject_id", "sample", "patient", "Batch"),
  stringsAsFactors = FALSE
)

save_plot_both <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 320,
         bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
         device = grDevices::pdf, bg = "white", limitsize = FALSE)
}

extract_integer_counts <- function(obj) {
  assays <- Seurat::Assays(obj)
  if (!"RNA" %in% assays) stop("No RNA assay. Available assays: ", paste(assays, collapse = ", "))
  assay <- obj[["RNA"]]
  counts <- tryCatch(assay@counts, error = function(e) NULL)
  if (is.null(counts) || !length(counts)) {
    counts <- tryCatch(Seurat::GetAssayData(obj, assay = "RNA", slot = "counts"), error = function(e) NULL)
  }
  if (is.null(counts) || !length(counts)) stop("RNA counts slot/layer is unavailable.")
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  counts
}

choose_sample_field <- function(metadata) {
  candidates <- c("Sample", "sample", "Sample_ID", "sample_id", "orig.ident", "Donor", "donor",
                  "donor_id", "individual", "Individual", "patient", "Patient", "Batch", "batch")
  found <- candidates[candidates %in% colnames(metadata)]
  if (length(found)) found[[1]] else NA_character_
}

overview_rows <- list()
sample_rows <- list()
field_rows <- list()
for (i in seq_len(nrow(input_manifest))) {
  study <- input_manifest$Study[[i]]
  source_file <- input_manifest$Source_object[[i]]
  requested_batch_field <- input_manifest$Batch_field[[i]]
  message("Preflight: ", study)
  if (!file.exists(source_file)) {
    overview_rows[[study]] <- data.frame(Study = study, Source_object = source_file, readable = FALSE,
      integer_RNA_counts = FALSE, n_genes = NA_integer_, n_cells = NA_integer_, sample_field = NA_character_,
      n_samples = NA_integer_, max_fractional_count_deviation = NA_real_, message = "Source file does not exist.")
    next
  }
  result <- tryCatch({
    obj <- qs::qread(source_file)
    counts <- extract_integer_counts(obj)
    metadata <- obj[[]]
    if (is.null(colnames(counts)) || anyDuplicated(colnames(counts))) stop("Count-matrix cell barcodes are missing or duplicated.")
    if (!all(colnames(counts) %in% rownames(metadata))) stop("Some RNA count cell barcodes are absent from metadata.")
    metadata <- metadata[colnames(counts), , drop = FALSE]
    sample_field <- if (requested_batch_field %in% colnames(metadata)) requested_batch_field else choose_sample_field(metadata)
    sample_values <- if (is.na(sample_field)) rep("study_level", ncol(counts)) else as.character(metadata[[sample_field]])
    sample_values[is.na(sample_values) | sample_values == ""] <- "missing_sample"
    deviation <- if (length(counts@x)) max(abs(counts@x - round(counts@x))) else 0
    overview <- data.frame(Study = study, Source_object = source_file, readable = TRUE,
      integer_RNA_counts = isTRUE(deviation < 1e-8), n_genes = nrow(counts), n_cells = ncol(counts),
      sample_field = ifelse(is.na(sample_field), "study_level_fallback", sample_field),
      n_samples = length(unique(sample_values)), max_fractional_count_deviation = deviation,
      message = if (is.na(sample_field)) "No sample field: study-level batch fallback required." else "OK")
    per_sample <- as.data.frame(table(Study = rep(study, length(sample_values)), sample_id = sample_values), stringsAsFactors = FALSE)
    colnames(per_sample)[3] <- "n_cells"
    fields <- data.frame(Study = study, metadata_field = colnames(metadata), stringsAsFactors = FALSE)
    rm(obj, counts, metadata)
    gc(verbose = FALSE)
    list(overview = overview, per_sample = per_sample, fields = fields)
  }, error = function(e) list(overview = data.frame(Study = study, Source_object = source_file, readable = FALSE,
    integer_RNA_counts = FALSE, n_genes = NA_integer_, n_cells = NA_integer_, sample_field = NA_character_,
    n_samples = NA_integer_, max_fractional_count_deviation = NA_real_, message = conditionMessage(e)),
    per_sample = NULL, fields = NULL))
  overview_rows[[study]] <- result$overview
  if (!is.null(result$per_sample)) sample_rows[[study]] <- result$per_sample
  if (!is.null(result$fields)) field_rows[[study]] <- result$fields
}

overview <- do.call(rbind, overview_rows)
per_sample <- do.call(rbind, sample_rows)
metadata_fields <- do.call(rbind, field_rows)
utils::write.csv(input_manifest, file.path(table_dir, "input_object_manifest.csv"), row.names = FALSE)
utils::write.csv(overview, file.path(table_dir, "integer_count_and_batch_preflight.csv"), row.names = FALSE)
utils::write.csv(per_sample, file.path(table_dir, "per_study_sample_cell_counts.csv"), row.names = FALSE)
utils::write.csv(metadata_fields, file.path(table_dir, "available_metadata_fields.csv"), row.names = FALSE)

plot_data <- overview[overview$readable, , drop = FALSE]
plot_data$Study <- factor(plot_data$Study, levels = input_manifest$Study)
p_preflight <- ggplot(plot_data, aes(Study, n_cells, fill = integer_RNA_counts)) +
  geom_col(width = 0.72, colour = "white") +
  geom_text(aes(label = paste0("cells=", format(n_cells, big.mark = ","), "\n", sample_field, ": ", n_samples)),
            hjust = -0.05, size = 3) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.38))) +
  scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02"), name = "Integer RNA counts") +
  labs(title = "Full-count input preflight for six studies", x = NULL, y = "Cells in source object",
       subtitle = "Labels report the selected sample/batch field and number of levels") +
  theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold"))
save_plot_both(p_preflight, file.path(figure_dir, "full_count_input_preflight"), 10, 5.5)
print(p_preflight)

writeLines(capture.output(sessionInfo()), file.path(runtime_dir, "preflight_R_sessionInfo.txt"))
if (any(!overview$readable) || any(!overview$integer_RNA_counts)) {
  message("Preflight completed with blocked inputs. Review integer_count_and_batch_preflight.csv before stage 1.")
} else {
  message("Preflight passed: all six source objects expose integer RNA counts. Stage 1 may proceed.")
}
