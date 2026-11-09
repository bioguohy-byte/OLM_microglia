#!/usr/bin/env Rscript
# Age-stratified MG0-MG5 ROE analysis.
# The ROE calculation and pheatmap rendering below are copied from Module 6
# of Reference/combine/scripts/04.module.r, with var2 changed to age_group.

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(readr)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
combine_root <- file.path(project_root, "combine")
input_object <- file.path(combine_root, "data/processed_data/mg_cca_meta_v2.qs")
input_age_table <- file.path(combine_root, "results/table/original_six_study_sample_age.csv")
output_dir <- file.path(combine_root, "results/age_MG0_MG5_ROE")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_object)) stop("Combined Seurat object was not found: ", input_object)
if (!file.exists(input_age_table)) stop("Original sample-age table was not found: ", input_age_table)

# Supplementary Table 1 (Grubman et al., 2019): paired case labels identify the
# two donors in the listed order, e.g. Ct1_Ct2 donor 1 -> Ct1 and donor 2 -> Ct2.
grubman_age_map <- tribble(
  ~Study, ~Sample, ~age_numeric, ~supplement_case, ~age_source,
  "Grubman et al., 2019", "AD1", 91.0, "AD1_AD2 donor 1", "Grubman 2019 Supplementary Table 1",
  "Grubman et al., 2019", "AD2", 83.8, "AD1_AD2 donor 2", "Grubman 2019 Supplementary Table 1",
  "Grubman et al., 2019", "AD3", 67.8, "AD3_AD4 donor 1", "Grubman 2019 Supplementary Table 1",
  "Grubman et al., 2019", "AD4", 83.0, "AD3_AD4 donor 2", "Grubman 2019 Supplementary Table 1",
  "Grubman et al., 2019", "AD5", 73.0, "AD5_AD6 donor 1", "Grubman 2019 Supplementary Table 1",
  "Grubman et al., 2019", "AD6", 74.6, "AD5_AD6 donor 2", "Grubman 2019 Supplementary Table 1",
  "Grubman et al., 2019", "Ct1", 67.3, "Ct1_Ct2 donor 1", "Grubman 2019 Supplementary Table 1",
  "Grubman et al., 2019", "Ct2", 82.7, "Ct1_Ct2 donor 2", "Grubman 2019 Supplementary Table 1",
  "Grubman et al., 2019", "Ct3", 72.6, "Ct3_Ct4 donor 1", "Grubman 2019 Supplementary Table 1",
  "Grubman et al., 2019", "Ct4", 75.6, "Ct3_Ct4 donor 2", "Grubman 2019 Supplementary Table 1",
  "Grubman et al., 2019", "Ct5", 77.5, "Ct5_Ct6 donor 1", "Grubman 2019 Supplementary Table 1",
  "Grubman et al., 2019", "Ct6", 82.7, "Ct5_Ct6 donor 2", "Grubman 2019 Supplementary Table 1"
)

write_csv(grubman_age_map, file.path(output_dir, "Grubman_Supplementary_Table1_age_mapping.csv"))

age_input <- read_csv(input_age_table, show_col_types = FALSE, na = c("", "NA", "Unknown"))
required_age_columns <- c("Study", "Sample", "Age", "Condition_source", "Pathology_status")
if (!all(required_age_columns %in% names(age_input))) {
  stop("Unexpected input age-table columns. Expected: ", paste(required_age_columns, collapse = ", "))
}
if (anyDuplicated(paste(age_input$Study, age_input$Sample, sep = "||"))) {
  stop("The input age table has duplicate Study/Sample keys.")
}

age_table <- age_input %>%
  mutate(
    Study = as.character(Study),
    Sample = as.character(Sample),
    age_source_text = as.character(Age),
    existing_age_numeric = suppressWarnings(parse_number(age_source_text))
  ) %>%
  select(-Age) %>%
  left_join(grubman_age_map, by = c("Study", "Sample")) %>%
  mutate(
    age_numeric = if_else(!is.na(age_numeric), age_numeric, existing_age_numeric),
    age_source = case_when(
      !is.na(supplement_case) ~ age_source,
      !is.na(age_numeric) ~ "Existing project age audit",
      !is.na(age_source_text) ~ "Existing qualitative/unparsed age label",
      TRUE ~ "Age not available in project age audit"
    ),
    # The user-specified bins leave no samples at 28-29. We include age 27 in
    # young so the known 27-year-old Su sample is not reclassified as unknown.
    age_group = case_when(
      !is.na(age_numeric) & age_numeric <= 27 ~ "young",
      !is.na(age_numeric) & age_numeric >= 30 & age_numeric < 60 ~ "middle age",
      !is.na(age_numeric) & age_numeric >= 60 ~ "aging",
      TRUE ~ "unknown"
    ),
    age_group_reason = case_when(
      is.na(age_numeric) ~ "numeric age unavailable",
      age_numeric <= 27 ~ "age <= 27",
      age_numeric >= 30 & age_numeric < 60 ~ "30 <= age < 60",
      age_numeric >= 60 ~ "age >= 60",
      TRUE ~ "known age outside prespecified bins (28-29); assigned unknown"
    ),
    age_group = factor(age_group, levels = c("young", "middle age", "aging", "unknown"))
  ) %>%
  select(Study, Sample, Condition_source, Pathology_status, age_source_text,
         age_numeric, age_group, age_group_reason, age_source, supplement_case)

write_csv(age_table, file.path(output_dir, "all_studies_sample_age_unified.csv"), na = "")

obj <- qread(input_object)
metadata <- obj@meta.data
study_col <- c("Study", "study")[c("Study", "study") %in% names(metadata)][1]
sample_col <- c("Sample", "sample", "Patient", "patient")[c("Sample", "sample", "Patient", "patient") %in% names(metadata)][1]
if (is.na(study_col) || is.na(sample_col)) {
  stop("Could not identify Study/Sample columns in object metadata. Available columns: ",
       paste(names(metadata), collapse = ", "))
}
if (!"CellType1" %in% names(metadata)) stop("CellType1 is absent from combined-object metadata.")

object_keys <- paste(as.character(metadata[[study_col]]), as.character(metadata[[sample_col]]), sep = "||")
age_keys <- paste(age_table$Study, age_table$Sample, sep = "||")
age_index <- match(object_keys, age_keys)

metadata$age_numeric <- age_table$age_numeric[age_index]
metadata$age_group <- as.character(age_table$age_group[age_index])
metadata$age_group_reason <- age_table$age_group_reason[age_index]
metadata$age_source <- age_table$age_source[age_index]
metadata$pathology_status <- age_table$Pathology_status[age_index]
metadata$condition_source_age_audit <- age_table$Condition_source[age_index]

mapping_audit <- metadata %>%
  as.data.frame() %>%
  transmute(
    Study = as.character(.data[[study_col]]),
    Sample = as.character(.data[[sample_col]]),
    pathology_status, age_numeric, age_group, age_group_reason, age_source
  ) %>%
  distinct() %>%
  arrange(Study, Sample)
write_csv(mapping_audit, file.path(output_dir, "combined_object_sample_age_mapping_audit.csv"), na = "")

unmapped <- mapping_audit %>% filter(is.na(pathology_status))
if (nrow(unmapped) > 0) {
  write_csv(unmapped, file.path(output_dir, "unmapped_combined_object_samples.csv"), na = "")
  stop("Some object Study/Sample keys are absent from the age audit; see unmapped_combined_object_samples.csv")
}

obj@meta.data <- metadata
qsave(obj, file.path(output_dir, "mg_cca_meta_v2_with_unified_age_metadata.qs"))

# Only pathology_status == Health enters the requested age comparison.
health_metadata <- obj@meta.data %>%
  as.data.frame() %>%
  filter(pathology_status == "Health", CellType1 %in% paste0("MG", 0:5)) %>%
  mutate(
    CellType1 = factor(as.character(CellType1), levels = paste0("MG", 0:5)),
    age_group = factor(age_group, levels = c("young", "middle age", "aging", "unknown"))
  )
if (nrow(health_metadata) == 0) stop("No Health MG0-MG5 cells were retained after age mapping.")

health_sample_audit <- health_metadata %>%
  transmute(
    Study = as.character(.data[[study_col]]),
    Sample = as.character(.data[[sample_col]]),
    age_numeric, age_group, age_group_reason, age_source
  ) %>%
  distinct() %>%
  arrange(age_group, Study, Sample)
write_csv(health_sample_audit, file.path(output_dir, "Health_samples_in_age_ROE.csv"), na = "")

health_counts <- health_metadata %>%
  count(age_group, CellType1, name = "observed_count") %>%
  tidyr::complete(age_group = factor(c("young", "middle age", "aging", "unknown"),
                                     levels = c("young", "middle age", "aging", "unknown")),
                  CellType1 = factor(paste0("MG", 0:5), levels = paste0("MG", 0:5)),
                  fill = list(observed_count = 0))
write_csv(health_counts, file.path(output_dir, "Health_age_group_by_MG0_MG5_observed_counts.csv"))

# Module 6 of the original 04.module.r: table -> expected under independence -> ROE -> pheatmap.
plot_roe <- function(plot_metadata, suffix) {
  var1 <- "CellType1"
  var2 <- "age_group"
  tab <- table(plot_metadata[[var1]], plot_metadata[[var2]])
  row_sum <- rowSums(tab)
  col_sum <- colSums(tab)
  total <- sum(tab)
  expected <- outer(row_sum, col_sum) / total
  roe <- as.matrix(tab / expected)
  rows_to_plot <- rownames(roe)
  cols_to_plot <- colnames(roe)
  plot_pal <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(100))

  output_pdf <- file.path(output_dir, paste0("ROE_CellType1_age_group_", suffix, ".pdf"))
  output_png <- file.path(output_dir, paste0("ROE_CellType1_age_group_", suffix, ".png"))
  output_csv <- file.path(output_dir, paste0("ROE_CellType1_age_group_", suffix, "_matrix.csv"))
  roe_export <- data.frame(CellType1 = rownames(roe), unclass(roe), check.names = FALSE)
  write_csv(roe_export, output_csv)

  draw_plot <- function() {
    pheatmap(
      roe[rows_to_plot, cols_to_plot],
      display_numbers = TRUE,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = plot_pal,
      number_format = "%.1f",
      number_color = "white",
      cellwidth = 35,
      cellheight = 25,
      main = paste0("ROE: ", var1, " vs ", var2, " (Health only)"),
      angle_col = "45"
    )
  }
  pdf(output_pdf, width = 6, height = 6)
  draw_plot()
  dev.off()
  png(output_png, width = 1800, height = 1800, res = 300)
  draw_plot()
  dev.off()
  invisible(roe)
}

# Requested complete Health analysis: unknown is retained as an explicitly labelled
# non-age-resolved group. A second matrix removes only that non-biological group.
roe_all_health <- plot_roe(health_metadata, "all_Health_including_unknown")
roe_known_age <- plot_roe(
  filter(health_metadata, age_group != "unknown") %>% mutate(age_group = droplevels(age_group)),
  "known_age_Health_only"
)

manifest <- c(
  "Input object: combine/data/processed_data/mg_cca_meta_v2.qs",
  "Age audit input: combine/results/table/original_six_study_sample_age.csv",
  "Grubman update: Supplementary Table 1, paired donor order mapped to AD1-AD6 and Ct1-Ct6",
  "Metadata added: age_numeric (numeric), age_group, age_group_reason, age_source, pathology_status",
  "Health filter: pathology_status == Health",
  "Age bins used: young <=27; middle age 30<=age<60; aging >=60; unknown when numeric age unavailable.",
  "ROE: observed count / expected count under independence, copied from Module 6 of combine/scripts/04.module.r."
)
writeLines(manifest, file.path(output_dir, "manifest.txt"))
message("Completed age metadata annotation and Health-only age ROE analysis.")
message("Output directory: ", output_dir)
