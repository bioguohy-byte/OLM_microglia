#!/usr/bin/env Rscript

# Extract one auditable age record per sample from the six original Seurat
# objects that were inputs to the historical MG0-MG5 integration.
suppressPackageStartupMessages(library(qs))

reference_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
output_file <- file.path(reference_root, "combine/results/table/original_six_study_sample_age.csv")

studies <- data.frame(
  Study = c(
    "Su et al., 2022",
    "Chia et al., 2025",
    "Fullard et al., 2021",
    "Jakel et al., 2019",
    "Grubman et al., 2019",
    "Yao et al., 2024"
  ),
  object_file = c(
    "Su_et_al./data/processed_data/Su_mg.qs",
    "Chia_et_al./data/processed_data/Chia_mg.qs",
    "Fullard_et_al./data/processed_data/Fullard_mg.qs",
    "Jakel_et_al./data/processed_data/Jakel_mg.qs",
    "Grubman_et_al./data/processed_data/Grubman_mg.qs",
    "Yao_et_al./data/processed_data/Yao_mg.qs"
  ),
  stringsAsFactors = FALSE
)

extract_sample_age <- function(study, object_file) {
  object <- qs::qread(file.path(reference_root, object_file))
  metadata <- object@meta.data
  required_fields <- c("Sample", "Age", "Condition")
  missing_fields <- setdiff(required_fields, names(metadata))
  if (length(missing_fields)) {
    stop(study, ": missing metadata field(s): ", paste(missing_fields, collapse = ", "))
  }

  # The status label is derived directly from the original input object's
  # Condition field. It follows the seven labels shown in the historical
  # cell-type-proportion panel, rather than relabelling the later harmonized
  # "Control" field. Condition_source preserves the exact input value.
  condition_source <- as.character(metadata$Condition)
  unexpected_conditions <- setdiff(
    unique(condition_source),
    c("CTRL", "Ctrl", "ct", "Healthy", "Control", "AD", "ADD", "PDD", "DLB", "MS", "COVID", "COVID_HIV", "Stroke")
  )
  if (length(unexpected_conditions)) {
    stop(study, ": unexpected original Condition value(s): ",
         paste(unexpected_conditions, collapse = ", "))
  }
  pathology_status <- dplyr::case_when(
    condition_source %in% c("CTRL", "Ctrl", "ct", "Healthy", "Control") ~ "Health",
    condition_source %in% c("AD", "ADD") ~ "AD",
    condition_source == "PDD" ~ "PDD",
    condition_source == "DLB" ~ "DLB",
    condition_source == "MS" ~ "MS",
    condition_source %in% c("COVID", "COVID_HIV") ~ "COVID",
    condition_source == "Stroke" ~ "Stroke"
  )
  sample_age <- unique(data.frame(
    Study = study,
    Sample = as.character(metadata$Sample),
    Age = as.character(metadata$Age),
    Condition_source = condition_source,
    Pathology_status = pathology_status,
    stringsAsFactors = FALSE
  ))
  if (anyNA(sample_age$Sample) || any(sample_age$Sample == "")) {
    stop(study, ": empty Sample values are not permitted in the export.")
  }
  age_conflicts <- aggregate(Age ~ Sample, data = sample_age, FUN = function(x) length(unique(x)))
  pathology_conflicts <- aggregate(
    Pathology_status ~ Sample,
    data = sample_age,
    FUN = function(x) length(unique(x))
  )
  source_condition_conflicts <- aggregate(
    Condition_source ~ Sample,
    data = sample_age,
    FUN = function(x) length(unique(x))
  )
  conflicting_samples <- unique(c(
    age_conflicts[age_conflicts$Age > 1L, "Sample"],
    pathology_conflicts[pathology_conflicts$Pathology_status > 1L, "Sample"],
    source_condition_conflicts[source_condition_conflicts$Condition_source > 1L, "Sample"]
  ))
  if (length(conflicting_samples)) {
    stop(study, ": conflicting metadata values for sample(s): ", paste(conflicting_samples, collapse = ", "))
  }
  sample_age[order(sample_age$Sample), c("Study", "Sample", "Age", "Condition_source", "Pathology_status")]
}

sample_age_table <- do.call(
  rbind,
  Map(extract_sample_age, studies$Study, studies$object_file)
)
sample_age_table$Study <- factor(sample_age_table$Study, levels = studies$Study)
sample_age_table <- sample_age_table[order(sample_age_table$Study, sample_age_table$Sample), ]
sample_age_table$Study <- as.character(sample_age_table$Study)

if (nrow(sample_age_table) != 93L || anyDuplicated(sample_age_table[c("Study", "Sample")])) {
  stop("Unexpected export cardinality or duplicated Study/Sample record.")
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(sample_age_table, output_file, row.names = FALSE, na = "")
message("Wrote ", nrow(sample_age_table), " unique Study/Sample age records to: ", output_file)
