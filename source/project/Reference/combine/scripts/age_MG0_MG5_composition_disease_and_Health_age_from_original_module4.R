#!/usr/bin/env Rscript
# Disease-condition and Health-age MG0-MG5 composition plot.
# The counting and geom_col(position = "fill") block below is retained from
# Module 4 (lines 381-409) of the original combine/scripts/04.module.r.

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(readr)
})

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
combine_root <- file.path(project_root, "combine")
input_file <- file.path(combine_root, "results/age_MG0_MG5_ROE/mg_cca_meta_v2_with_unified_age_metadata.qs")
output_dir <- file.path(combine_root, "results/age_MG0_MG5_ROE")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Age-annotated combined object was not found: ", input_file)
tmp_obj <- qread(input_file)

required_columns <- c("CellType1", "age_group", "pathology_status")
if (!all(required_columns %in% names(tmp_obj@meta.data))) {
  stop("Required metadata are missing from the age-annotated object.")
}

# The original disease groups are retained. Health cells are represented twice
# by design: once as the original aggregate Health group and once in mutually
# exclusive age strata. This duplication is only in plot_data, never in the
# Seurat object or its underlying metadata.
plot_metadata <- tmp_obj@meta.data %>%
  mutate(CellType1 = factor(as.character(CellType1), levels = paste0("MG", 0:5))) %>%
  filter(!is.na(CellType1), !is.na(pathology_status)) %>%
  filter(pathology_status != "Health") %>%
  transmute(CellType1, Group = as.character(pathology_status)) %>%
  bind_rows(
    tmp_obj@meta.data %>%
      filter(pathology_status == "Health", !is.na(CellType1)) %>%
      transmute(CellType1 = factor(as.character(CellType1), levels = paste0("MG", 0:5)),
                Group = "Health"),
    tmp_obj@meta.data %>%
      filter(pathology_status == "Health", !is.na(CellType1), !is.na(age_group)) %>%
      transmute(CellType1 = factor(as.character(CellType1), levels = paste0("MG", 0:5)),
                Group = paste0("Health: ", as.character(age_group)))
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = c("Health: unknown", "Health: aging", "Health: middle age",
                 "Health: young", "Health", "AD", "PDD", "DLB", "MS", "Stroke", "COVID")
    )
  )

# MODULE: Single-cell Sample & Composition Analysis -- original Module 4 pattern
group_col <- "Group"
type_col <- "CellType1"
manual_colors <- c(
  "MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A",
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628"
)

# Original 04.module.r lines 381-394: single-cell count by Group x CellType,
# then normalize within every displayed group. No sample-level weighting.
plot_df_pct <- plot_metadata %>%
  filter(!is.na(.data[[group_col]]), !is.na(.data[[type_col]])) %>%
  group_by(Group = .data[[group_col]], CellType = .data[[type_col]]) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(
    Proportion = Count / sum(Count),
    Label = ifelse(Proportion > 0.03, paste0(round(Proportion * 100, 0), "%"), ""),
    Group = factor(Group, levels = levels(plot_metadata[[group_col]]))
  )

# Original 04.module.r lines 396-407: 100% stacked bars and historical palette.
p_stack <- ggplot(plot_df_pct, aes(x = Group, y = Count, fill = CellType)) +
  geom_col(position = "fill", color = "black", width = 0.8) +
  coord_flip() +
  geom_text(aes(label = Label), position = position_fill(vjust = 0.5), size = 3.5) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = manual_colors) +
  labs(
    x = NULL, y = "Proportion", fill = "Cell Type",
    title = "Cell Type Composition per Disease Condition and Health Age Group"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "right"
  )

file_stub <- "Bar_Composition_disease_and_Health_age"
# Requested canvas change: retain 5-inch height and double the width from 8 to 16 inches.
ggsave(file.path(output_dir, paste0(file_stub, ".pdf")), p_stack, width = 16, height = 5)
ggsave(file.path(output_dir, paste0(file_stub, ".png")), p_stack, width = 16, height = 5, dpi = 300)
write_csv(plot_df_pct, file.path(output_dir, paste0(file_stub, "_source_data.csv")))

message("Completed disease-condition + Health-age composition plot.")
