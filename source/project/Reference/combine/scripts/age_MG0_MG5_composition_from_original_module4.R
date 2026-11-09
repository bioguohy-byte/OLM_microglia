#!/usr/bin/env Rscript
# Health-only age composition plot. The calculation and geom_col(position = "fill")
# block is retained from Module 4 (lines 381-409) of combine/scripts/04.module.r.

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

if (!all(c("CellType1", "age_group", "pathology_status") %in% names(tmp_obj@meta.data))) {
  stop("Required metadata are missing from age-annotated object.")
}

# Requested analysis universe: only Health cells in the historical MG0-MG5 object.
tmp_obj <- subset(tmp_obj, subset = pathology_status == "Health")
tmp_obj@meta.data <- tmp_obj@meta.data %>%
  mutate(
    CellType1 = factor(as.character(CellType1), levels = paste0("MG", 0:5)),
    age_group = factor(as.character(age_group),
                       levels = c("young", "middle age", "aging", "unknown"))
  )

# MODULE: Single-cell Sample & Composition Analysis -- original Module 4 pattern
group_col <- "age_group"
type_col <- "CellType1"
manual_colors <- c(
  "MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A",
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628"
)

# Original 04.module.r lines 381-394: count cells by Group x CellType and calculate
# within-group percentages. No sample-level weighting is introduced.
plot_df_pct <- tmp_obj@meta.data %>%
  filter(!is.na(.data[[group_col]]), !is.na(.data[[type_col]])) %>%
  group_by(Group = .data[[group_col]], CellType = .data[[type_col]]) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(
    Proportion = Count / sum(Count),
    Label = ifelse(Proportion > 0.03, paste0(round(Proportion * 100, 0), "%"), ""),
    Group = factor(Group, levels = levels(tmp_obj@meta.data[[group_col]])),
    Group = forcats::fct_relevel(Group, intersect("unknown", levels(Group)))
  )

# Original 04.module.r lines 396-407: 100% stacked bars with the historical MG0-MG5 palette.
p_stack <- ggplot(plot_df_pct, aes(x = Group, y = Count, fill = CellType)) +
  geom_col(position = "fill", color = "black", width = 0.8) +
  coord_flip() +
  geom_text(aes(label = Label), position = position_fill(vjust = 0.5), size = 3.5) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = manual_colors) +
  labs(x = NULL, y = "Proportion", fill = "Cell Type",
       title = "Cell Type Composition per Age Group (Health only)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "right")

pdf_file <- file.path(output_dir, "Bar_Composition_age_group_Health.pdf")
png_file <- file.path(output_dir, "Bar_Composition_age_group_Health.png")
data_file <- file.path(output_dir, "Bar_Composition_age_group_Health_source_data.csv")
ggsave(pdf_file, p_stack, width = 8, height = 3)
ggsave(png_file, p_stack, width = 8, height = 3, dpi = 300)
write_csv(plot_df_pct, data_file)
message("Completed Health-only age-group composition plot: ", pdf_file)
