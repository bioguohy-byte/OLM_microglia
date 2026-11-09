#!/usr/bin/env Rscript
# Figure-specific grouping of the original human microglia plotting code.
suppressPackageStartupMessages({library(Seurat);library(qs);library(dplyr);library(ggplot2);library(forcats);library(patchwork);library(ggplotify);library(tidydr);library(RColorBrewer);library(pheatmap)})
input_obj_v2 <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs")
output_root <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "S3C", "plots")
figure_dir <- file.path(output_root, "pdf")
png_dir <- file.path(output_root, "png")
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(png_dir, recursive=TRUE, showWarnings=FALSE)
group_col <- "Study"
type_col <- "CellType1"
tmp_obj <- qread(input_obj_v2)
plot_df_num <- tmp_obj@meta.data %>%
  group_by(Group = .data[[group_col]]) %>%
  summarise(
    Count = n(),
    Mean_Genes = round(mean(nFeature_RNA, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    Label = paste0(Count, " (", Mean_Genes, ")"),
    Group = reorder(Group, Count),
    Group = forcats::fct_relevel(Group, intersect("Unknown", levels(Group)))
  )

p_stat <- ggplot(plot_df_num, aes(x = Group, y = Count, fill = Group)) +
  geom_col(color = "black", width = 0.8) +
  coord_flip() +
  geom_text(aes(label = Label), hjust = 0, nudge_y = 200, size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = NULL, y = NULL, title = "Number of nuclei (mean genes per nucleus)") +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "plain"),
    axis.text = element_text(size = 12, color = "black")
  )

ggsave(
  file.path(figure_dir, paste0("Bar_Numbers_", group_col, ".pdf")),
  p_stat,
  width = 8,
  height = 6
)
ggsave(file.path(png_dir, paste0("Bar_Numbers_", group_col, ".png")), p_stat, width=8, height=6, dpi=300)

