library(Seurat)
library(dplyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(tibble)
library(qs)

data <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))

target_cts <- paste0("MG", 0:5)
outdir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/compare/MG_boxplots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

meta <- data@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    Sample = as.character(Sample),
    Condition = as.character(Condition),
    CellType1 = as.character(CellType1)
  )

# 横轴顺序：疾病在前，Control 在后
all_cond <- unique(meta$Condition)
disease_first <- c("AD", "PDD", "DLB", "MS", "Stroke", "COVID")
disease_present <- setdiff(all_cond, "Control")
cond_order <- c(
  intersect(disease_first, disease_present),
  setdiff(disease_present, disease_first),
  "Control"
)
cond_order <- cond_order[cond_order %in% all_cond]

for (target_ct in target_cts) {
  # 每个样本中目标细胞比例（分母=该样本在当前对象中的全部细胞）
  df_prop <- meta %>%
    group_by(Sample, Condition) %>%
    summarise(
      n_all = n(),
      n_target = sum(CellType1 == target_ct),
      prop = n_target / n_all,
      .groups = "drop"
    ) %>%
    filter(Condition %in% cond_order) %>%
    mutate(Condition = factor(Condition, levels = cond_order))
  
  # 每个疾病 vs Control（BH 校正）
  stat_df <- tryCatch({
    df_prop %>%
      wilcox_test(prop ~ Condition, ref.group = "Control", p.adjust.method = "BH") %>%
      mutate(
        label = paste0("p.adj=", signif(p.adj, 2)),
        y.position = max(df_prop$prop, na.rm = TRUE) * (1.05 + 0.1 * row_number())
      )
  }, error = function(e) NULL)
  
  p <- ggplot(df_prop, aes(x = Condition, y = prop, color = Condition)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, fill = "white") +
    geom_jitter(width = 0.15, size = 2, alpha = 0.85) +
    labs(
      x = NULL,
      y = paste0(target_ct, " proportion per sample"),
      title = paste0(target_ct, " proportion by condition")
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
  
  if (!is.null(stat_df) && nrow(stat_df) > 0) {
    p <- p + stat_pvalue_manual(
      stat_df, label = "label", tip.length = 0.01, hide.ns = FALSE
    )
  }
  
  ggsave(file.path(outdir, paste0(target_ct, "_fraction_by_condition.png")),
         plot = p, width = 6, height = 5, dpi = 300)
  ggsave(file.path(outdir, paste0(target_ct, "_fraction_by_condition.pdf")),
         plot = p, width = 6, height = 5)
}
