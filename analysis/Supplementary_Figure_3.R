#!/usr/bin/env Rscript
# Figure-ordered code from the reviewed manuscript analysis package.
# Each panel contains its original quoted R code.
# Source this chapter, then use mg_list_panels() or mg_show_panel("PANEL_ID").
# Source paths in panel records identify the original analysis scripts.
.mg_panels <- new.env(parent=emptyenv())
.mg_order <- character()
mg_register_panel <- function(id, code) {
  if (!id %in% .mg_order) .mg_order <<- c(.mg_order,id)
  assign(id, code, envir=.mg_panels)
}
mg_list_panels <- function() .mg_order
mg_panel_code <- function(id) {
  if (!exists(id,envir=.mg_panels,inherits=FALSE)) stop("Unknown panel: ",id)
  get(id,envir=.mg_panels,inherits=FALSE)
}
mg_show_panel <- function(id) {
  blocks <- mg_panel_code(id)
  if (is.null(blocks)) {cat(id, ": No R block is included for this panel.\n");return(invisible(NULL))}
  for (block in blocks) {
    cat("# Source: ",block$source," [actual ",block$ranges,"; requested ",block$requested_ranges,"]\n# Role: ",block$role,"\n",sep="")
    cat(paste(deparse(block$code,width.cutoff=120L),collapse="\n"),"\n\n")
  }
  invisible(blocks)
}


# Supplementary Figure 3

# S3A  Human MG UMAP split by Condition
# Kind: computational
mg_register_panel("S3A", list(
  list(source="source/figure_settings/S3A_Condition.R", ranges="1-51", requested_ranges="1-51", role="Original plotting expressions with the figure-specific grouping", code=quote({
#!/usr/bin/env Rscript
# Figure-specific grouping of the original human microglia plotting code.
suppressPackageStartupMessages({library(Seurat);library(qs);library(dplyr);library(ggplot2);library(forcats);library(patchwork);library(ggplotify);library(tidydr);library(RColorBrewer);library(pheatmap)})
input_obj_v2 <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs")
output_root <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "S3A", "plots")
figure_dir <- file.path(output_root, "pdf")
png_dir <- file.path(output_root, "png")
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(png_dir, recursive=TRUE, showWarnings=FALSE)
group_col <- "Condition"
type_col <- "CellType1"
tmp_obj <- qread(input_obj_v2)

p_group_dim <- DimPlot(
  tmp_obj,
  reduction = "umap",
  group.by = group_col,
  raster = FALSE,
  label = FALSE,
  pt.size = 0.1
) +
  ggtitle(group_col) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2.5))) +
  tidydr::theme_dr(
    xlength = 0.2,
    ylength = 0.2,
    arrow = grid::arrow(length = grid::unit(0.15, "inches"), type = "closed")
  ) &
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none",
    axis.title.x = element_text(hjust = 0),
    axis.title.y = element_text(hjust = 0)
  )

p_group_split <- DimPlot(tmp_obj, group.by = group_col, split.by = group_col) +
  NoAxes() +
  theme(plot.title = element_blank())

combined_group_dim <- (as.ggplot(p_group_dim) | as.ggplot(p_group_split)) +
  plot_layout(widths = c(1, 6))

ggsave(
  file.path(figure_dir, paste0("UMAP_Split_", group_col, ".pdf")),
  combined_group_dim,
  width = 18,
  height = 3
)
ggsave(file.path(png_dir, paste0("UMAP_Split_", group_col, ".png")), combined_group_dim, width=18, height=3, dpi=300)

  }))
))

# S3B  Human MG UMAP split by Region
# Kind: computational
mg_register_panel("S3B", list(
  list(source="source/figure_settings/S3B_Region.R", ranges="1-51", requested_ranges="1-51", role="Original plotting expressions with the figure-specific grouping", code=quote({
#!/usr/bin/env Rscript
# Figure-specific grouping of the original human microglia plotting code.
suppressPackageStartupMessages({library(Seurat);library(qs);library(dplyr);library(ggplot2);library(forcats);library(patchwork);library(ggplotify);library(tidydr);library(RColorBrewer);library(pheatmap)})
input_obj_v2 <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs")
output_root <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "S3B", "plots")
figure_dir <- file.path(output_root, "pdf")
png_dir <- file.path(output_root, "png")
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(png_dir, recursive=TRUE, showWarnings=FALSE)
group_col <- "Region"
type_col <- "CellType1"
tmp_obj <- qread(input_obj_v2)

p_group_dim <- DimPlot(
  tmp_obj,
  reduction = "umap",
  group.by = group_col,
  raster = FALSE,
  label = FALSE,
  pt.size = 0.1
) +
  ggtitle(group_col) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2.5))) +
  tidydr::theme_dr(
    xlength = 0.2,
    ylength = 0.2,
    arrow = grid::arrow(length = grid::unit(0.15, "inches"), type = "closed")
  ) &
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none",
    axis.title.x = element_text(hjust = 0),
    axis.title.y = element_text(hjust = 0)
  )

p_group_split <- DimPlot(tmp_obj, group.by = group_col, split.by = group_col) +
  NoAxes() +
  theme(plot.title = element_blank())

combined_group_dim <- (as.ggplot(p_group_dim) | as.ggplot(p_group_split)) +
  plot_layout(widths = c(1, 6))

ggsave(
  file.path(figure_dir, paste0("UMAP_Split_", group_col, ".pdf")),
  combined_group_dim,
  width = 18,
  height = 3
)
ggsave(file.path(png_dir, paste0("UMAP_Split_", group_col, ".png")), combined_group_dim, width=18, height=3, dpi=300)

  }))
))

# S3C  Nuclei and mean genes by Study
# Kind: computational
mg_register_panel("S3C", list(
  list(source="source/figure_settings/S3C_Study.R", ranges="1-46", requested_ranges="1-46", role="Original plotting expressions with the figure-specific grouping", code=quote({
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

  }))
))

# S3D  Nuclei and mean genes by Region
# Kind: computational
mg_register_panel("S3D", list(
  list(source="source/figure_settings/S3D_Region.R", ranges="1-46", requested_ranges="1-46", role="Original plotting expressions with the figure-specific grouping", code=quote({
#!/usr/bin/env Rscript
# Figure-specific grouping of the original human microglia plotting code.
suppressPackageStartupMessages({library(Seurat);library(qs);library(dplyr);library(ggplot2);library(forcats);library(patchwork);library(ggplotify);library(tidydr);library(RColorBrewer);library(pheatmap)})
input_obj_v2 <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs")
output_root <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "S3D", "plots")
figure_dir <- file.path(output_root, "pdf")
png_dir <- file.path(output_root, "png")
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(png_dir, recursive=TRUE, showWarnings=FALSE)
group_col <- "Region"
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

  }))
))

# S3E  Condition x MG0-MG5 ROE heatmap
# Kind: computational
mg_register_panel("S3E", list(
  list(source="source/figure_settings/S3E_Condition.R", ranges="1-63", requested_ranges="1-63", role="Original plotting expressions with the figure-specific grouping", code=quote({
#!/usr/bin/env Rscript
# Figure-specific grouping of the original human microglia plotting code.
suppressPackageStartupMessages({library(Seurat);library(qs);library(dplyr);library(ggplot2);library(forcats);library(patchwork);library(ggplotify);library(tidydr);library(RColorBrewer);library(pheatmap)})
input_obj_v2 <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs")
output_root <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "S3E", "plots")
figure_dir <- file.path(output_root, "pdf")
png_dir <- file.path(output_root, "png")
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(png_dir, recursive=TRUE, showWarnings=FALSE)
# Fig.09: ROE Heatmap ----------------------------------------------------------
var1 <- "CellType1"
var2 <- "Condition"
roe_output_path <- file.path(figure_dir, paste0("ROE_", var1, "_", var2, ".pdf"))

tmp_obj <- qread(input_obj_v2)

tab <- table(tmp_obj[[var1, drop = TRUE]], tmp_obj[[var2, drop = TRUE]])

if ("Unknown" %in% colnames(tab)) {
  tab <- tab[, colnames(tab) != "Unknown", drop = FALSE]
}

row_sum <- rowSums(tab)
col_sum <- colSums(tab)
total <- sum(tab)
expected <- outer(row_sum, col_sum) / total
roe <- as.matrix(tab / expected)

roe_plot_pal <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(100))

pdf(roe_output_path, width = 6, height = 6)
pheatmap(
  roe,
  display_numbers = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = roe_plot_pal,
  number_format = "%.1f",
  number_color = "white",
  cellwidth = 35,
  cellheight = 25,
  main = paste0("ROE: ", var1, " vs ", var2),
  angle_col = "45"
)
dev.off()

message("Fig.09 done.")
png(file.path(png_dir, paste0("ROE_", var1, "_", var2, ".png")),width=6,height=6,units="in",res=300)
pheatmap(
  roe,
  display_numbers = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = roe_plot_pal,
  number_format = "%.1f",
  number_color = "white",
  cellwidth = 35,
  cellheight = 25,
  main = paste0("ROE: ", var1, " vs ", var2),
  angle_col = "45"
)
dev.off()

  }))
))

# S3F  Region x MG0-MG5 ROE heatmap
# Kind: computational
mg_register_panel("S3F", list(
  list(source="source/figure_settings/S3F_Region.R", ranges="1-63", requested_ranges="1-63", role="Original plotting expressions with the figure-specific grouping", code=quote({
#!/usr/bin/env Rscript
# Figure-specific grouping of the original human microglia plotting code.
suppressPackageStartupMessages({library(Seurat);library(qs);library(dplyr);library(ggplot2);library(forcats);library(patchwork);library(ggplotify);library(tidydr);library(RColorBrewer);library(pheatmap)})
input_obj_v2 <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs")
output_root <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "S3F", "plots")
figure_dir <- file.path(output_root, "pdf")
png_dir <- file.path(output_root, "png")
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(png_dir, recursive=TRUE, showWarnings=FALSE)
# Fig.09: ROE Heatmap ----------------------------------------------------------
var1 <- "CellType1"
var2 <- "Region"
roe_output_path <- file.path(figure_dir, paste0("ROE_", var1, "_", var2, ".pdf"))

tmp_obj <- qread(input_obj_v2)

tab <- table(tmp_obj[[var1, drop = TRUE]], tmp_obj[[var2, drop = TRUE]])

if ("Unknown" %in% colnames(tab)) {
  tab <- tab[, colnames(tab) != "Unknown", drop = FALSE]
}

row_sum <- rowSums(tab)
col_sum <- colSums(tab)
total <- sum(tab)
expected <- outer(row_sum, col_sum) / total
roe <- as.matrix(tab / expected)

roe_plot_pal <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(100))

pdf(roe_output_path, width = 6, height = 6)
pheatmap(
  roe,
  display_numbers = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = roe_plot_pal,
  number_format = "%.1f",
  number_color = "white",
  cellwidth = 35,
  cellheight = 25,
  main = paste0("ROE: ", var1, " vs ", var2),
  angle_col = "45"
)
dev.off()

message("Fig.09 done.")
png(file.path(png_dir, paste0("ROE_", var1, "_", var2, ".png")),width=6,height=6,units="in",res=300)
pheatmap(
  roe,
  display_numbers = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = roe_plot_pal,
  number_format = "%.1f",
  number_color = "white",
  cellwidth = 35,
  cellheight = 25,
  main = paste0("ROE: ", var1, " vs ", var2),
  angle_col = "45"
)
dev.off()

  }))
))

# S3G  NFIX SOX2 SOX10 SREBF1 expression and AUCell scores
# Kind: computational
mg_register_panel("S3G", list(
  list(source="source/project/code_upload/SCENIC.r", ranges="136-192", requested_ranges="136-191", role="", code=quote({
gene_list <- c("NFIX", "SOX10", "SOX2", "SREBF1", "ZKSCAN1")
pal <- viridis(n = 10, option = "D" )

# Fig.04 ---------------
for (i in gene_list) {
  j <- paste0(i,"(+)")
  
  tmp.SeuObj <- tmp_obj
  feature <- i
  tmp.regulon <- j
  
  plot.data <- FeaturePlot(tmp.SeuObj, features = feature, order = T, cols = pal, raster = F)$data
  fdata <- FetchData(object = tmp.SeuObj, vars = feature)
  plot.data$feature <- fdata[match(rownames(plot.data), rownames(fdata)), ]
  
  plot1 <- ggplot(plot.data, aes(UMAP_1, UMAP_2, color = feature)) + 
    geom_point(size = 0.001) +
    scale_color_gradientn(colors = pal) +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) + 
    NoAxes() +
    ggtitle(paste0(gsub("\\(\\+\\)", "", tmp.regulon), " gene expression"))
  
  tmp.df <- data.frame(tmp.regulon = sub_regulonAUC@assays@data@listData[["AUC"]][tmp.regulon, ]) 
  tmp.SeuObj@meta.data$tmp.regulon <- tmp.df$tmp.regulon
  
  plot.data_auc <- FeaturePlot(tmp.SeuObj, features = "tmp.regulon", order = T, cols = pal, raster = F)$data
  plot.data_auc$tmp.regulon <- tmp.df[match(rownames(plot.data_auc), rownames(tmp.df)), ]
  
  plot2 <- ggplot(plot.data_auc, aes(UMAP_1, UMAP_2, color = tmp.regulon)) +
    geom_point(size = 0.001) +
    scale_color_gradientn(colors = pal) +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) + 
    NoAxes() +
    ggtitle(paste0("AUCell Score of the ", gsub("\\(\\+\\)", "", tmp.regulon), " regulon"))
  
  plot_final <- plot1 + plot2
  
  filename_pdf <- paste0(base_path, "Fig.4_AUCell_", gsub("\\(\\+\\)", "", tmp.regulon), ".v2.pdf")
  ggsave(filename = filename_pdf, plot_final, width = 9, height = 3.5) 
  
  filename_png <- paste0(base_path, "Fig.4_AUCell_", gsub("\\(\\+\\)", "", tmp.regulon), ".v2.png")
  ggsave(filename = filename_png, plot_final, width = 9, height = 3.5) 
}
  }))
))
