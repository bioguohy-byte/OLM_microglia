#!/usr/bin/env Rscript

# Run this script from the project root so relative paths resolve correctly.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(qs)
  library(openxlsx)
  library(pheatmap)
  library(ggplotify)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(patchwork)
  library(stringr)
  library(RColorBrewer)
  library(scales)
  library(forcats)
  library(cowplot)
  library(GseaVis)
  library(tidydr)
})

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
}

figure_dir <- file.path("results", "figure")
table_dir <- file.path("results", "table")
saved_dir <- file.path("data", "saved_data")
processed_dir <- file.path("data", "processed_data")
gsea_data_dir <- file.path("scripts", "GSEA", "data")
gsea_result_dir <- file.path("scripts", "GSEA", "results")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(saved_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_result_dir, recursive = TRUE, showWarnings = FALSE)

input_obj_raw <- file.path(processed_dir, "mg_cca_meta.qs")
input_obj_v2 <- file.path(processed_dir, "mg_cca_meta_v2.qs")

manual_colors <- c(
  "MG2" = "#377EB8",
  "MG5" = "#8B0000",
  "MG0" = "#4DAF4A",
  "MG3" = "#984EA3",
  "MG4" = "#E6AB02",
  "MG1" = "#A65628"
)

stop_if_missing(input_obj_raw)

# Fig.01: Prepare CellType1 labels ---------------------------------------------
tmp_obj <- qread(input_obj_raw)

tmp_obj@meta.data <- tmp_obj@meta.data %>%
  mutate(
    CellType1 = case_when(
      integrated_snn_res.0.3 %in% c("4", "6") ~ "MG0",
      integrated_snn_res.0.3 %in% c("8", "9") ~ "MG1",
      integrated_snn_res.0.3 %in% c("0", "2") ~ "MG2",
      integrated_snn_res.0.3 %in% c("5") ~ "MG3",
      integrated_snn_res.0.3 %in% c("1", "7") ~ "MG4",
      integrated_snn_res.0.3 %in% c("3") ~ "MG5",
      TRUE ~ integrated_snn_res.0.3
    )
  )

tmp_obj$CellType1 <- as.factor(tmp_obj$CellType1)
qsave(tmp_obj, input_obj_v2)

# Fig.02: CellType1 DimPlot ----------------------------------------------------
cluster_col <- "CellType1"

p_dim_label <- DimPlot(
  tmp_obj,
  group.by = cluster_col,
  cols = manual_colors,
  label = TRUE
)

p_dim <- DimPlot(
  tmp_obj,
  group.by = cluster_col,
  cols = manual_colors
)

ggsave(
  file.path(figure_dir, paste0(cluster_col, "_dimplot.pdf")),
  plot = p_dim_label,
  width = 6,
  height = 5
)

ggsave(
  file.path(figure_dir, paste0(cluster_col, "_dimplot.png")),
  plot = p_dim_label,
  width = 6,
  height = 5
)

message("Fig.02 done.")

# Fig.03: FindAllMarkers -------------------------------------------------------
cluster_col <- "CellType1"
min_pct <- 0.25
logfc_thr <- 0.25

res_suffix <- gsub("integrated_snn_", "", cluster_col)
marker_qs_path <- file.path(saved_dir, paste0(res_suffix, "_markers.qs"))
marker_xlsx_path <- file.path(table_dir, paste0(res_suffix, "_markers.xlsx"))

tmp_obj <- qread(input_obj_v2)
Idents(tmp_obj) <- tmp_obj[[cluster_col, drop = TRUE]]

message("Calculating markers for: ", cluster_col)
all_markers <- FindAllMarkers(
  tmp_obj,
  only.pos = TRUE,
  min.pct = min_pct,
  logfc.threshold = logfc_thr
)

qsave(all_markers, marker_qs_path)
write.xlsx(all_markers, file = marker_xlsx_path, overwrite = TRUE)

message("Fig.03 done.")

# Fig.04: Marker DotPlot -------------------------------------------------------
cluster_col <- "CellType1"
top_n <- 5

res_suffix <- gsub("integrated_snn_", "", cluster_col)
marker_path <- file.path(saved_dir, paste0(res_suffix, "_markers.qs"))
dotplot_png_path <- file.path(figure_dir, paste0(res_suffix, "_top", top_n, "_dotplot.png"))
dotplot_pdf_path <- file.path(figure_dir, paste0(res_suffix, "_top", top_n, "_dotplot.pdf"))

tmp_obj <- qread(input_obj_v2)
stop_if_missing(marker_path)
markers <- qread(marker_path)

top_markers <- markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup()

genes_list <- unique(top_markers$gene)

DefaultAssay(tmp_obj) <- "RNA"
p_base <- DotPlot(tmp_obj, features = genes_list, group.by = cluster_col)

p_dot <- p_base$data %>%
  rename(
    pct_exp = pct.exp,
    avg_exp = avg.exp.scaled,
    Gene = features.plot,
    Cluster = id
  ) %>%
  ggplot(aes(x = Cluster, y = Gene, size = pct_exp, color = avg_exp)) +
  geom_point() +
  scale_color_gradient2(
    low = "#000000",
    mid = "#007965",
    high = "#ffcc29",
    midpoint = 0
  ) +
  scale_y_discrete(limits = rev(genes_list)) +
  theme_classic() +
  labs(
    title = paste0("Top ", top_n, " Genes in ", res_suffix),
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.margin = margin(t = 50, r = 10, b = 10, l = 10)
  )

ggsave(dotplot_png_path, p_dot, width = 3.5, height = 9, limitsize = FALSE)
ggsave(dotplot_pdf_path, p_dot, width = 3.5, height = 9, limitsize = FALSE)

message("Fig.04 done.")

# Fig.05: Marker Heatmap -------------------------------------------------------
cluster_col <- "CellType1"
top_n <- 50

res_suffix <- gsub("integrated_snn_", "", cluster_col)
marker_path <- file.path(saved_dir, paste0(res_suffix, "_markers.qs"))
heatmap_png_path <- file.path(figure_dir, paste0(res_suffix, "_top", top_n, "_heatmap.png"))
heatmap_pdf_path <- file.path(figure_dir, paste0(res_suffix, "_top", top_n, "_heatmap.pdf"))

tmp_obj <- qread(input_obj_v2)
stop_if_missing(marker_path)
markers <- qread(marker_path)

top_markers <- markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup()

DefaultAssay(tmp_obj) <- "RNA"
avg_exp_obj <- AverageExpression(tmp_obj, group.by = cluster_col, assays = "RNA")
heatmap_data <- avg_exp_obj$RNA[top_markers$gene, ]

p_heat_obj <- pheatmap(
  heatmap_data,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  show_rownames = FALSE,
  legend = TRUE,
  silent = TRUE,
  angle_col = "0",
  border_color = NA
)

p_heat_gg <- as.ggplot(p_heat_obj)

ggsave(heatmap_png_path, p_heat_gg, width = 5, height = 18, limitsize = FALSE)
ggsave(heatmap_pdf_path, p_heat_gg, width = 5, height = 18, limitsize = FALSE)

message("Fig.05 done.")

# Fig.06: GO Enrichment --------------------------------------------------------
cluster_col <- "CellType1"
top_n_go <- 10

res_suffix <- gsub("integrated_snn_", "", cluster_col)
marker_path <- file.path(saved_dir, paste0(res_suffix, "_markers.qs"))
go_save_path <- file.path(saved_dir, paste0(res_suffix, "_res_go_list.qs"))
go_excel_path <- file.path(table_dir, paste0(res_suffix, "_GO_enrichment.xlsx"))
go_combined_pdf <- file.path(figure_dir, paste0(res_suffix, "_GO_Combined_Text.pdf"))
go_barplot_dir <- file.path(figure_dir, paste0(res_suffix, "_individual_GO_barplots"))

dir.create(go_barplot_dir, recursive = TRUE, showWarnings = FALSE)
stop_if_missing(marker_path)
markers <- qread(marker_path)

all_clusters <- as.character(sort(unique(markers$cluster)))
cluster_colors <- manual_colors[all_clusters]
missing_clusters <- all_clusters[is.na(cluster_colors)]

if (length(missing_clusters) > 0) {
  auto_cols <- if (length(missing_clusters) <= 12) {
    brewer.pal(max(length(missing_clusters), 3), "Set3")[seq_along(missing_clusters)]
  } else {
    colorRampPalette(brewer.pal(12, "Set3"))(length(missing_clusters))
  }
  cluster_colors[missing_clusters] <- auto_cols
}

if (file.exists(go_save_path)) {
  res_go <- qread(go_save_path)
} else {
  res_go <- list()
  for (i in all_clusters) {
    genes <- markers$gene[
      markers$cluster == i &
        markers$p_val_adj < 0.05 &
        markers$avg_log2FC > 0.5
    ]
    
    if (length(genes) > 0) {
      gene_entrez <- bitr(
        genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = "org.Hs.eg.db"
      )
      
      if (nrow(gene_entrez) > 0) {
        ego <- enrichGO(
          gene = gene_entrez$ENTREZID,
          OrgDb = org.Hs.eg.db,
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          readable = TRUE
        )
        
        if (!is.null(ego) && nrow(ego) > 0) {
          res_go[[i]] <- ego
        }
      }
    }
  }
  qsave(res_go, go_save_path)
}

if (length(res_go) > 0) {
  go_df <- as.data.frame(merge_result(res_go))
  write.xlsx(go_df, file = go_excel_path, overwrite = TRUE)
} else {
  go_df <- data.frame(
    Cluster = character(),
    Description = character(),
    p.adjust = numeric(),
    stringsAsFactors = FALSE
  )
}

if (nrow(go_df) > 0) {
  go_df$Cluster <- as.character(go_df$Cluster)
  valid_clusters <- unique(go_df$Cluster)
  
  for (i in valid_clusters) {
    plot_data <- go_df %>%
      filter(Cluster == i) %>%
      arrange(p.adjust) %>%
      head(top_n_go) %>%
      mutate(
        logP = -log10(p.adjust),
        Description = str_wrap(Description, width = 60)
      )
    
    plot_data$Description <- factor(plot_data$Description, levels = rev(plot_data$Description))
    
    p_go_bar <- ggplot(plot_data, aes(x = logP, y = Description)) +
      geom_bar(
        stat = "identity",
        width = 0.8,
        alpha = 0.8,
        fill = cluster_colors[[i]]
      ) +
      geom_text(
        aes(x = 0.05, label = Description),
        size = 5,
        hjust = 0,
        color = "black"
      ) +
      theme_minimal() +
      labs(
        title = paste0("Cluster ", i, " Top ", top_n_go, " GO Terms"),
        x = "-log10(p.adjust)",
        y = NULL
      ) +
      theme(
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.4)))
    
    safe_name <- gsub("[^a-zA-Z0-9]", "_", i)
    
    ggsave(
      file.path(go_barplot_dir, paste0("Cluster_", safe_name, "_barplot.pdf")),
      plot = p_go_bar,
      width = 8,
      height = 5
    )
  }
}

plot_list <- list()

for (i in all_clusters) {
  plot_data <- go_df %>%
    filter(Cluster == i) %>%
    arrange(p.adjust) %>%
    head(top_n_go)
  
  if (nrow(plot_data) > 0) {
    plot_data <- plot_data %>%
      mutate(y_pos = rev(seq_len(n())), x_pos = 0.05)
    
    p_block <- ggplot(plot_data) +
      annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf,
        fill = cluster_colors[[i]],
        alpha = 0.7
      ) +
      geom_text(
        aes(x = x_pos, y = y_pos, label = Description),
        size = 5,
        hjust = 0,
        color = "black",
        fontface = "bold"
      ) +
      theme_void() +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0.5, top_n_go + 0.5), expand = c(0, 0)) +
      theme(plot.margin = margin(t = 2, r = 0, b = 2, l = 0, unit = "pt"))
  } else {
    p_block <- ggplot() +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = top_n_go + 1, fill = "grey95") +
      annotate("rect", xmin = 0, xmax = 0.02, ymin = 0, ymax = top_n_go + 1, fill = cluster_colors[[i]]) +
      annotate(
        "text",
        x = 0.5,
        y = (top_n_go + 1) / 2,
        label = paste(i, ": No Enrichment"),
        size = 5,
        color = "grey50"
      ) +
      theme_void() +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0.5, top_n_go + 0.5), expand = c(0, 0))
  }
  
  plot_list[[i]] <- p_block
}

if (length(plot_list) > 0) {
  combined_plot <- wrap_plots(plot_list, ncol = 1) &
    theme(plot.margin = margin(0, 0, 0, 0))
  
  ggsave(go_combined_pdf, combined_plot, width = 6, height = 18, limitsize = FALSE)
}

message("Fig.06 done.")

# Fig.07: Sample Split and Composition ----------------------------------------
group_col <- "Study"
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

plot_df_pct <- tmp_obj@meta.data %>%
  group_by(Group = .data[[group_col]], CellType = .data[[type_col]]) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(
    Proportion = Count / sum(Count),
    Label = ifelse(Proportion > 0.03, paste0(round(Proportion * 100, 0), "%"), ""),
    Group = factor(Group, levels = levels(tmp_obj@meta.data[[group_col]])),
    Group = forcats::fct_relevel(Group, intersect("Unknown", levels(Group)))
  )

p_stack <- ggplot(plot_df_pct, aes(x = Group, y = Count, fill = CellType)) +
  geom_col(position = "fill", color = "black", width = 0.8) +
  coord_flip() +
  geom_text(aes(label = Label), position = position_fill(vjust = 0.5), size = 3.5) +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  scale_fill_manual(values = manual_colors) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = "Cell Type",
    title = paste0("Cell Type Composition per ", group_col)
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right"
  )

ggsave(
  file.path(figure_dir, paste0("Bar_Composition_", group_col, ".pdf")),
  p_stack,
  width = 8,
  height = 3
)

message("Fig.07 done.")

# Fig.08: Microglia Module Scores ---------------------------------------------
combined_module_score_path <- file.path(figure_dir, "combined_MG_ModuleScore_FeaturePlot.pdf")
plot_pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(10)

mg_marker_list <- list(
  Myelination = c(
    "PLP1", "NFASC", "BCAS1", "UGT8", "ABCA2", "MBP", "CNTN2", "ANK2", "FA2H",
    "MYRF", "DLG1", "JAM3", "MAG", "SH3TC2", "MAL", "PLLP", "ASPA", "NDRG1",
    "PTN", "PMP22", "TPPP", "NKX6-2", "QKI", "SERINC5", "MTMR2", "CD9",
    "NTRK2", "SIRT2", "DICER1"
  ),
  Axonogenesis = c(
    "NRXN1", "CNTN1", "PAK3", "DCLK1", "ROBO2", "KALRN", "MAP2", "CHL1", "SPTBN4",
    "CTNNA2", "ATP8A2", "DCC", "EPHA6", "NRCAM", "NTRK2", "PARD3", "CDH2",
    "LRRC4C", "CHN1", "MAP1B", "NRXN3", "NFIB", "ROBO1", "ANK3", "DAB1",
    "NCAM1", "PTPRO", "APP"
  ),
  Mononuclear_Cell_Proliferation = c(
    "MEF2C", "DOCK8", "ST6GAL1", "IL6ST", "CCND3", "SYK", "CD86", "TGFBR2", "CARD11"
  ),
  Signal_Transduction_By_P53_Class_Mediator = c(
    "RPL11", "RPS20", "RPS15", "RPL23", "RPL37", "RPL26", "RPS7", "RPL5", "CD74"
  ),
  Receptor_Mediated_Endocytosis = c(
    "RAB31", "STAB1", "GRB2", "SDCBP", "IGF2R", "AP2A2", "FCHO2", "ITSN1"
  ),
  Macrophage_Activation = c(
    "MFHAS1", "GPR137B", "VSIG4", "HSPD1", "CTSC", "HAVCR2"
  )
)

tmp_obj <- qread(input_obj_v2)
score_features <- character(0)

for (id in names(mg_marker_list)) {
  score_name <- paste0(id, "_Score")
  tmp_obj <- AddModuleScore(
    object = tmp_obj,
    features = list(mg_marker_list[[id]]),
    name = score_name
  )
  score_features <- c(score_features, paste0(score_name, "1"))
}

plots <- lapply(names(mg_marker_list), function(id) {
  feature_name <- paste0(id, "_Score1")
  
  FeaturePlot(
    tmp_obj,
    features = feature_name,
    order = TRUE,
    cols = plot_pal,
    raster = FALSE
  ) +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    labs(title = paste0(gsub("_", " ", id), " Score")) +
    NoAxes()
})

names(plots) <- names(mg_marker_list)

if (length(plots) > 0) {
  temp_plot <- FeaturePlot(tmp_obj, features = score_features[1], cols = plot_pal) +
    theme(legend.text = element_text(color = "transparent"))
  
  common_legend <- cowplot::get_legend(temp_plot + theme(legend.position = "right"))
  combined_module_plot <- wrap_plots(plots, ncol = 3)
  final_module_plot <- cowplot::plot_grid(
    combined_module_plot,
    common_legend,
    rel_widths = c(1, 0.08)
  )
  
  ggsave(combined_module_score_path, final_module_plot, height = 10, width = 14, bg = "#ffffff")
}

for (id in names(plots)) {
  p_single <- plots[[id]] +
    theme(legend.position = "right", legend.text = element_text(size = 8))
  
  ggsave(
    file.path(figure_dir, paste0("ModuleScore_", id, ".pdf")),
    p_single,
    height = 5,
    width = 6,
    bg = "#ffffff"
  )
}

message("Fig.08 done.")

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

# Fig.10: Literature Signature DotPlot ----------------------------------------
geneset_file <- file.path(saved_dir, "Genesets_literature.xlsx")
literature_dotplot_path <- file.path(figure_dir, "Literature_Signature_DotPlot.pdf")
plot_colors <- c("darkblue", "white", "darkred")

target_signatures_map <- c(
  "DAM.ALL" = "up_DAM_Keren_shaul_Thrupp",
  "Human_DAM_Silvin" = "Human_Adult_DAMs_Silvin",
  "Human_YAM_Silvin" = "Human_Adult_YAMs_Silvin",
  "Human_DIM_Silvin" = "Human_Adult_DIMs_Silvin"
)

stop_if_missing(geneset_file)
tmp_obj <- qread(input_obj_v2)

genesets_df <- read.xlsx(geneset_file)
genesets_df$genesets <- paste0(genesets_df$Population, "_", genesets_df$Study)

literature_marker_list <- lapply(target_signatures_map, function(sig_id) {
  genesets_df$Gene[genesets_df$genesets == sig_id]
})

score_features <- character(0)
for (name in names(literature_marker_list)) {
  tmp_obj <- AddModuleScore(
    object = tmp_obj,
    features = list(literature_marker_list[[name]]),
    name = name
  )
  score_features <- c(score_features, paste0(name, "1"))
}

p_literature <- DotPlot(
  tmp_obj,
  features = score_features,
  dot.scale = 10,
  group.by = "CellType1"
) +
  scale_colour_gradient2(low = plot_colors[1], mid = plot_colors[2], high = plot_colors[3]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = NULL, y = NULL)

ggsave(literature_dotplot_path, p_literature, height = 5, width = 4.5, bg = "#ffffff")

message("Fig.10 done.")

# Fig.11: Myelination Gene Heatmap --------------------------------------------
myelination_heatmap_path <- file.path(figure_dir, "MG_myelination_gene_heatmap.pdf")

myelination_genes <- c(
  "ST18", "SLC44A1", "PHLPP1", "PTPRD", "NCAM2", "PLP1", "NCAM1", "KCNMB4",
  "GNAO1", "UGT8", "MOBP", "MYRF", "BCAS1", "MBP", "MOG", "RNF220",
  "TRIM2", "CLDN11", "FA2H", "DLG1", "SH3TC2", "PTN"
)

tmp_obj <- qread(input_obj_v2)
avg_exp_list <- AverageExpression(tmp_obj, group.by = "CellType1", assays = "RNA")
avg_mat <- as.data.frame(avg_exp_list$RNA)

valid_genes <- intersect(myelination_genes, rownames(avg_mat))
plot_data <- avg_mat[valid_genes, , drop = FALSE]
plot_data <- plot_data[match(valid_genes, rownames(plot_data)), ]

pdf(myelination_heatmap_path, width = 3, height = 3)
pheatmap(
  plot_data,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  show_rownames = TRUE,
  legend = FALSE
)
dev.off()

message("Fig.11 done.")

# Fig.12: Other Function Gene Heatmap -----------------------------------------
other_function_heatmap_path <- file.path(figure_dir, "MG_other_function_gene_heatmap.pdf")

target_genes <- c(
  "CTNNA3", "CNTNAP4", "SPOCK1", "CADM2", "TJP1", "CNTN2", "CTNND2",
  "RNF220", "C10orf90", "FBXL7", "TRIM2", "PPM1H", "RNF144A",
  "UNC5C", "TTLL7", "MAGI2", "NTM", "NFASC", "ENPP2",
  "DNAJC6", "DLC1", "LPGAT1"
)

tmp_obj <- qread(input_obj_v2)
avg_exp_list <- AverageExpression(tmp_obj, group.by = "CellType1", assays = "RNA")
avg_mat <- as.data.frame(avg_exp_list$RNA)

valid_genes <- intersect(target_genes, rownames(avg_mat))
plot_data <- avg_mat[valid_genes, , drop = FALSE]
plot_data <- plot_data[match(valid_genes, rownames(plot_data)), ]

pdf(other_function_heatmap_path, width = 3, height = 3)
pheatmap(
  plot_data,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  show_rownames = TRUE,
  legend = FALSE
)
dev.off()

message("Fig.12 done.")

# Fig.13: Single-Cluster GSEA --------------------------------------------------
target_group_by <- "CellType1"
target_cluster <- "MG0"

custom_set_file <- file.path(gsea_data_dir, "my_gene_set.xlsx")
go_gmt_file <- file.path(gsea_data_dir, "c5.go.v2023.2.Hs.symbols.gmt")
temp_gmt_path <- file.path(gsea_data_dir, "temp_my_genesets.gmt")

gsea_qs_path <- file.path(gsea_data_dir, paste0("GSEA_kk_", target_cluster, ".qs"))
gsea_txt_path <- file.path(gsea_result_dir, paste0(target_cluster, "_combined_genesets_GSEA_result.txt"))
gsea_pdf_p_path <- file.path(gsea_result_dir, paste0(target_cluster, "_combined_GSEA_select_plot.pdf"))
gsea_pdf_clean_path <- file.path(
  gsea_result_dir,
  paste0(target_cluster, "_combined_GSEA_select_plot_no_pvalue.pdf")
)

target_go_terms <- c(
  "GOBP_POSITIVE_REGULATION_OF_MYELINATION",
  "GOBP_REGULATION_OF_MYELINATION",
  "GOBP_OLIGODENDROCYTE_DIFFERENTIATION",
  "GOBP_CELL_JUNCTION_ORGANIZATION",
  "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
  "GOBP_FATTY_ACID_TRANSPORT",
  "GOBP_MYELIN_ASSEMBLY",
  "GOBP_OLIGODENDROCYTE_DEVELOPMENT",
  "GOBP_LIPID_CATABOLIC_PROCESS"
)

selected_plot_ids <- c(
  "oligodendrocyte development",
  "oligodendrocyte differentiation",
  "myelin assembly",
  "cell junction organization"
)

plot_colors <- c("#D51F26", "#272E6A", "#208A42", "#89288F")

stop_if_missing(custom_set_file)
stop_if_missing(go_gmt_file)
tmp_obj <- qread(input_obj_v2)

gene_sets_raw <- read.xlsx(custom_set_file, colNames = FALSE, na.strings = "")
write.table(
  gene_sets_raw,
  temp_gmt_path,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

custom_gmt <- read.gmt(temp_gmt_path)
custom_gmt$gene <- toupper(custom_gmt$gene)

gmt_go <- read.gmt(go_gmt_file)
filtered_go <- gmt_go %>%
  filter(term %in% target_go_terms) %>%
  mutate(
    term = gsub("GOBP_", "", term),
    term = tolower(gsub("_", " ", term))
  )

combined_gmt <- rbind(custom_gmt, filtered_go)

Idents(tmp_obj) <- target_group_by

deg_genes <- FindMarkers(
  tmp_obj,
  ident.1 = target_cluster,
  min.pct = 0.1,
  logfc.threshold = 0,
  only.pos = FALSE
)

deg_genes$gene <- rownames(deg_genes)

ranked_genes <- deg_genes %>%
  arrange(desc(avg_log2FC)) %>%
  distinct(gene, .keep_all = TRUE)

logFC_vector <- ranked_genes$avg_log2FC
names(logFC_vector) <- toupper(ranked_genes$gene)

kk <- GSEA(
  geneList = logFC_vector,
  TERM2GENE = combined_gmt,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 1,
  maxGSSize = 500,
  scoreType = "pos",
  eps = 0,
  verbose = FALSE
)

p_gsea_pval <- gseaNb(
  object = kk,
  geneSetID = selected_plot_ids,
  subPlot = 3,
  curveCol = plot_colors,
  addPval = TRUE,
  pvalSize = 5,
  base_size = 15
)

p_gsea_clean <- gseaNb(
  object = kk,
  geneSetID = selected_plot_ids,
  subPlot = 3,
  curveCol = plot_colors,
  addPval = FALSE,
  pvalSize = 5,
  base_size = 15
)

qsave(kk, file = gsea_qs_path)

write.table(
  as.data.frame(kk),
  file = gsea_txt_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf(file = gsea_pdf_p_path, width = 12, height = 8)
print(p_gsea_pval)
dev.off()

pdf(file = gsea_pdf_clean_path, width = 12, height = 8)
print(p_gsea_clean)
dev.off()

if (file.exists(temp_gmt_path)) {
  unlink(temp_gmt_path)
}

message("Fig.13 done for cluster: ", target_cluster)

# Fig.14: Cell Cycle -----------------------------------------------------------
cell_cycle_path <- file.path(figure_dir, "cellcycle_dimplot.pdf")
phase_colors <- c("G1" = "#1F77B4", "S" = "#FF7F0E", "G2M" = "#2CA02C")

tmp_obj <- qread(input_obj_v2)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

DefaultAssay(tmp_obj) <- "RNA"
tmp_obj <- CellCycleScoring(
  tmp_obj,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)

p_cycle <- DimPlot(tmp_obj, group.by = "Phase", cols = phase_colors)
ggsave(cell_cycle_path, p_cycle, width = 6, height = 5)

message("Fig.14 done.")