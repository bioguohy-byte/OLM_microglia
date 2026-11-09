#!/usr/bin/env Rscript
# Oligodendrocyte-like microglia: figure-ordered analysis code.
# Run the relevant sections with the input objects and libraries from source/.


# Figure 1 ------------------------------------------------------------

# F1A  Study-design schematic and embedded human microglia UMAP
# Source: source/project/code_upload/Fig.01.r
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

# F1B  Nuclei and genes per nucleus by Condition
# Source: source/figure_settings/F1B_Condition.R
#!/usr/bin/env Rscript
# Figure-specific grouping of the original human microglia plotting code.
suppressPackageStartupMessages({library(Seurat);library(qs);library(dplyr);library(ggplot2);library(forcats);library(patchwork);library(ggplotify);library(tidydr);library(RColorBrewer);library(pheatmap)})
input_obj_v2 <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs")
output_root <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "F1B", "plots")
figure_dir <- file.path(output_root, "pdf")
png_dir <- file.path(output_root, "png")
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(png_dir, recursive=TRUE, showWarnings=FALSE)
group_col <- "Condition"
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


# F1C  Six-study dataset summary table

# F1D  Six functional AddModuleScore UMAPs
# Source: source/project/code_upload/Fig.01.r
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

# F1E_dot  36 selected genes (six per MG state), dot plot component
# Source: source/project/Reference/combine/scripts/20260519_dotplot.r
## 2.1 Parameters --------------------------------------------------------------
cluster_col <- 'CellType1'      
top_n       <- 5                            

# Auto-filenames
res_suffix  <- gsub("integrated_snn_", "", cluster_col)
png_path    <- paste0("results/figure/20260519_dotplot.png")
pdf_path    <- paste0("results/figure/20260519_dotplot.pdf")

## 2.2 Process Data ------------------------------------------------------------
tmp_obj <- qread('data/processed_data/mg_cca_meta_v2.qs')

genes_list <-  c("IL1RAPL1", "PCDH9", "CTNNA3", "ST18", "RNF220", "TRIM2", 
                 "NRG3", "NRXN1", "SYT1", "OPCML", "LSAMP", "NTRK3", 
                 "ANKRD44", "ITPR2", "ARHGAP24", "MEF2A", "SORL1", "P2RY12", 
                 "TMSB4X", "TPT1", "TMSB10", "FTL", "AIF1", "B2M", 
                 "SLC11A1", "CD163", "ACSL1", "NEAT1", "NAMPT", "SAT1", 
                 "SPP1", "LGMN", "EEPD1", "CTSL", "CD83", "PPARG")

## 2.3 Plotting ----------------------------------------------------------------
DefaultAssay(tmp_obj) <- 'RNA'
p_base <- DotPlot(tmp_obj, features = genes_list, group.by = cluster_col)

p_dot <- p_base$data %>%
  dplyr::rename(pct_exp = pct.exp, avg_exp = avg.exp.scaled, Gene = features.plot, Cluster = id) %>%
  ggplot(aes(x = Cluster, y = Gene, size = pct_exp, color = avg_exp)) +
  geom_point() +
  scale_color_gradient2(low = "#000000", mid = "#007965", high = "#ffcc29", midpoint = 0) +
  scale_y_discrete(limits = rev(genes_list)) +
  theme_classic() +
  labs(title = paste0("Genes Dotplot"), x = NULL, y = NULL) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.margin = margin(t = 50, r = 10, b = 10, l = 10)
  )

## 2.4 Save --------------------------------------------------------------------
ggsave(png_path, p_dot, width = 3.5, height = 9, limitsize = FALSE)
ggsave(pdf_path, p_dot, width = 3.5, height = 9, limitsize = FALSE)
message("Module 2 Done.")

# F1E_heatmap  Top50 marker heatmap component
# Source: source/project/code_upload/Fig.01.r
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

# F1E_GO  MG-state GO text blocks
# Source: source/project/code_upload/Fig.01.r
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

# F1E_SCENIC  Regulon activity heatmap component
# Source: source/project/code_upload/SCENIC.r
# SCENIC-------------------
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(qs)
library(openxlsx)
library(ggplot2)
library(pheatmap)
library(viridis)
library(patchwork)

annotation_colors <- list(
  cluster = c("MG0" = "#4DAF4A", "MG1" = "#A65628", "MG2" = "#377EB8", 
              "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG5" = "#8B0000")
)

loom <- open_loom('./results/out_SCENIC.loom') 
tmp_obj <- qread('../../data/processed_data/mg_cca_meta_v2.qs')

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)  
close_loom(loom)

sub_regulonAUC <- regulonAUC[, match(colnames(tmp_obj), colnames(regulonAUC))]

cellClusters <- data.frame(row.names = colnames(tmp_obj), 
                           sub_clusters = as.character(tmp_obj$CellType1))

qsavem(sub_regulonAUC, cellClusters, tmp_obj, file = "./results/for_rss_and_visual.qsm")

selectedResolution <- "sub_clusters"
cellsPerGroup <- split(rownames(cellClusters), cellClusters[,selectedResolution]) 

sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 

regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 

regulonActivity_byGroup_Scaled <- na.omit(regulonActivity_byGroup_Scaled)

# Fig.01 ---------------
pheatmap(regulonActivity_byGroup_Scaled)

rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellClusters[colnames(sub_regulonAUC), selectedResolution]) 
rss <- na.omit(rss) 
rssPlot <- plotRSS(rss)

# Fig.02 ---------------
p1 <- rssPlot$plot

rss_scaled <- regulonActivity_byGroup_Scaled

df <- do.call(rbind,
              lapply(1:ncol(rss_scaled), function(i){
                dat <- data.frame(
                  path  = rownames(rss_scaled),
                  cluster = colnames(rss_scaled)[i],
                  sd.1 = rss_scaled[,i],
                  sd.2 = apply(rss_scaled[,-i], 1, median)  
                )
              }))
df$fc <- df$sd.1 - df$sd.2

top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
top5 <- top5[!duplicated(top5$path), ]

rowcn <- data.frame(cluster = top5$cluster) 
n <- rss_scaled[top5$path,] 
rownames(rowcn) <- top5$path  
rownames(n) <- top5$path

write.xlsx(df, './results/SCENIC.xlsx', overwrite = T)

base_path <- "./results/figure/"
pdf_file <- paste0(base_path, "SCENIC_heatmap.pdf")
png_file <- paste0(base_path, "SCENIC_heatmap.png")

# Fig.03 ---------------
p2 <- pheatmap(n,
               color = colorRampPalette(rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", 
                                              "#ffe6b7", "#aadce0", "#72bcd5", "#528fad", 
                                              "#376795", "#1e466e")))(100),
               annotation_row = rowcn,
               annotation_colors = annotation_colors,
               show_rownames = TRUE,
               cluster_cols = FALSE,
               cluster_rows = FALSE,
               main = "MG Clusters SCENIC",
               filename = pdf_file, 
               width = 6, 
               height = 12)

pdf(file = pdf_file, width = 5, height = 8, bg = "white")
pheatmap(n,
         color = colorRampPalette(rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", 
                                        "#ffe6b7", "#aadce0", "#72bcd5", "#528fad", 
                                        "#376795", "#1e466e")))(100),
         annotation_row = rowcn,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         main = "MG Clusters SCENIC")
dev.off() 

png(file = png_file, width = 6, height = 12, units = "in", res = 300, bg = "white")
pheatmap(n,
         color = colorRampPalette(rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", 
                                        "#ffe6b7", "#aadce0", "#72bcd5", "#528fad", 
                                        "#376795", "#1e466e")))(100),
         annotation_row = rowcn,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         main = "MG Clusters SCENIC")
dev.off()


# Supplementary Figure 1 ------------------------------------------------------------

# S1A  Chia cell-type UMAP
# Source: source/project/Reference/combine/results/source_celltype_umap_20260920/revisions/20261022_Arial_square_celltype_and_marker_UMAPs/code/plot_uniform_umaps.R
#!/usr/bin/env Rscript
# Display-only revision: exact selected inputs, selected marker bars, Arial,
# one shared square view per study and an identical physical data panel.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(scattermore)
  library(scales); library(grid); library(gtable); library(digest)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task <- normalizePath(args[[1]], mustWork = TRUE)
for (d in c("output/pdf", "output/png", "audit", "runtime"))
  dir.create(file.path(task, d), recursive = TRUE, showWarnings = FALSE)
studies <- c("Chia", "Fullard", "Grubman", "Jakel", "Su", "Yao")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
annotation <- c(Chia="cell.type", Fullard="celltype", Grubman="cellType",
                Jakel="CellType", Su="CellType", Yao="AuthorMajor")
expected_n <- c(Chia=87470L, Fullard=68557L, Grubman=13214L, Jakel=17799L,
                Su=224464L, Yao=73721L)
font <- systemfonts::match_fonts("Arial")
stopifnot(grepl("Arial", basename(font$path), fixed = TRUE))
writeLines(font$path, file.path(task, "audit/Arial_font_path.txt"))
base_colours <- read.csv(file.path(task,"inputs/source_tables/base_celltype_colours.csv"), stringsAsFactors=FALSE)
base_anchors <- read.csv(file.path(task,"inputs/source_tables/base_celltype_anchors.csv"), stringsAsFactors=FALSE)
yao_colours <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_colours.csv"), stringsAsFactors=FALSE)
yao_anchors <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_anchors.csv"), stringsAsFactors=FALSE)
bars <- read.csv(file.path(task,"audit/selected_marker_bar_provenance.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(bars)==18L, all(toupper(as.character(bars$exact_match))=="TRUE"),
          all(bars$upper_limit[bars$study %in% c("Grubman","Jakel")] == 1),
          all(bars$upper_limit[bars$study == "Fullard"] == 2))
theme_uniform <- theme_classic(base_size=10, base_family="Arial") + theme(
  text=element_text(family="Arial"),
  plot.title=element_text(face="bold",size=12,margin=margin(b=3)),
  plot.subtitle=element_text(size=8,colour="#333333",margin=margin(b=7)),
  axis.title=element_text(size=10), axis.text=element_text(size=9,colour="black"),
  axis.line=element_line(linewidth=.45), axis.ticks=element_line(linewidth=.35),
  legend.title=element_text(size=8),legend.text=element_text(size=8),
  legend.position="right",legend.box.margin=margin(0,0,0,5),
  plot.margin=margin(8,8,8,8),plot.title.position="panel")
plots <- list(); cell_ids <- character(); marker_ids <- character()
coordinate_audit <- list(); scale_audit <- list(); label_audit <- list()
palette_audit <- list(); views <- list()
for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Preparing ",study)
  ct <- readRDS(file.path(task,"inputs/celltype",paste0(study,".rds")))
  mb <- readRDS(file.path(task,"inputs/marker",paste0(study,".rds")))
  stopifnot(length(ct$cells)==expected_n[[study]], !anyDuplicated(ct$cells),
    identical(rownames(ct$umap),ct$cells),identical(rownames(ct$metadata),ct$cells),
    identical(rownames(mb$umap),mb$cells),identical(colnames(mb$expression),mb$cells),
    identical(colnames(mb$counts),mb$cells),setequal(ct$cells,mb$cells))
  align <- match(mb$cells,ct$cells)
  ct_aligned <- ct$umap[align,,drop=FALSE]
  # Matrix subsetting drops optional scaled:center attributes. Compare every
  # coordinate exactly after validating dimensions and cell alignment.
  stopifnot(identical(dim(ct_aligned),dim(mb$umap)),
    identical(as.numeric(ct_aligned),as.numeric(mb$umap)))
  coords <- ct$umap
  xr <- range(coords[,1]); yr <- range(coords[,2])
  span <- max(diff(xr),diff(yr))*1.16
  xlim <- mean(xr)+c(-.5,.5)*span; ylim <- mean(yr)+c(-.5,.5)*span
  views[[study]] <- list(xlim=xlim,ylim=ylim)
  values <- as.character(ct$metadata[[annotation[[study]]]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  if (study=="Yao") {
    pc <- subset(yao_colours,plot_id=="Yao_current_UMAP_author_labels")
    anchors <- subset(yao_anchors,plot_id=="Yao_current_UMAP_author_labels")
    names(pc)[names(pc)=="author_Major"] <- "celltype"
    names(anchors)[names(anchors)=="type"] <- "celltype"
  } else {
    pc <- base_colours[base_colours$plot_id==study,,drop=FALSE]
    anchors <- base_anchors[base_anchors$plot_id==study,,drop=FALSE]
  }
  types <- pc$celltype; palette <- setNames(pc$colour,types)
  counts <- table(factor(values,levels=types))
  stopifnot(setequal(unique(values),types),identical(as.integer(counts),as.integer(pc$n_cells)),
    identical(as.integer(counts[anchors$celltype]),as.integer(anchors$n_cells)),
    all(anchors$anchor_cell %in% ct$cells))
  anchor_pos <- match(anchors$anchor_cell,ct$cells)
  stopifnot(max(abs(coords[anchor_pos,1]-anchors$x))<1e-12,
            max(abs(coords[anchor_pos,2]-anchors$y))<1e-12)
  data <- data.frame(x=coords[,1],y=coords[,2],celltype=factor(values,levels=types),row.names=ct$cells)
  subtitle <- if (study=="Yao") "73,721 cells; 16 author cell types; 174 unmatched" else
    paste0(comma(length(values))," cells; ",length(types)," cell types")
  title <- if(study=="Yao") "Yao | author cell types" else paste(study,"| cell types")
  p <- ggplot(data,aes(x=x,y=y,colour=celltype))+
    geom_scattermore(pointsize=2.2,pixels=c(1400,1400),show.legend=FALSE)+
    scale_colour_manual(values=palette,drop=FALSE)+
    geom_text_repel(data=anchors,aes(x=x,y=y,label=label_text),inherit.aes=FALSE,
      family="Arial",colour="black",size=if(length(types)>20) 2.35 else 2.75,
      lineheight=.95,seed=20261016,nudge_x=anchors$nudge_x,nudge_y=anchors$nudge_y,
      box.padding=unit(.35,"lines"),point.padding=unit(.10,"lines"),
      min.segment.length=0,segment.colour="#454545",segment.size=.20,
      force=2,force_pull=.3,max.overlaps=Inf,max.time=Inf,max.iter=20000,
      show.legend=FALSE)+
    coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
    labs(title=title,subtitle=subtitle,x="UMAP_1",y="UMAP_2")+
    theme_uniform+theme(legend.position="none")
  id <- paste0(study,"_celltype_UMAP")
  plots[[id]] <- p; cell_ids <- c(cell_ids,id)
  anchors$study <- study
  label_audit[[study]] <- anchors[,c("study","celltype","n_cells","x","y","anchor_cell","label_text","nudge_x","nudge_y")]
  palette_audit[[study]] <- data.frame(study=study,celltype=types,n_cells=as.integer(counts),colour=pc$colour)
  coordinate_audit[[study]] <- data.frame(study=study,n_cells=length(values),annotation_column=annotation[[study]],
    n_celltype_labels=length(types),celltype_and_marker_barcode_sets_identical=TRUE,
    celltype_and_marker_coordinates_identical_by_barcode=TRUE,
    coordinates_sha256=digest(coords,algo="sha256"),
    celltype_values_sha256=digest(values,algo="sha256"),
    marker_expression_sha256=digest(mb$expression,algo="sha256"),
    original_x_min=xr[1],original_x_max=xr[2],original_y_min=yr[1],original_y_max=yr[2],
    view_x_min=xlim[1],view_x_max=xlim[2],view_y_min=ylim[1],view_y_max=ylim[2],
    coordinate_ratio=1,umap_recomputed=FALSE,annotation_modified=FALSE)
  thin <- CreateSeuratObject(counts=mb$counts,min.cells=0,min.features=0,project=study)
  thin <- SetAssayData(thin,assay="RNA",layer="data",new.data=mb$expression)
  mc <- mb$umap; colnames(mc) <- c("UMAP_1","UMAP_2")
  thin[["umap"]] <- CreateDimReducObject(embeddings=mc,key="UMAP_",assay="RNA")
  stopifnot(identical(dim(Embeddings(thin,"umap")),dim(mb$umap)),
    identical(as.numeric(Embeddings(thin,"umap")),as.numeric(mb$umap)),
    identical(as.matrix(GetAssayData(thin,assay="RNA",layer="data")),as.matrix(mb$expression)))
  for (gene in markers) {
    bar <- bars[bars$study==study & bars$marker==gene,,drop=FALSE]
    stopifnot(nrow(bar)==1L)
    feature <- mb$features[[gene]]; upper <- bar$upper_limit
    tick_values <- as.numeric(strsplit(bar$tick_values,";",fixed=TRUE)[[1]])
    tick_labels <- strsplit(bar$tick_labels,";",fixed=TRUE)[[1]]
    vals <- as.numeric(mb$expression[feature,])
    p <- FeaturePlot(thin,features=feature,reduction="umap",dims=c(1,2),
      cols=c("#D3D3D3","#E31A1C"),order=TRUE,pt.size=1.3,
      min.cutoff=0,max.cutoff=upper,slot="data",raster=TRUE,raster.dpi=c(1200,1200),combine=FALSE)[[1]]+
      scale_colour_gradient(low=bar$colour_low,high=bar$colour_high,
        limits=c(0,upper),oob=scales::squish,breaks=tick_values,labels=tick_labels,
        guide=guide_colourbar(title="Log-normalized\nRNA",barheight=unit(2.6,"cm"),barwidth=unit(.30,"cm")))+
      coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
      labs(title=paste(study,"|",gene),subtitle=paste0(comma(length(vals))," cells"),x="UMAP_1",y="UMAP_2")+
      theme_uniform
    # FeaturePlot must retain the cached values after applying only the chosen cap.
    j <- match(rownames(p$data),mb$cells)
    stopifnot(!anyNA(j),nrow(p$data)==length(mb$cells),!anyDuplicated(rownames(p$data)),
      setequal(rownames(p$data),mb$cells),
      identical(as.numeric(p$data[[feature]]),pmin(pmax(vals[j],0),upper)))
    id <- paste0(study,"_",gene,"_FeaturePlot")
    plots[[id]] <- p; marker_ids <- c(marker_ids,id)
    scale_audit[[id]] <- data.frame(study=study,marker=gene,lower_limit=0,upper_limit=upper,
      visible_breaks=bar$tick_values,visible_tick_labels=bar$tick_labels,
      colour_low=bar$colour_low,colour_high=bar$colour_high,n_cells=length(vals),
      n_positive=sum(vals>0),n_saturated=sum(vals>upper),selected_pdf_sha256=bar$sha256,
      expression_identical_before_display_cap=TRUE,sorted_ascending=TRUE,
      panel_width_in=4.5,panel_height_in=4.5,coordinate_ratio=1)
  }
  rm(ct,mb,thin); gc(verbose=FALSE)
}

# Measure on the actual Cairo/Arial PDF device. Fix the panel row and column,
# then equalize every surrounding row/column so no legend can shrink a panel.
measurement <- file.path(task,"runtime/layout_measurement.pdf")
cairo_pdf(measurement,width=9,height=9,family="Arial")
grobs <- lapply(plots,ggplotGrob)
panel_mm <- 114.3 # 4.5 inches, identical in all 24 panels
shape <- vapply(grobs,function(g) paste(length(g$widths),length(g$heights)),character(1))
stopifnot(length(unique(shape))==1L)
panel_locations <- vapply(grobs,function(g) {
  loc <- g$layout[g$layout$name=="panel",]
  paste(loc$l,loc$r,loc$t,loc$b,collapse=";")
},character(1))
stopifnot(length(unique(panel_locations))==1L)
dimensions <- list()
for (id in names(grobs)) {
  g <- grobs[[id]]; loc <- g$layout[g$layout$name=="panel",]
  stopifnot(nrow(loc)==1L,loc$l==loc$r,loc$t==loc$b)
  g$widths[loc$l] <- unit(panel_mm,"mm");g$heights[loc$t] <- unit(panel_mm,"mm")
  grobs[[id]] <- g
}
all_widths <- lapply(grobs,function(g) convertWidth(g$widths,"in",valueOnly=TRUE))
all_heights <- lapply(grobs,function(g) convertHeight(g$heights,"in",valueOnly=TRUE))
common_widths <- do.call(pmax,all_widths)
common_heights <- do.call(pmax,all_heights)
for (id in names(grobs)) {
  g <- grobs[[id]];g$widths <- unit(common_widths,"in");g$heights <- unit(common_heights,"in")
  grobs[[id]] <- g
  loc <- g$layout[g$layout$name=="panel",]
  pw <- convertWidth(g$widths[loc$l],"mm",valueOnly=TRUE)
  ph <- convertHeight(g$heights[loc$t],"mm",valueOnly=TRUE)
  stopifnot(abs(pw-panel_mm)<1e-9,abs(ph-panel_mm)<1e-9)
  dimensions[[id]] <- data.frame(plot_id=id,panel_width_mm=pw,panel_height_mm=ph,
    outer_width_in=sum(common_widths),outer_height_in=sum(common_heights),font="Arial")
}
dev.off();unlink(measurement)
slot_width <- sum(common_widths);slot_height <- sum(common_heights)
draw_page <- function(ids,ncol,nrow) {
  stopifnot(length(ids)==ncol*nrow)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol,
    widths=unit(rep(slot_width,ncol),"in"),heights=unit(rep(slot_height,nrow),"in"))))
  for (i in seq_along(ids)) {
    pushViewport(viewport(layout.pos.row=(i-1L)%/%ncol+1L,layout.pos.col=(i-1L)%%ncol+1L))
    grid.draw(grobs[[ids[[i]]]])
    upViewport()
  }
  upViewport()
}
outputs <- list()
save_page <- function(ids,name,ncol=1L,nrow=1L) {
  w <- slot_width*ncol;h <- slot_height*nrow
  message(format(Sys.time(),"%F %T")," Saving ",name," (",round(w,2)," x ",round(h,2)," in)")
  cairo_pdf(file.path(task,"output/pdf",paste0(name,".pdf")),width=w,height=h,family="Arial")
  draw_page(ids,ncol,nrow);dev.off()
  ragg::agg_png(file.path(task,"output/png",paste0(name,".png")),width=w,height=h,units="in",res=240,background="white")
  draw_page(ids,ncol,nrow);dev.off()
  outputs[[name]] <<- data.frame(name=name,n_panels=length(ids),n_columns=ncol,n_rows=nrow,
    width_in=w,height_in=h,pdf=paste0("output/pdf/",name,".pdf"),png=paste0("output/png/",name,".png"))
}
for (id in c(cell_ids,marker_ids)) save_page(id,id)
save_page(cell_ids,"Six_study_celltype_UMAPs_Arial_square",3L,2L)
save_page(marker_ids[1:9],"Marker_UMAPs_3x3_page1_Chia_Fullard_Grubman",3L,3L)
save_page(marker_ids[10:18],"Marker_UMAPs_3x3_page2_Jakel_Su_Yao",3L,3L)
write.csv(do.call(rbind,coordinate_audit),file.path(task,"audit/coordinate_and_annotation_preservation.csv"),row.names=FALSE)
write.csv(do.call(rbind,scale_audit),file.path(task,"audit/marker_scale_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,label_audit),file.path(task,"audit/celltype_labels_and_anchors.csv"),row.names=FALSE)
write.csv(do.call(rbind,palette_audit),file.path(task,"audit/celltype_counts_and_colours.csv"),row.names=FALSE)
write.csv(do.call(rbind,dimensions),file.path(task,"audit/physical_panel_dimensions.csv"),row.names=FALSE)
write.csv(do.call(rbind,outputs),file.path(task,"audit/output_manifest.csv"),row.names=FALSE)
saveRDS(views,file.path(task,"audit/shared_square_view_limits.rds"))
writeLines(capture.output(sessionInfo()),file.path(task,"audit/plotting_sessionInfo.txt"))
message("UNIFORM_UMAP_PLOTTING_COMPLETE: 27 PDF and PNG files; 24 equal 114.3-mm square panels")

# S1B  Fullard cell-type UMAP
# Source: source/project/Reference/combine/results/source_celltype_umap_20260920/revisions/20261022_Arial_square_celltype_and_marker_UMAPs/code/plot_uniform_umaps.R
#!/usr/bin/env Rscript
# Display-only revision: exact selected inputs, selected marker bars, Arial,
# one shared square view per study and an identical physical data panel.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(scattermore)
  library(scales); library(grid); library(gtable); library(digest)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task <- normalizePath(args[[1]], mustWork = TRUE)
for (d in c("output/pdf", "output/png", "audit", "runtime"))
  dir.create(file.path(task, d), recursive = TRUE, showWarnings = FALSE)
studies <- c("Chia", "Fullard", "Grubman", "Jakel", "Su", "Yao")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
annotation <- c(Chia="cell.type", Fullard="celltype", Grubman="cellType",
                Jakel="CellType", Su="CellType", Yao="AuthorMajor")
expected_n <- c(Chia=87470L, Fullard=68557L, Grubman=13214L, Jakel=17799L,
                Su=224464L, Yao=73721L)
font <- systemfonts::match_fonts("Arial")
stopifnot(grepl("Arial", basename(font$path), fixed = TRUE))
writeLines(font$path, file.path(task, "audit/Arial_font_path.txt"))
base_colours <- read.csv(file.path(task,"inputs/source_tables/base_celltype_colours.csv"), stringsAsFactors=FALSE)
base_anchors <- read.csv(file.path(task,"inputs/source_tables/base_celltype_anchors.csv"), stringsAsFactors=FALSE)
yao_colours <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_colours.csv"), stringsAsFactors=FALSE)
yao_anchors <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_anchors.csv"), stringsAsFactors=FALSE)
bars <- read.csv(file.path(task,"audit/selected_marker_bar_provenance.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(bars)==18L, all(toupper(as.character(bars$exact_match))=="TRUE"),
          all(bars$upper_limit[bars$study %in% c("Grubman","Jakel")] == 1),
          all(bars$upper_limit[bars$study == "Fullard"] == 2))
theme_uniform <- theme_classic(base_size=10, base_family="Arial") + theme(
  text=element_text(family="Arial"),
  plot.title=element_text(face="bold",size=12,margin=margin(b=3)),
  plot.subtitle=element_text(size=8,colour="#333333",margin=margin(b=7)),
  axis.title=element_text(size=10), axis.text=element_text(size=9,colour="black"),
  axis.line=element_line(linewidth=.45), axis.ticks=element_line(linewidth=.35),
  legend.title=element_text(size=8),legend.text=element_text(size=8),
  legend.position="right",legend.box.margin=margin(0,0,0,5),
  plot.margin=margin(8,8,8,8),plot.title.position="panel")
plots <- list(); cell_ids <- character(); marker_ids <- character()
coordinate_audit <- list(); scale_audit <- list(); label_audit <- list()
palette_audit <- list(); views <- list()
for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Preparing ",study)
  ct <- readRDS(file.path(task,"inputs/celltype",paste0(study,".rds")))
  mb <- readRDS(file.path(task,"inputs/marker",paste0(study,".rds")))
  stopifnot(length(ct$cells)==expected_n[[study]], !anyDuplicated(ct$cells),
    identical(rownames(ct$umap),ct$cells),identical(rownames(ct$metadata),ct$cells),
    identical(rownames(mb$umap),mb$cells),identical(colnames(mb$expression),mb$cells),
    identical(colnames(mb$counts),mb$cells),setequal(ct$cells,mb$cells))
  align <- match(mb$cells,ct$cells)
  ct_aligned <- ct$umap[align,,drop=FALSE]
  # Matrix subsetting drops optional scaled:center attributes. Compare every
  # coordinate exactly after validating dimensions and cell alignment.
  stopifnot(identical(dim(ct_aligned),dim(mb$umap)),
    identical(as.numeric(ct_aligned),as.numeric(mb$umap)))
  coords <- ct$umap
  xr <- range(coords[,1]); yr <- range(coords[,2])
  span <- max(diff(xr),diff(yr))*1.16
  xlim <- mean(xr)+c(-.5,.5)*span; ylim <- mean(yr)+c(-.5,.5)*span
  views[[study]] <- list(xlim=xlim,ylim=ylim)
  values <- as.character(ct$metadata[[annotation[[study]]]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  if (study=="Yao") {
    pc <- subset(yao_colours,plot_id=="Yao_current_UMAP_author_labels")
    anchors <- subset(yao_anchors,plot_id=="Yao_current_UMAP_author_labels")
    names(pc)[names(pc)=="author_Major"] <- "celltype"
    names(anchors)[names(anchors)=="type"] <- "celltype"
  } else {
    pc <- base_colours[base_colours$plot_id==study,,drop=FALSE]
    anchors <- base_anchors[base_anchors$plot_id==study,,drop=FALSE]
  }
  types <- pc$celltype; palette <- setNames(pc$colour,types)
  counts <- table(factor(values,levels=types))
  stopifnot(setequal(unique(values),types),identical(as.integer(counts),as.integer(pc$n_cells)),
    identical(as.integer(counts[anchors$celltype]),as.integer(anchors$n_cells)),
    all(anchors$anchor_cell %in% ct$cells))
  anchor_pos <- match(anchors$anchor_cell,ct$cells)
  stopifnot(max(abs(coords[anchor_pos,1]-anchors$x))<1e-12,
            max(abs(coords[anchor_pos,2]-anchors$y))<1e-12)
  data <- data.frame(x=coords[,1],y=coords[,2],celltype=factor(values,levels=types),row.names=ct$cells)
  subtitle <- if (study=="Yao") "73,721 cells; 16 author cell types; 174 unmatched" else
    paste0(comma(length(values))," cells; ",length(types)," cell types")
  title <- if(study=="Yao") "Yao | author cell types" else paste(study,"| cell types")
  p <- ggplot(data,aes(x=x,y=y,colour=celltype))+
    geom_scattermore(pointsize=2.2,pixels=c(1400,1400),show.legend=FALSE)+
    scale_colour_manual(values=palette,drop=FALSE)+
    geom_text_repel(data=anchors,aes(x=x,y=y,label=label_text),inherit.aes=FALSE,
      family="Arial",colour="black",size=if(length(types)>20) 2.35 else 2.75,
      lineheight=.95,seed=20261016,nudge_x=anchors$nudge_x,nudge_y=anchors$nudge_y,
      box.padding=unit(.35,"lines"),point.padding=unit(.10,"lines"),
      min.segment.length=0,segment.colour="#454545",segment.size=.20,
      force=2,force_pull=.3,max.overlaps=Inf,max.time=Inf,max.iter=20000,
      show.legend=FALSE)+
    coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
    labs(title=title,subtitle=subtitle,x="UMAP_1",y="UMAP_2")+
    theme_uniform+theme(legend.position="none")
  id <- paste0(study,"_celltype_UMAP")
  plots[[id]] <- p; cell_ids <- c(cell_ids,id)
  anchors$study <- study
  label_audit[[study]] <- anchors[,c("study","celltype","n_cells","x","y","anchor_cell","label_text","nudge_x","nudge_y")]
  palette_audit[[study]] <- data.frame(study=study,celltype=types,n_cells=as.integer(counts),colour=pc$colour)
  coordinate_audit[[study]] <- data.frame(study=study,n_cells=length(values),annotation_column=annotation[[study]],
    n_celltype_labels=length(types),celltype_and_marker_barcode_sets_identical=TRUE,
    celltype_and_marker_coordinates_identical_by_barcode=TRUE,
    coordinates_sha256=digest(coords,algo="sha256"),
    celltype_values_sha256=digest(values,algo="sha256"),
    marker_expression_sha256=digest(mb$expression,algo="sha256"),
    original_x_min=xr[1],original_x_max=xr[2],original_y_min=yr[1],original_y_max=yr[2],
    view_x_min=xlim[1],view_x_max=xlim[2],view_y_min=ylim[1],view_y_max=ylim[2],
    coordinate_ratio=1,umap_recomputed=FALSE,annotation_modified=FALSE)
  thin <- CreateSeuratObject(counts=mb$counts,min.cells=0,min.features=0,project=study)
  thin <- SetAssayData(thin,assay="RNA",layer="data",new.data=mb$expression)
  mc <- mb$umap; colnames(mc) <- c("UMAP_1","UMAP_2")
  thin[["umap"]] <- CreateDimReducObject(embeddings=mc,key="UMAP_",assay="RNA")
  stopifnot(identical(dim(Embeddings(thin,"umap")),dim(mb$umap)),
    identical(as.numeric(Embeddings(thin,"umap")),as.numeric(mb$umap)),
    identical(as.matrix(GetAssayData(thin,assay="RNA",layer="data")),as.matrix(mb$expression)))
  for (gene in markers) {
    bar <- bars[bars$study==study & bars$marker==gene,,drop=FALSE]
    stopifnot(nrow(bar)==1L)
    feature <- mb$features[[gene]]; upper <- bar$upper_limit
    tick_values <- as.numeric(strsplit(bar$tick_values,";",fixed=TRUE)[[1]])
    tick_labels <- strsplit(bar$tick_labels,";",fixed=TRUE)[[1]]
    vals <- as.numeric(mb$expression[feature,])
    p <- FeaturePlot(thin,features=feature,reduction="umap",dims=c(1,2),
      cols=c("#D3D3D3","#E31A1C"),order=TRUE,pt.size=1.3,
      min.cutoff=0,max.cutoff=upper,slot="data",raster=TRUE,raster.dpi=c(1200,1200),combine=FALSE)[[1]]+
      scale_colour_gradient(low=bar$colour_low,high=bar$colour_high,
        limits=c(0,upper),oob=scales::squish,breaks=tick_values,labels=tick_labels,
        guide=guide_colourbar(title="Log-normalized\nRNA",barheight=unit(2.6,"cm"),barwidth=unit(.30,"cm")))+
      coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
      labs(title=paste(study,"|",gene),subtitle=paste0(comma(length(vals))," cells"),x="UMAP_1",y="UMAP_2")+
      theme_uniform
    # FeaturePlot must retain the cached values after applying only the chosen cap.
    j <- match(rownames(p$data),mb$cells)
    stopifnot(!anyNA(j),nrow(p$data)==length(mb$cells),!anyDuplicated(rownames(p$data)),
      setequal(rownames(p$data),mb$cells),
      identical(as.numeric(p$data[[feature]]),pmin(pmax(vals[j],0),upper)))
    id <- paste0(study,"_",gene,"_FeaturePlot")
    plots[[id]] <- p; marker_ids <- c(marker_ids,id)
    scale_audit[[id]] <- data.frame(study=study,marker=gene,lower_limit=0,upper_limit=upper,
      visible_breaks=bar$tick_values,visible_tick_labels=bar$tick_labels,
      colour_low=bar$colour_low,colour_high=bar$colour_high,n_cells=length(vals),
      n_positive=sum(vals>0),n_saturated=sum(vals>upper),selected_pdf_sha256=bar$sha256,
      expression_identical_before_display_cap=TRUE,sorted_ascending=TRUE,
      panel_width_in=4.5,panel_height_in=4.5,coordinate_ratio=1)
  }
  rm(ct,mb,thin); gc(verbose=FALSE)
}

# Measure on the actual Cairo/Arial PDF device. Fix the panel row and column,
# then equalize every surrounding row/column so no legend can shrink a panel.
measurement <- file.path(task,"runtime/layout_measurement.pdf")
cairo_pdf(measurement,width=9,height=9,family="Arial")
grobs <- lapply(plots,ggplotGrob)
panel_mm <- 114.3 # 4.5 inches, identical in all 24 panels
shape <- vapply(grobs,function(g) paste(length(g$widths),length(g$heights)),character(1))
stopifnot(length(unique(shape))==1L)
panel_locations <- vapply(grobs,function(g) {
  loc <- g$layout[g$layout$name=="panel",]
  paste(loc$l,loc$r,loc$t,loc$b,collapse=";")
},character(1))
stopifnot(length(unique(panel_locations))==1L)
dimensions <- list()
for (id in names(grobs)) {
  g <- grobs[[id]]; loc <- g$layout[g$layout$name=="panel",]
  stopifnot(nrow(loc)==1L,loc$l==loc$r,loc$t==loc$b)
  g$widths[loc$l] <- unit(panel_mm,"mm");g$heights[loc$t] <- unit(panel_mm,"mm")
  grobs[[id]] <- g
}
all_widths <- lapply(grobs,function(g) convertWidth(g$widths,"in",valueOnly=TRUE))
all_heights <- lapply(grobs,function(g) convertHeight(g$heights,"in",valueOnly=TRUE))
common_widths <- do.call(pmax,all_widths)
common_heights <- do.call(pmax,all_heights)
for (id in names(grobs)) {
  g <- grobs[[id]];g$widths <- unit(common_widths,"in");g$heights <- unit(common_heights,"in")
  grobs[[id]] <- g
  loc <- g$layout[g$layout$name=="panel",]
  pw <- convertWidth(g$widths[loc$l],"mm",valueOnly=TRUE)
  ph <- convertHeight(g$heights[loc$t],"mm",valueOnly=TRUE)
  stopifnot(abs(pw-panel_mm)<1e-9,abs(ph-panel_mm)<1e-9)
  dimensions[[id]] <- data.frame(plot_id=id,panel_width_mm=pw,panel_height_mm=ph,
    outer_width_in=sum(common_widths),outer_height_in=sum(common_heights),font="Arial")
}
dev.off();unlink(measurement)
slot_width <- sum(common_widths);slot_height <- sum(common_heights)
draw_page <- function(ids,ncol,nrow) {
  stopifnot(length(ids)==ncol*nrow)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol,
    widths=unit(rep(slot_width,ncol),"in"),heights=unit(rep(slot_height,nrow),"in"))))
  for (i in seq_along(ids)) {
    pushViewport(viewport(layout.pos.row=(i-1L)%/%ncol+1L,layout.pos.col=(i-1L)%%ncol+1L))
    grid.draw(grobs[[ids[[i]]]])
    upViewport()
  }
  upViewport()
}
outputs <- list()
save_page <- function(ids,name,ncol=1L,nrow=1L) {
  w <- slot_width*ncol;h <- slot_height*nrow
  message(format(Sys.time(),"%F %T")," Saving ",name," (",round(w,2)," x ",round(h,2)," in)")
  cairo_pdf(file.path(task,"output/pdf",paste0(name,".pdf")),width=w,height=h,family="Arial")
  draw_page(ids,ncol,nrow);dev.off()
  ragg::agg_png(file.path(task,"output/png",paste0(name,".png")),width=w,height=h,units="in",res=240,background="white")
  draw_page(ids,ncol,nrow);dev.off()
  outputs[[name]] <<- data.frame(name=name,n_panels=length(ids),n_columns=ncol,n_rows=nrow,
    width_in=w,height_in=h,pdf=paste0("output/pdf/",name,".pdf"),png=paste0("output/png/",name,".png"))
}
for (id in c(cell_ids,marker_ids)) save_page(id,id)
save_page(cell_ids,"Six_study_celltype_UMAPs_Arial_square",3L,2L)
save_page(marker_ids[1:9],"Marker_UMAPs_3x3_page1_Chia_Fullard_Grubman",3L,3L)
save_page(marker_ids[10:18],"Marker_UMAPs_3x3_page2_Jakel_Su_Yao",3L,3L)
write.csv(do.call(rbind,coordinate_audit),file.path(task,"audit/coordinate_and_annotation_preservation.csv"),row.names=FALSE)
write.csv(do.call(rbind,scale_audit),file.path(task,"audit/marker_scale_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,label_audit),file.path(task,"audit/celltype_labels_and_anchors.csv"),row.names=FALSE)
write.csv(do.call(rbind,palette_audit),file.path(task,"audit/celltype_counts_and_colours.csv"),row.names=FALSE)
write.csv(do.call(rbind,dimensions),file.path(task,"audit/physical_panel_dimensions.csv"),row.names=FALSE)
write.csv(do.call(rbind,outputs),file.path(task,"audit/output_manifest.csv"),row.names=FALSE)
saveRDS(views,file.path(task,"audit/shared_square_view_limits.rds"))
writeLines(capture.output(sessionInfo()),file.path(task,"audit/plotting_sessionInfo.txt"))
message("UNIFORM_UMAP_PLOTTING_COMPLETE: 27 PDF and PNG files; 24 equal 114.3-mm square panels")

# S1C  Grubman cell-type UMAP
# Source: source/project/Reference/combine/results/source_celltype_umap_20260920/revisions/20261022_Arial_square_celltype_and_marker_UMAPs/code/plot_uniform_umaps.R
#!/usr/bin/env Rscript
# Display-only revision: exact selected inputs, selected marker bars, Arial,
# one shared square view per study and an identical physical data panel.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(scattermore)
  library(scales); library(grid); library(gtable); library(digest)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task <- normalizePath(args[[1]], mustWork = TRUE)
for (d in c("output/pdf", "output/png", "audit", "runtime"))
  dir.create(file.path(task, d), recursive = TRUE, showWarnings = FALSE)
studies <- c("Chia", "Fullard", "Grubman", "Jakel", "Su", "Yao")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
annotation <- c(Chia="cell.type", Fullard="celltype", Grubman="cellType",
                Jakel="CellType", Su="CellType", Yao="AuthorMajor")
expected_n <- c(Chia=87470L, Fullard=68557L, Grubman=13214L, Jakel=17799L,
                Su=224464L, Yao=73721L)
font <- systemfonts::match_fonts("Arial")
stopifnot(grepl("Arial", basename(font$path), fixed = TRUE))
writeLines(font$path, file.path(task, "audit/Arial_font_path.txt"))
base_colours <- read.csv(file.path(task,"inputs/source_tables/base_celltype_colours.csv"), stringsAsFactors=FALSE)
base_anchors <- read.csv(file.path(task,"inputs/source_tables/base_celltype_anchors.csv"), stringsAsFactors=FALSE)
yao_colours <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_colours.csv"), stringsAsFactors=FALSE)
yao_anchors <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_anchors.csv"), stringsAsFactors=FALSE)
bars <- read.csv(file.path(task,"audit/selected_marker_bar_provenance.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(bars)==18L, all(toupper(as.character(bars$exact_match))=="TRUE"),
          all(bars$upper_limit[bars$study %in% c("Grubman","Jakel")] == 1),
          all(bars$upper_limit[bars$study == "Fullard"] == 2))
theme_uniform <- theme_classic(base_size=10, base_family="Arial") + theme(
  text=element_text(family="Arial"),
  plot.title=element_text(face="bold",size=12,margin=margin(b=3)),
  plot.subtitle=element_text(size=8,colour="#333333",margin=margin(b=7)),
  axis.title=element_text(size=10), axis.text=element_text(size=9,colour="black"),
  axis.line=element_line(linewidth=.45), axis.ticks=element_line(linewidth=.35),
  legend.title=element_text(size=8),legend.text=element_text(size=8),
  legend.position="right",legend.box.margin=margin(0,0,0,5),
  plot.margin=margin(8,8,8,8),plot.title.position="panel")
plots <- list(); cell_ids <- character(); marker_ids <- character()
coordinate_audit <- list(); scale_audit <- list(); label_audit <- list()
palette_audit <- list(); views <- list()
for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Preparing ",study)
  ct <- readRDS(file.path(task,"inputs/celltype",paste0(study,".rds")))
  mb <- readRDS(file.path(task,"inputs/marker",paste0(study,".rds")))
  stopifnot(length(ct$cells)==expected_n[[study]], !anyDuplicated(ct$cells),
    identical(rownames(ct$umap),ct$cells),identical(rownames(ct$metadata),ct$cells),
    identical(rownames(mb$umap),mb$cells),identical(colnames(mb$expression),mb$cells),
    identical(colnames(mb$counts),mb$cells),setequal(ct$cells,mb$cells))
  align <- match(mb$cells,ct$cells)
  ct_aligned <- ct$umap[align,,drop=FALSE]
  # Matrix subsetting drops optional scaled:center attributes. Compare every
  # coordinate exactly after validating dimensions and cell alignment.
  stopifnot(identical(dim(ct_aligned),dim(mb$umap)),
    identical(as.numeric(ct_aligned),as.numeric(mb$umap)))
  coords <- ct$umap
  xr <- range(coords[,1]); yr <- range(coords[,2])
  span <- max(diff(xr),diff(yr))*1.16
  xlim <- mean(xr)+c(-.5,.5)*span; ylim <- mean(yr)+c(-.5,.5)*span
  views[[study]] <- list(xlim=xlim,ylim=ylim)
  values <- as.character(ct$metadata[[annotation[[study]]]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  if (study=="Yao") {
    pc <- subset(yao_colours,plot_id=="Yao_current_UMAP_author_labels")
    anchors <- subset(yao_anchors,plot_id=="Yao_current_UMAP_author_labels")
    names(pc)[names(pc)=="author_Major"] <- "celltype"
    names(anchors)[names(anchors)=="type"] <- "celltype"
  } else {
    pc <- base_colours[base_colours$plot_id==study,,drop=FALSE]
    anchors <- base_anchors[base_anchors$plot_id==study,,drop=FALSE]
  }
  types <- pc$celltype; palette <- setNames(pc$colour,types)
  counts <- table(factor(values,levels=types))
  stopifnot(setequal(unique(values),types),identical(as.integer(counts),as.integer(pc$n_cells)),
    identical(as.integer(counts[anchors$celltype]),as.integer(anchors$n_cells)),
    all(anchors$anchor_cell %in% ct$cells))
  anchor_pos <- match(anchors$anchor_cell,ct$cells)
  stopifnot(max(abs(coords[anchor_pos,1]-anchors$x))<1e-12,
            max(abs(coords[anchor_pos,2]-anchors$y))<1e-12)
  data <- data.frame(x=coords[,1],y=coords[,2],celltype=factor(values,levels=types),row.names=ct$cells)
  subtitle <- if (study=="Yao") "73,721 cells; 16 author cell types; 174 unmatched" else
    paste0(comma(length(values))," cells; ",length(types)," cell types")
  title <- if(study=="Yao") "Yao | author cell types" else paste(study,"| cell types")
  p <- ggplot(data,aes(x=x,y=y,colour=celltype))+
    geom_scattermore(pointsize=2.2,pixels=c(1400,1400),show.legend=FALSE)+
    scale_colour_manual(values=palette,drop=FALSE)+
    geom_text_repel(data=anchors,aes(x=x,y=y,label=label_text),inherit.aes=FALSE,
      family="Arial",colour="black",size=if(length(types)>20) 2.35 else 2.75,
      lineheight=.95,seed=20261016,nudge_x=anchors$nudge_x,nudge_y=anchors$nudge_y,
      box.padding=unit(.35,"lines"),point.padding=unit(.10,"lines"),
      min.segment.length=0,segment.colour="#454545",segment.size=.20,
      force=2,force_pull=.3,max.overlaps=Inf,max.time=Inf,max.iter=20000,
      show.legend=FALSE)+
    coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
    labs(title=title,subtitle=subtitle,x="UMAP_1",y="UMAP_2")+
    theme_uniform+theme(legend.position="none")
  id <- paste0(study,"_celltype_UMAP")
  plots[[id]] <- p; cell_ids <- c(cell_ids,id)
  anchors$study <- study
  label_audit[[study]] <- anchors[,c("study","celltype","n_cells","x","y","anchor_cell","label_text","nudge_x","nudge_y")]
  palette_audit[[study]] <- data.frame(study=study,celltype=types,n_cells=as.integer(counts),colour=pc$colour)
  coordinate_audit[[study]] <- data.frame(study=study,n_cells=length(values),annotation_column=annotation[[study]],
    n_celltype_labels=length(types),celltype_and_marker_barcode_sets_identical=TRUE,
    celltype_and_marker_coordinates_identical_by_barcode=TRUE,
    coordinates_sha256=digest(coords,algo="sha256"),
    celltype_values_sha256=digest(values,algo="sha256"),
    marker_expression_sha256=digest(mb$expression,algo="sha256"),
    original_x_min=xr[1],original_x_max=xr[2],original_y_min=yr[1],original_y_max=yr[2],
    view_x_min=xlim[1],view_x_max=xlim[2],view_y_min=ylim[1],view_y_max=ylim[2],
    coordinate_ratio=1,umap_recomputed=FALSE,annotation_modified=FALSE)
  thin <- CreateSeuratObject(counts=mb$counts,min.cells=0,min.features=0,project=study)
  thin <- SetAssayData(thin,assay="RNA",layer="data",new.data=mb$expression)
  mc <- mb$umap; colnames(mc) <- c("UMAP_1","UMAP_2")
  thin[["umap"]] <- CreateDimReducObject(embeddings=mc,key="UMAP_",assay="RNA")
  stopifnot(identical(dim(Embeddings(thin,"umap")),dim(mb$umap)),
    identical(as.numeric(Embeddings(thin,"umap")),as.numeric(mb$umap)),
    identical(as.matrix(GetAssayData(thin,assay="RNA",layer="data")),as.matrix(mb$expression)))
  for (gene in markers) {
    bar <- bars[bars$study==study & bars$marker==gene,,drop=FALSE]
    stopifnot(nrow(bar)==1L)
    feature <- mb$features[[gene]]; upper <- bar$upper_limit
    tick_values <- as.numeric(strsplit(bar$tick_values,";",fixed=TRUE)[[1]])
    tick_labels <- strsplit(bar$tick_labels,";",fixed=TRUE)[[1]]
    vals <- as.numeric(mb$expression[feature,])
    p <- FeaturePlot(thin,features=feature,reduction="umap",dims=c(1,2),
      cols=c("#D3D3D3","#E31A1C"),order=TRUE,pt.size=1.3,
      min.cutoff=0,max.cutoff=upper,slot="data",raster=TRUE,raster.dpi=c(1200,1200),combine=FALSE)[[1]]+
      scale_colour_gradient(low=bar$colour_low,high=bar$colour_high,
        limits=c(0,upper),oob=scales::squish,breaks=tick_values,labels=tick_labels,
        guide=guide_colourbar(title="Log-normalized\nRNA",barheight=unit(2.6,"cm"),barwidth=unit(.30,"cm")))+
      coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
      labs(title=paste(study,"|",gene),subtitle=paste0(comma(length(vals))," cells"),x="UMAP_1",y="UMAP_2")+
      theme_uniform
    # FeaturePlot must retain the cached values after applying only the chosen cap.
    j <- match(rownames(p$data),mb$cells)
    stopifnot(!anyNA(j),nrow(p$data)==length(mb$cells),!anyDuplicated(rownames(p$data)),
      setequal(rownames(p$data),mb$cells),
      identical(as.numeric(p$data[[feature]]),pmin(pmax(vals[j],0),upper)))
    id <- paste0(study,"_",gene,"_FeaturePlot")
    plots[[id]] <- p; marker_ids <- c(marker_ids,id)
    scale_audit[[id]] <- data.frame(study=study,marker=gene,lower_limit=0,upper_limit=upper,
      visible_breaks=bar$tick_values,visible_tick_labels=bar$tick_labels,
      colour_low=bar$colour_low,colour_high=bar$colour_high,n_cells=length(vals),
      n_positive=sum(vals>0),n_saturated=sum(vals>upper),selected_pdf_sha256=bar$sha256,
      expression_identical_before_display_cap=TRUE,sorted_ascending=TRUE,
      panel_width_in=4.5,panel_height_in=4.5,coordinate_ratio=1)
  }
  rm(ct,mb,thin); gc(verbose=FALSE)
}

# Measure on the actual Cairo/Arial PDF device. Fix the panel row and column,
# then equalize every surrounding row/column so no legend can shrink a panel.
measurement <- file.path(task,"runtime/layout_measurement.pdf")
cairo_pdf(measurement,width=9,height=9,family="Arial")
grobs <- lapply(plots,ggplotGrob)
panel_mm <- 114.3 # 4.5 inches, identical in all 24 panels
shape <- vapply(grobs,function(g) paste(length(g$widths),length(g$heights)),character(1))
stopifnot(length(unique(shape))==1L)
panel_locations <- vapply(grobs,function(g) {
  loc <- g$layout[g$layout$name=="panel",]
  paste(loc$l,loc$r,loc$t,loc$b,collapse=";")
},character(1))
stopifnot(length(unique(panel_locations))==1L)
dimensions <- list()
for (id in names(grobs)) {
  g <- grobs[[id]]; loc <- g$layout[g$layout$name=="panel",]
  stopifnot(nrow(loc)==1L,loc$l==loc$r,loc$t==loc$b)
  g$widths[loc$l] <- unit(panel_mm,"mm");g$heights[loc$t] <- unit(panel_mm,"mm")
  grobs[[id]] <- g
}
all_widths <- lapply(grobs,function(g) convertWidth(g$widths,"in",valueOnly=TRUE))
all_heights <- lapply(grobs,function(g) convertHeight(g$heights,"in",valueOnly=TRUE))
common_widths <- do.call(pmax,all_widths)
common_heights <- do.call(pmax,all_heights)
for (id in names(grobs)) {
  g <- grobs[[id]];g$widths <- unit(common_widths,"in");g$heights <- unit(common_heights,"in")
  grobs[[id]] <- g
  loc <- g$layout[g$layout$name=="panel",]
  pw <- convertWidth(g$widths[loc$l],"mm",valueOnly=TRUE)
  ph <- convertHeight(g$heights[loc$t],"mm",valueOnly=TRUE)
  stopifnot(abs(pw-panel_mm)<1e-9,abs(ph-panel_mm)<1e-9)
  dimensions[[id]] <- data.frame(plot_id=id,panel_width_mm=pw,panel_height_mm=ph,
    outer_width_in=sum(common_widths),outer_height_in=sum(common_heights),font="Arial")
}
dev.off();unlink(measurement)
slot_width <- sum(common_widths);slot_height <- sum(common_heights)
draw_page <- function(ids,ncol,nrow) {
  stopifnot(length(ids)==ncol*nrow)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol,
    widths=unit(rep(slot_width,ncol),"in"),heights=unit(rep(slot_height,nrow),"in"))))
  for (i in seq_along(ids)) {
    pushViewport(viewport(layout.pos.row=(i-1L)%/%ncol+1L,layout.pos.col=(i-1L)%%ncol+1L))
    grid.draw(grobs[[ids[[i]]]])
    upViewport()
  }
  upViewport()
}
outputs <- list()
save_page <- function(ids,name,ncol=1L,nrow=1L) {
  w <- slot_width*ncol;h <- slot_height*nrow
  message(format(Sys.time(),"%F %T")," Saving ",name," (",round(w,2)," x ",round(h,2)," in)")
  cairo_pdf(file.path(task,"output/pdf",paste0(name,".pdf")),width=w,height=h,family="Arial")
  draw_page(ids,ncol,nrow);dev.off()
  ragg::agg_png(file.path(task,"output/png",paste0(name,".png")),width=w,height=h,units="in",res=240,background="white")
  draw_page(ids,ncol,nrow);dev.off()
  outputs[[name]] <<- data.frame(name=name,n_panels=length(ids),n_columns=ncol,n_rows=nrow,
    width_in=w,height_in=h,pdf=paste0("output/pdf/",name,".pdf"),png=paste0("output/png/",name,".png"))
}
for (id in c(cell_ids,marker_ids)) save_page(id,id)
save_page(cell_ids,"Six_study_celltype_UMAPs_Arial_square",3L,2L)
save_page(marker_ids[1:9],"Marker_UMAPs_3x3_page1_Chia_Fullard_Grubman",3L,3L)
save_page(marker_ids[10:18],"Marker_UMAPs_3x3_page2_Jakel_Su_Yao",3L,3L)
write.csv(do.call(rbind,coordinate_audit),file.path(task,"audit/coordinate_and_annotation_preservation.csv"),row.names=FALSE)
write.csv(do.call(rbind,scale_audit),file.path(task,"audit/marker_scale_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,label_audit),file.path(task,"audit/celltype_labels_and_anchors.csv"),row.names=FALSE)
write.csv(do.call(rbind,palette_audit),file.path(task,"audit/celltype_counts_and_colours.csv"),row.names=FALSE)
write.csv(do.call(rbind,dimensions),file.path(task,"audit/physical_panel_dimensions.csv"),row.names=FALSE)
write.csv(do.call(rbind,outputs),file.path(task,"audit/output_manifest.csv"),row.names=FALSE)
saveRDS(views,file.path(task,"audit/shared_square_view_limits.rds"))
writeLines(capture.output(sessionInfo()),file.path(task,"audit/plotting_sessionInfo.txt"))
message("UNIFORM_UMAP_PLOTTING_COMPLETE: 27 PDF and PNG files; 24 equal 114.3-mm square panels")

# S1D  Jakel cell-type UMAP
# Source: source/project/Reference/combine/results/source_celltype_umap_20260920/revisions/20261022_Arial_square_celltype_and_marker_UMAPs/code/plot_uniform_umaps.R
#!/usr/bin/env Rscript
# Display-only revision: exact selected inputs, selected marker bars, Arial,
# one shared square view per study and an identical physical data panel.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(scattermore)
  library(scales); library(grid); library(gtable); library(digest)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task <- normalizePath(args[[1]], mustWork = TRUE)
for (d in c("output/pdf", "output/png", "audit", "runtime"))
  dir.create(file.path(task, d), recursive = TRUE, showWarnings = FALSE)
studies <- c("Chia", "Fullard", "Grubman", "Jakel", "Su", "Yao")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
annotation <- c(Chia="cell.type", Fullard="celltype", Grubman="cellType",
                Jakel="CellType", Su="CellType", Yao="AuthorMajor")
expected_n <- c(Chia=87470L, Fullard=68557L, Grubman=13214L, Jakel=17799L,
                Su=224464L, Yao=73721L)
font <- systemfonts::match_fonts("Arial")
stopifnot(grepl("Arial", basename(font$path), fixed = TRUE))
writeLines(font$path, file.path(task, "audit/Arial_font_path.txt"))
base_colours <- read.csv(file.path(task,"inputs/source_tables/base_celltype_colours.csv"), stringsAsFactors=FALSE)
base_anchors <- read.csv(file.path(task,"inputs/source_tables/base_celltype_anchors.csv"), stringsAsFactors=FALSE)
yao_colours <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_colours.csv"), stringsAsFactors=FALSE)
yao_anchors <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_anchors.csv"), stringsAsFactors=FALSE)
bars <- read.csv(file.path(task,"audit/selected_marker_bar_provenance.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(bars)==18L, all(toupper(as.character(bars$exact_match))=="TRUE"),
          all(bars$upper_limit[bars$study %in% c("Grubman","Jakel")] == 1),
          all(bars$upper_limit[bars$study == "Fullard"] == 2))
theme_uniform <- theme_classic(base_size=10, base_family="Arial") + theme(
  text=element_text(family="Arial"),
  plot.title=element_text(face="bold",size=12,margin=margin(b=3)),
  plot.subtitle=element_text(size=8,colour="#333333",margin=margin(b=7)),
  axis.title=element_text(size=10), axis.text=element_text(size=9,colour="black"),
  axis.line=element_line(linewidth=.45), axis.ticks=element_line(linewidth=.35),
  legend.title=element_text(size=8),legend.text=element_text(size=8),
  legend.position="right",legend.box.margin=margin(0,0,0,5),
  plot.margin=margin(8,8,8,8),plot.title.position="panel")
plots <- list(); cell_ids <- character(); marker_ids <- character()
coordinate_audit <- list(); scale_audit <- list(); label_audit <- list()
palette_audit <- list(); views <- list()
for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Preparing ",study)
  ct <- readRDS(file.path(task,"inputs/celltype",paste0(study,".rds")))
  mb <- readRDS(file.path(task,"inputs/marker",paste0(study,".rds")))
  stopifnot(length(ct$cells)==expected_n[[study]], !anyDuplicated(ct$cells),
    identical(rownames(ct$umap),ct$cells),identical(rownames(ct$metadata),ct$cells),
    identical(rownames(mb$umap),mb$cells),identical(colnames(mb$expression),mb$cells),
    identical(colnames(mb$counts),mb$cells),setequal(ct$cells,mb$cells))
  align <- match(mb$cells,ct$cells)
  ct_aligned <- ct$umap[align,,drop=FALSE]
  # Matrix subsetting drops optional scaled:center attributes. Compare every
  # coordinate exactly after validating dimensions and cell alignment.
  stopifnot(identical(dim(ct_aligned),dim(mb$umap)),
    identical(as.numeric(ct_aligned),as.numeric(mb$umap)))
  coords <- ct$umap
  xr <- range(coords[,1]); yr <- range(coords[,2])
  span <- max(diff(xr),diff(yr))*1.16
  xlim <- mean(xr)+c(-.5,.5)*span; ylim <- mean(yr)+c(-.5,.5)*span
  views[[study]] <- list(xlim=xlim,ylim=ylim)
  values <- as.character(ct$metadata[[annotation[[study]]]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  if (study=="Yao") {
    pc <- subset(yao_colours,plot_id=="Yao_current_UMAP_author_labels")
    anchors <- subset(yao_anchors,plot_id=="Yao_current_UMAP_author_labels")
    names(pc)[names(pc)=="author_Major"] <- "celltype"
    names(anchors)[names(anchors)=="type"] <- "celltype"
  } else {
    pc <- base_colours[base_colours$plot_id==study,,drop=FALSE]
    anchors <- base_anchors[base_anchors$plot_id==study,,drop=FALSE]
  }
  types <- pc$celltype; palette <- setNames(pc$colour,types)
  counts <- table(factor(values,levels=types))
  stopifnot(setequal(unique(values),types),identical(as.integer(counts),as.integer(pc$n_cells)),
    identical(as.integer(counts[anchors$celltype]),as.integer(anchors$n_cells)),
    all(anchors$anchor_cell %in% ct$cells))
  anchor_pos <- match(anchors$anchor_cell,ct$cells)
  stopifnot(max(abs(coords[anchor_pos,1]-anchors$x))<1e-12,
            max(abs(coords[anchor_pos,2]-anchors$y))<1e-12)
  data <- data.frame(x=coords[,1],y=coords[,2],celltype=factor(values,levels=types),row.names=ct$cells)
  subtitle <- if (study=="Yao") "73,721 cells; 16 author cell types; 174 unmatched" else
    paste0(comma(length(values))," cells; ",length(types)," cell types")
  title <- if(study=="Yao") "Yao | author cell types" else paste(study,"| cell types")
  p <- ggplot(data,aes(x=x,y=y,colour=celltype))+
    geom_scattermore(pointsize=2.2,pixels=c(1400,1400),show.legend=FALSE)+
    scale_colour_manual(values=palette,drop=FALSE)+
    geom_text_repel(data=anchors,aes(x=x,y=y,label=label_text),inherit.aes=FALSE,
      family="Arial",colour="black",size=if(length(types)>20) 2.35 else 2.75,
      lineheight=.95,seed=20261016,nudge_x=anchors$nudge_x,nudge_y=anchors$nudge_y,
      box.padding=unit(.35,"lines"),point.padding=unit(.10,"lines"),
      min.segment.length=0,segment.colour="#454545",segment.size=.20,
      force=2,force_pull=.3,max.overlaps=Inf,max.time=Inf,max.iter=20000,
      show.legend=FALSE)+
    coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
    labs(title=title,subtitle=subtitle,x="UMAP_1",y="UMAP_2")+
    theme_uniform+theme(legend.position="none")
  id <- paste0(study,"_celltype_UMAP")
  plots[[id]] <- p; cell_ids <- c(cell_ids,id)
  anchors$study <- study
  label_audit[[study]] <- anchors[,c("study","celltype","n_cells","x","y","anchor_cell","label_text","nudge_x","nudge_y")]
  palette_audit[[study]] <- data.frame(study=study,celltype=types,n_cells=as.integer(counts),colour=pc$colour)
  coordinate_audit[[study]] <- data.frame(study=study,n_cells=length(values),annotation_column=annotation[[study]],
    n_celltype_labels=length(types),celltype_and_marker_barcode_sets_identical=TRUE,
    celltype_and_marker_coordinates_identical_by_barcode=TRUE,
    coordinates_sha256=digest(coords,algo="sha256"),
    celltype_values_sha256=digest(values,algo="sha256"),
    marker_expression_sha256=digest(mb$expression,algo="sha256"),
    original_x_min=xr[1],original_x_max=xr[2],original_y_min=yr[1],original_y_max=yr[2],
    view_x_min=xlim[1],view_x_max=xlim[2],view_y_min=ylim[1],view_y_max=ylim[2],
    coordinate_ratio=1,umap_recomputed=FALSE,annotation_modified=FALSE)
  thin <- CreateSeuratObject(counts=mb$counts,min.cells=0,min.features=0,project=study)
  thin <- SetAssayData(thin,assay="RNA",layer="data",new.data=mb$expression)
  mc <- mb$umap; colnames(mc) <- c("UMAP_1","UMAP_2")
  thin[["umap"]] <- CreateDimReducObject(embeddings=mc,key="UMAP_",assay="RNA")
  stopifnot(identical(dim(Embeddings(thin,"umap")),dim(mb$umap)),
    identical(as.numeric(Embeddings(thin,"umap")),as.numeric(mb$umap)),
    identical(as.matrix(GetAssayData(thin,assay="RNA",layer="data")),as.matrix(mb$expression)))
  for (gene in markers) {
    bar <- bars[bars$study==study & bars$marker==gene,,drop=FALSE]
    stopifnot(nrow(bar)==1L)
    feature <- mb$features[[gene]]; upper <- bar$upper_limit
    tick_values <- as.numeric(strsplit(bar$tick_values,";",fixed=TRUE)[[1]])
    tick_labels <- strsplit(bar$tick_labels,";",fixed=TRUE)[[1]]
    vals <- as.numeric(mb$expression[feature,])
    p <- FeaturePlot(thin,features=feature,reduction="umap",dims=c(1,2),
      cols=c("#D3D3D3","#E31A1C"),order=TRUE,pt.size=1.3,
      min.cutoff=0,max.cutoff=upper,slot="data",raster=TRUE,raster.dpi=c(1200,1200),combine=FALSE)[[1]]+
      scale_colour_gradient(low=bar$colour_low,high=bar$colour_high,
        limits=c(0,upper),oob=scales::squish,breaks=tick_values,labels=tick_labels,
        guide=guide_colourbar(title="Log-normalized\nRNA",barheight=unit(2.6,"cm"),barwidth=unit(.30,"cm")))+
      coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
      labs(title=paste(study,"|",gene),subtitle=paste0(comma(length(vals))," cells"),x="UMAP_1",y="UMAP_2")+
      theme_uniform
    # FeaturePlot must retain the cached values after applying only the chosen cap.
    j <- match(rownames(p$data),mb$cells)
    stopifnot(!anyNA(j),nrow(p$data)==length(mb$cells),!anyDuplicated(rownames(p$data)),
      setequal(rownames(p$data),mb$cells),
      identical(as.numeric(p$data[[feature]]),pmin(pmax(vals[j],0),upper)))
    id <- paste0(study,"_",gene,"_FeaturePlot")
    plots[[id]] <- p; marker_ids <- c(marker_ids,id)
    scale_audit[[id]] <- data.frame(study=study,marker=gene,lower_limit=0,upper_limit=upper,
      visible_breaks=bar$tick_values,visible_tick_labels=bar$tick_labels,
      colour_low=bar$colour_low,colour_high=bar$colour_high,n_cells=length(vals),
      n_positive=sum(vals>0),n_saturated=sum(vals>upper),selected_pdf_sha256=bar$sha256,
      expression_identical_before_display_cap=TRUE,sorted_ascending=TRUE,
      panel_width_in=4.5,panel_height_in=4.5,coordinate_ratio=1)
  }
  rm(ct,mb,thin); gc(verbose=FALSE)
}

# Measure on the actual Cairo/Arial PDF device. Fix the panel row and column,
# then equalize every surrounding row/column so no legend can shrink a panel.
measurement <- file.path(task,"runtime/layout_measurement.pdf")
cairo_pdf(measurement,width=9,height=9,family="Arial")
grobs <- lapply(plots,ggplotGrob)
panel_mm <- 114.3 # 4.5 inches, identical in all 24 panels
shape <- vapply(grobs,function(g) paste(length(g$widths),length(g$heights)),character(1))
stopifnot(length(unique(shape))==1L)
panel_locations <- vapply(grobs,function(g) {
  loc <- g$layout[g$layout$name=="panel",]
  paste(loc$l,loc$r,loc$t,loc$b,collapse=";")
},character(1))
stopifnot(length(unique(panel_locations))==1L)
dimensions <- list()
for (id in names(grobs)) {
  g <- grobs[[id]]; loc <- g$layout[g$layout$name=="panel",]
  stopifnot(nrow(loc)==1L,loc$l==loc$r,loc$t==loc$b)
  g$widths[loc$l] <- unit(panel_mm,"mm");g$heights[loc$t] <- unit(panel_mm,"mm")
  grobs[[id]] <- g
}
all_widths <- lapply(grobs,function(g) convertWidth(g$widths,"in",valueOnly=TRUE))
all_heights <- lapply(grobs,function(g) convertHeight(g$heights,"in",valueOnly=TRUE))
common_widths <- do.call(pmax,all_widths)
common_heights <- do.call(pmax,all_heights)
for (id in names(grobs)) {
  g <- grobs[[id]];g$widths <- unit(common_widths,"in");g$heights <- unit(common_heights,"in")
  grobs[[id]] <- g
  loc <- g$layout[g$layout$name=="panel",]
  pw <- convertWidth(g$widths[loc$l],"mm",valueOnly=TRUE)
  ph <- convertHeight(g$heights[loc$t],"mm",valueOnly=TRUE)
  stopifnot(abs(pw-panel_mm)<1e-9,abs(ph-panel_mm)<1e-9)
  dimensions[[id]] <- data.frame(plot_id=id,panel_width_mm=pw,panel_height_mm=ph,
    outer_width_in=sum(common_widths),outer_height_in=sum(common_heights),font="Arial")
}
dev.off();unlink(measurement)
slot_width <- sum(common_widths);slot_height <- sum(common_heights)
draw_page <- function(ids,ncol,nrow) {
  stopifnot(length(ids)==ncol*nrow)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol,
    widths=unit(rep(slot_width,ncol),"in"),heights=unit(rep(slot_height,nrow),"in"))))
  for (i in seq_along(ids)) {
    pushViewport(viewport(layout.pos.row=(i-1L)%/%ncol+1L,layout.pos.col=(i-1L)%%ncol+1L))
    grid.draw(grobs[[ids[[i]]]])
    upViewport()
  }
  upViewport()
}
outputs <- list()
save_page <- function(ids,name,ncol=1L,nrow=1L) {
  w <- slot_width*ncol;h <- slot_height*nrow
  message(format(Sys.time(),"%F %T")," Saving ",name," (",round(w,2)," x ",round(h,2)," in)")
  cairo_pdf(file.path(task,"output/pdf",paste0(name,".pdf")),width=w,height=h,family="Arial")
  draw_page(ids,ncol,nrow);dev.off()
  ragg::agg_png(file.path(task,"output/png",paste0(name,".png")),width=w,height=h,units="in",res=240,background="white")
  draw_page(ids,ncol,nrow);dev.off()
  outputs[[name]] <<- data.frame(name=name,n_panels=length(ids),n_columns=ncol,n_rows=nrow,
    width_in=w,height_in=h,pdf=paste0("output/pdf/",name,".pdf"),png=paste0("output/png/",name,".png"))
}
for (id in c(cell_ids,marker_ids)) save_page(id,id)
save_page(cell_ids,"Six_study_celltype_UMAPs_Arial_square",3L,2L)
save_page(marker_ids[1:9],"Marker_UMAPs_3x3_page1_Chia_Fullard_Grubman",3L,3L)
save_page(marker_ids[10:18],"Marker_UMAPs_3x3_page2_Jakel_Su_Yao",3L,3L)
write.csv(do.call(rbind,coordinate_audit),file.path(task,"audit/coordinate_and_annotation_preservation.csv"),row.names=FALSE)
write.csv(do.call(rbind,scale_audit),file.path(task,"audit/marker_scale_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,label_audit),file.path(task,"audit/celltype_labels_and_anchors.csv"),row.names=FALSE)
write.csv(do.call(rbind,palette_audit),file.path(task,"audit/celltype_counts_and_colours.csv"),row.names=FALSE)
write.csv(do.call(rbind,dimensions),file.path(task,"audit/physical_panel_dimensions.csv"),row.names=FALSE)
write.csv(do.call(rbind,outputs),file.path(task,"audit/output_manifest.csv"),row.names=FALSE)
saveRDS(views,file.path(task,"audit/shared_square_view_limits.rds"))
writeLines(capture.output(sessionInfo()),file.path(task,"audit/plotting_sessionInfo.txt"))
message("UNIFORM_UMAP_PLOTTING_COMPLETE: 27 PDF and PNG files; 24 equal 114.3-mm square panels")

# S1E  Su cell-type UMAP
# Source: source/project/Reference/combine/results/source_celltype_umap_20260920/revisions/20261022_Arial_square_celltype_and_marker_UMAPs/code/plot_uniform_umaps.R
#!/usr/bin/env Rscript
# Display-only revision: exact selected inputs, selected marker bars, Arial,
# one shared square view per study and an identical physical data panel.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(scattermore)
  library(scales); library(grid); library(gtable); library(digest)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task <- normalizePath(args[[1]], mustWork = TRUE)
for (d in c("output/pdf", "output/png", "audit", "runtime"))
  dir.create(file.path(task, d), recursive = TRUE, showWarnings = FALSE)
studies <- c("Chia", "Fullard", "Grubman", "Jakel", "Su", "Yao")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
annotation <- c(Chia="cell.type", Fullard="celltype", Grubman="cellType",
                Jakel="CellType", Su="CellType", Yao="AuthorMajor")
expected_n <- c(Chia=87470L, Fullard=68557L, Grubman=13214L, Jakel=17799L,
                Su=224464L, Yao=73721L)
font <- systemfonts::match_fonts("Arial")
stopifnot(grepl("Arial", basename(font$path), fixed = TRUE))
writeLines(font$path, file.path(task, "audit/Arial_font_path.txt"))
base_colours <- read.csv(file.path(task,"inputs/source_tables/base_celltype_colours.csv"), stringsAsFactors=FALSE)
base_anchors <- read.csv(file.path(task,"inputs/source_tables/base_celltype_anchors.csv"), stringsAsFactors=FALSE)
yao_colours <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_colours.csv"), stringsAsFactors=FALSE)
yao_anchors <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_anchors.csv"), stringsAsFactors=FALSE)
bars <- read.csv(file.path(task,"audit/selected_marker_bar_provenance.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(bars)==18L, all(toupper(as.character(bars$exact_match))=="TRUE"),
          all(bars$upper_limit[bars$study %in% c("Grubman","Jakel")] == 1),
          all(bars$upper_limit[bars$study == "Fullard"] == 2))
theme_uniform <- theme_classic(base_size=10, base_family="Arial") + theme(
  text=element_text(family="Arial"),
  plot.title=element_text(face="bold",size=12,margin=margin(b=3)),
  plot.subtitle=element_text(size=8,colour="#333333",margin=margin(b=7)),
  axis.title=element_text(size=10), axis.text=element_text(size=9,colour="black"),
  axis.line=element_line(linewidth=.45), axis.ticks=element_line(linewidth=.35),
  legend.title=element_text(size=8),legend.text=element_text(size=8),
  legend.position="right",legend.box.margin=margin(0,0,0,5),
  plot.margin=margin(8,8,8,8),plot.title.position="panel")
plots <- list(); cell_ids <- character(); marker_ids <- character()
coordinate_audit <- list(); scale_audit <- list(); label_audit <- list()
palette_audit <- list(); views <- list()
for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Preparing ",study)
  ct <- readRDS(file.path(task,"inputs/celltype",paste0(study,".rds")))
  mb <- readRDS(file.path(task,"inputs/marker",paste0(study,".rds")))
  stopifnot(length(ct$cells)==expected_n[[study]], !anyDuplicated(ct$cells),
    identical(rownames(ct$umap),ct$cells),identical(rownames(ct$metadata),ct$cells),
    identical(rownames(mb$umap),mb$cells),identical(colnames(mb$expression),mb$cells),
    identical(colnames(mb$counts),mb$cells),setequal(ct$cells,mb$cells))
  align <- match(mb$cells,ct$cells)
  ct_aligned <- ct$umap[align,,drop=FALSE]
  # Matrix subsetting drops optional scaled:center attributes. Compare every
  # coordinate exactly after validating dimensions and cell alignment.
  stopifnot(identical(dim(ct_aligned),dim(mb$umap)),
    identical(as.numeric(ct_aligned),as.numeric(mb$umap)))
  coords <- ct$umap
  xr <- range(coords[,1]); yr <- range(coords[,2])
  span <- max(diff(xr),diff(yr))*1.16
  xlim <- mean(xr)+c(-.5,.5)*span; ylim <- mean(yr)+c(-.5,.5)*span
  views[[study]] <- list(xlim=xlim,ylim=ylim)
  values <- as.character(ct$metadata[[annotation[[study]]]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  if (study=="Yao") {
    pc <- subset(yao_colours,plot_id=="Yao_current_UMAP_author_labels")
    anchors <- subset(yao_anchors,plot_id=="Yao_current_UMAP_author_labels")
    names(pc)[names(pc)=="author_Major"] <- "celltype"
    names(anchors)[names(anchors)=="type"] <- "celltype"
  } else {
    pc <- base_colours[base_colours$plot_id==study,,drop=FALSE]
    anchors <- base_anchors[base_anchors$plot_id==study,,drop=FALSE]
  }
  types <- pc$celltype; palette <- setNames(pc$colour,types)
  counts <- table(factor(values,levels=types))
  stopifnot(setequal(unique(values),types),identical(as.integer(counts),as.integer(pc$n_cells)),
    identical(as.integer(counts[anchors$celltype]),as.integer(anchors$n_cells)),
    all(anchors$anchor_cell %in% ct$cells))
  anchor_pos <- match(anchors$anchor_cell,ct$cells)
  stopifnot(max(abs(coords[anchor_pos,1]-anchors$x))<1e-12,
            max(abs(coords[anchor_pos,2]-anchors$y))<1e-12)
  data <- data.frame(x=coords[,1],y=coords[,2],celltype=factor(values,levels=types),row.names=ct$cells)
  subtitle <- if (study=="Yao") "73,721 cells; 16 author cell types; 174 unmatched" else
    paste0(comma(length(values))," cells; ",length(types)," cell types")
  title <- if(study=="Yao") "Yao | author cell types" else paste(study,"| cell types")
  p <- ggplot(data,aes(x=x,y=y,colour=celltype))+
    geom_scattermore(pointsize=2.2,pixels=c(1400,1400),show.legend=FALSE)+
    scale_colour_manual(values=palette,drop=FALSE)+
    geom_text_repel(data=anchors,aes(x=x,y=y,label=label_text),inherit.aes=FALSE,
      family="Arial",colour="black",size=if(length(types)>20) 2.35 else 2.75,
      lineheight=.95,seed=20261016,nudge_x=anchors$nudge_x,nudge_y=anchors$nudge_y,
      box.padding=unit(.35,"lines"),point.padding=unit(.10,"lines"),
      min.segment.length=0,segment.colour="#454545",segment.size=.20,
      force=2,force_pull=.3,max.overlaps=Inf,max.time=Inf,max.iter=20000,
      show.legend=FALSE)+
    coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
    labs(title=title,subtitle=subtitle,x="UMAP_1",y="UMAP_2")+
    theme_uniform+theme(legend.position="none")
  id <- paste0(study,"_celltype_UMAP")
  plots[[id]] <- p; cell_ids <- c(cell_ids,id)
  anchors$study <- study
  label_audit[[study]] <- anchors[,c("study","celltype","n_cells","x","y","anchor_cell","label_text","nudge_x","nudge_y")]
  palette_audit[[study]] <- data.frame(study=study,celltype=types,n_cells=as.integer(counts),colour=pc$colour)
  coordinate_audit[[study]] <- data.frame(study=study,n_cells=length(values),annotation_column=annotation[[study]],
    n_celltype_labels=length(types),celltype_and_marker_barcode_sets_identical=TRUE,
    celltype_and_marker_coordinates_identical_by_barcode=TRUE,
    coordinates_sha256=digest(coords,algo="sha256"),
    celltype_values_sha256=digest(values,algo="sha256"),
    marker_expression_sha256=digest(mb$expression,algo="sha256"),
    original_x_min=xr[1],original_x_max=xr[2],original_y_min=yr[1],original_y_max=yr[2],
    view_x_min=xlim[1],view_x_max=xlim[2],view_y_min=ylim[1],view_y_max=ylim[2],
    coordinate_ratio=1,umap_recomputed=FALSE,annotation_modified=FALSE)
  thin <- CreateSeuratObject(counts=mb$counts,min.cells=0,min.features=0,project=study)
  thin <- SetAssayData(thin,assay="RNA",layer="data",new.data=mb$expression)
  mc <- mb$umap; colnames(mc) <- c("UMAP_1","UMAP_2")
  thin[["umap"]] <- CreateDimReducObject(embeddings=mc,key="UMAP_",assay="RNA")
  stopifnot(identical(dim(Embeddings(thin,"umap")),dim(mb$umap)),
    identical(as.numeric(Embeddings(thin,"umap")),as.numeric(mb$umap)),
    identical(as.matrix(GetAssayData(thin,assay="RNA",layer="data")),as.matrix(mb$expression)))
  for (gene in markers) {
    bar <- bars[bars$study==study & bars$marker==gene,,drop=FALSE]
    stopifnot(nrow(bar)==1L)
    feature <- mb$features[[gene]]; upper <- bar$upper_limit
    tick_values <- as.numeric(strsplit(bar$tick_values,";",fixed=TRUE)[[1]])
    tick_labels <- strsplit(bar$tick_labels,";",fixed=TRUE)[[1]]
    vals <- as.numeric(mb$expression[feature,])
    p <- FeaturePlot(thin,features=feature,reduction="umap",dims=c(1,2),
      cols=c("#D3D3D3","#E31A1C"),order=TRUE,pt.size=1.3,
      min.cutoff=0,max.cutoff=upper,slot="data",raster=TRUE,raster.dpi=c(1200,1200),combine=FALSE)[[1]]+
      scale_colour_gradient(low=bar$colour_low,high=bar$colour_high,
        limits=c(0,upper),oob=scales::squish,breaks=tick_values,labels=tick_labels,
        guide=guide_colourbar(title="Log-normalized\nRNA",barheight=unit(2.6,"cm"),barwidth=unit(.30,"cm")))+
      coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
      labs(title=paste(study,"|",gene),subtitle=paste0(comma(length(vals))," cells"),x="UMAP_1",y="UMAP_2")+
      theme_uniform
    # FeaturePlot must retain the cached values after applying only the chosen cap.
    j <- match(rownames(p$data),mb$cells)
    stopifnot(!anyNA(j),nrow(p$data)==length(mb$cells),!anyDuplicated(rownames(p$data)),
      setequal(rownames(p$data),mb$cells),
      identical(as.numeric(p$data[[feature]]),pmin(pmax(vals[j],0),upper)))
    id <- paste0(study,"_",gene,"_FeaturePlot")
    plots[[id]] <- p; marker_ids <- c(marker_ids,id)
    scale_audit[[id]] <- data.frame(study=study,marker=gene,lower_limit=0,upper_limit=upper,
      visible_breaks=bar$tick_values,visible_tick_labels=bar$tick_labels,
      colour_low=bar$colour_low,colour_high=bar$colour_high,n_cells=length(vals),
      n_positive=sum(vals>0),n_saturated=sum(vals>upper),selected_pdf_sha256=bar$sha256,
      expression_identical_before_display_cap=TRUE,sorted_ascending=TRUE,
      panel_width_in=4.5,panel_height_in=4.5,coordinate_ratio=1)
  }
  rm(ct,mb,thin); gc(verbose=FALSE)
}

# Measure on the actual Cairo/Arial PDF device. Fix the panel row and column,
# then equalize every surrounding row/column so no legend can shrink a panel.
measurement <- file.path(task,"runtime/layout_measurement.pdf")
cairo_pdf(measurement,width=9,height=9,family="Arial")
grobs <- lapply(plots,ggplotGrob)
panel_mm <- 114.3 # 4.5 inches, identical in all 24 panels
shape <- vapply(grobs,function(g) paste(length(g$widths),length(g$heights)),character(1))
stopifnot(length(unique(shape))==1L)
panel_locations <- vapply(grobs,function(g) {
  loc <- g$layout[g$layout$name=="panel",]
  paste(loc$l,loc$r,loc$t,loc$b,collapse=";")
},character(1))
stopifnot(length(unique(panel_locations))==1L)
dimensions <- list()
for (id in names(grobs)) {
  g <- grobs[[id]]; loc <- g$layout[g$layout$name=="panel",]
  stopifnot(nrow(loc)==1L,loc$l==loc$r,loc$t==loc$b)
  g$widths[loc$l] <- unit(panel_mm,"mm");g$heights[loc$t] <- unit(panel_mm,"mm")
  grobs[[id]] <- g
}
all_widths <- lapply(grobs,function(g) convertWidth(g$widths,"in",valueOnly=TRUE))
all_heights <- lapply(grobs,function(g) convertHeight(g$heights,"in",valueOnly=TRUE))
common_widths <- do.call(pmax,all_widths)
common_heights <- do.call(pmax,all_heights)
for (id in names(grobs)) {
  g <- grobs[[id]];g$widths <- unit(common_widths,"in");g$heights <- unit(common_heights,"in")
  grobs[[id]] <- g
  loc <- g$layout[g$layout$name=="panel",]
  pw <- convertWidth(g$widths[loc$l],"mm",valueOnly=TRUE)
  ph <- convertHeight(g$heights[loc$t],"mm",valueOnly=TRUE)
  stopifnot(abs(pw-panel_mm)<1e-9,abs(ph-panel_mm)<1e-9)
  dimensions[[id]] <- data.frame(plot_id=id,panel_width_mm=pw,panel_height_mm=ph,
    outer_width_in=sum(common_widths),outer_height_in=sum(common_heights),font="Arial")
}
dev.off();unlink(measurement)
slot_width <- sum(common_widths);slot_height <- sum(common_heights)
draw_page <- function(ids,ncol,nrow) {
  stopifnot(length(ids)==ncol*nrow)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol,
    widths=unit(rep(slot_width,ncol),"in"),heights=unit(rep(slot_height,nrow),"in"))))
  for (i in seq_along(ids)) {
    pushViewport(viewport(layout.pos.row=(i-1L)%/%ncol+1L,layout.pos.col=(i-1L)%%ncol+1L))
    grid.draw(grobs[[ids[[i]]]])
    upViewport()
  }
  upViewport()
}
outputs <- list()
save_page <- function(ids,name,ncol=1L,nrow=1L) {
  w <- slot_width*ncol;h <- slot_height*nrow
  message(format(Sys.time(),"%F %T")," Saving ",name," (",round(w,2)," x ",round(h,2)," in)")
  cairo_pdf(file.path(task,"output/pdf",paste0(name,".pdf")),width=w,height=h,family="Arial")
  draw_page(ids,ncol,nrow);dev.off()
  ragg::agg_png(file.path(task,"output/png",paste0(name,".png")),width=w,height=h,units="in",res=240,background="white")
  draw_page(ids,ncol,nrow);dev.off()
  outputs[[name]] <<- data.frame(name=name,n_panels=length(ids),n_columns=ncol,n_rows=nrow,
    width_in=w,height_in=h,pdf=paste0("output/pdf/",name,".pdf"),png=paste0("output/png/",name,".png"))
}
for (id in c(cell_ids,marker_ids)) save_page(id,id)
save_page(cell_ids,"Six_study_celltype_UMAPs_Arial_square",3L,2L)
save_page(marker_ids[1:9],"Marker_UMAPs_3x3_page1_Chia_Fullard_Grubman",3L,3L)
save_page(marker_ids[10:18],"Marker_UMAPs_3x3_page2_Jakel_Su_Yao",3L,3L)
write.csv(do.call(rbind,coordinate_audit),file.path(task,"audit/coordinate_and_annotation_preservation.csv"),row.names=FALSE)
write.csv(do.call(rbind,scale_audit),file.path(task,"audit/marker_scale_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,label_audit),file.path(task,"audit/celltype_labels_and_anchors.csv"),row.names=FALSE)
write.csv(do.call(rbind,palette_audit),file.path(task,"audit/celltype_counts_and_colours.csv"),row.names=FALSE)
write.csv(do.call(rbind,dimensions),file.path(task,"audit/physical_panel_dimensions.csv"),row.names=FALSE)
write.csv(do.call(rbind,outputs),file.path(task,"audit/output_manifest.csv"),row.names=FALSE)
saveRDS(views,file.path(task,"audit/shared_square_view_limits.rds"))
writeLines(capture.output(sessionInfo()),file.path(task,"audit/plotting_sessionInfo.txt"))
message("UNIFORM_UMAP_PLOTTING_COMPLETE: 27 PDF and PNG files; 24 equal 114.3-mm square panels")

# S1F  Yao cell-type UMAP
# Source: source/project/Reference/combine/results/source_celltype_umap_20260920/revisions/20261022_Arial_square_celltype_and_marker_UMAPs/code/plot_uniform_umaps.R
#!/usr/bin/env Rscript
# Display-only revision: exact selected inputs, selected marker bars, Arial,
# one shared square view per study and an identical physical data panel.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(scattermore)
  library(scales); library(grid); library(gtable); library(digest)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task <- normalizePath(args[[1]], mustWork = TRUE)
for (d in c("output/pdf", "output/png", "audit", "runtime"))
  dir.create(file.path(task, d), recursive = TRUE, showWarnings = FALSE)
studies <- c("Chia", "Fullard", "Grubman", "Jakel", "Su", "Yao")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
annotation <- c(Chia="cell.type", Fullard="celltype", Grubman="cellType",
                Jakel="CellType", Su="CellType", Yao="AuthorMajor")
expected_n <- c(Chia=87470L, Fullard=68557L, Grubman=13214L, Jakel=17799L,
                Su=224464L, Yao=73721L)
font <- systemfonts::match_fonts("Arial")
stopifnot(grepl("Arial", basename(font$path), fixed = TRUE))
writeLines(font$path, file.path(task, "audit/Arial_font_path.txt"))
base_colours <- read.csv(file.path(task,"inputs/source_tables/base_celltype_colours.csv"), stringsAsFactors=FALSE)
base_anchors <- read.csv(file.path(task,"inputs/source_tables/base_celltype_anchors.csv"), stringsAsFactors=FALSE)
yao_colours <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_colours.csv"), stringsAsFactors=FALSE)
yao_anchors <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_anchors.csv"), stringsAsFactors=FALSE)
bars <- read.csv(file.path(task,"audit/selected_marker_bar_provenance.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(bars)==18L, all(toupper(as.character(bars$exact_match))=="TRUE"),
          all(bars$upper_limit[bars$study %in% c("Grubman","Jakel")] == 1),
          all(bars$upper_limit[bars$study == "Fullard"] == 2))
theme_uniform <- theme_classic(base_size=10, base_family="Arial") + theme(
  text=element_text(family="Arial"),
  plot.title=element_text(face="bold",size=12,margin=margin(b=3)),
  plot.subtitle=element_text(size=8,colour="#333333",margin=margin(b=7)),
  axis.title=element_text(size=10), axis.text=element_text(size=9,colour="black"),
  axis.line=element_line(linewidth=.45), axis.ticks=element_line(linewidth=.35),
  legend.title=element_text(size=8),legend.text=element_text(size=8),
  legend.position="right",legend.box.margin=margin(0,0,0,5),
  plot.margin=margin(8,8,8,8),plot.title.position="panel")
plots <- list(); cell_ids <- character(); marker_ids <- character()
coordinate_audit <- list(); scale_audit <- list(); label_audit <- list()
palette_audit <- list(); views <- list()
for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Preparing ",study)
  ct <- readRDS(file.path(task,"inputs/celltype",paste0(study,".rds")))
  mb <- readRDS(file.path(task,"inputs/marker",paste0(study,".rds")))
  stopifnot(length(ct$cells)==expected_n[[study]], !anyDuplicated(ct$cells),
    identical(rownames(ct$umap),ct$cells),identical(rownames(ct$metadata),ct$cells),
    identical(rownames(mb$umap),mb$cells),identical(colnames(mb$expression),mb$cells),
    identical(colnames(mb$counts),mb$cells),setequal(ct$cells,mb$cells))
  align <- match(mb$cells,ct$cells)
  ct_aligned <- ct$umap[align,,drop=FALSE]
  # Matrix subsetting drops optional scaled:center attributes. Compare every
  # coordinate exactly after validating dimensions and cell alignment.
  stopifnot(identical(dim(ct_aligned),dim(mb$umap)),
    identical(as.numeric(ct_aligned),as.numeric(mb$umap)))
  coords <- ct$umap
  xr <- range(coords[,1]); yr <- range(coords[,2])
  span <- max(diff(xr),diff(yr))*1.16
  xlim <- mean(xr)+c(-.5,.5)*span; ylim <- mean(yr)+c(-.5,.5)*span
  views[[study]] <- list(xlim=xlim,ylim=ylim)
  values <- as.character(ct$metadata[[annotation[[study]]]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  if (study=="Yao") {
    pc <- subset(yao_colours,plot_id=="Yao_current_UMAP_author_labels")
    anchors <- subset(yao_anchors,plot_id=="Yao_current_UMAP_author_labels")
    names(pc)[names(pc)=="author_Major"] <- "celltype"
    names(anchors)[names(anchors)=="type"] <- "celltype"
  } else {
    pc <- base_colours[base_colours$plot_id==study,,drop=FALSE]
    anchors <- base_anchors[base_anchors$plot_id==study,,drop=FALSE]
  }
  types <- pc$celltype; palette <- setNames(pc$colour,types)
  counts <- table(factor(values,levels=types))
  stopifnot(setequal(unique(values),types),identical(as.integer(counts),as.integer(pc$n_cells)),
    identical(as.integer(counts[anchors$celltype]),as.integer(anchors$n_cells)),
    all(anchors$anchor_cell %in% ct$cells))
  anchor_pos <- match(anchors$anchor_cell,ct$cells)
  stopifnot(max(abs(coords[anchor_pos,1]-anchors$x))<1e-12,
            max(abs(coords[anchor_pos,2]-anchors$y))<1e-12)
  data <- data.frame(x=coords[,1],y=coords[,2],celltype=factor(values,levels=types),row.names=ct$cells)
  subtitle <- if (study=="Yao") "73,721 cells; 16 author cell types; 174 unmatched" else
    paste0(comma(length(values))," cells; ",length(types)," cell types")
  title <- if(study=="Yao") "Yao | author cell types" else paste(study,"| cell types")
  p <- ggplot(data,aes(x=x,y=y,colour=celltype))+
    geom_scattermore(pointsize=2.2,pixels=c(1400,1400),show.legend=FALSE)+
    scale_colour_manual(values=palette,drop=FALSE)+
    geom_text_repel(data=anchors,aes(x=x,y=y,label=label_text),inherit.aes=FALSE,
      family="Arial",colour="black",size=if(length(types)>20) 2.35 else 2.75,
      lineheight=.95,seed=20261016,nudge_x=anchors$nudge_x,nudge_y=anchors$nudge_y,
      box.padding=unit(.35,"lines"),point.padding=unit(.10,"lines"),
      min.segment.length=0,segment.colour="#454545",segment.size=.20,
      force=2,force_pull=.3,max.overlaps=Inf,max.time=Inf,max.iter=20000,
      show.legend=FALSE)+
    coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
    labs(title=title,subtitle=subtitle,x="UMAP_1",y="UMAP_2")+
    theme_uniform+theme(legend.position="none")
  id <- paste0(study,"_celltype_UMAP")
  plots[[id]] <- p; cell_ids <- c(cell_ids,id)
  anchors$study <- study
  label_audit[[study]] <- anchors[,c("study","celltype","n_cells","x","y","anchor_cell","label_text","nudge_x","nudge_y")]
  palette_audit[[study]] <- data.frame(study=study,celltype=types,n_cells=as.integer(counts),colour=pc$colour)
  coordinate_audit[[study]] <- data.frame(study=study,n_cells=length(values),annotation_column=annotation[[study]],
    n_celltype_labels=length(types),celltype_and_marker_barcode_sets_identical=TRUE,
    celltype_and_marker_coordinates_identical_by_barcode=TRUE,
    coordinates_sha256=digest(coords,algo="sha256"),
    celltype_values_sha256=digest(values,algo="sha256"),
    marker_expression_sha256=digest(mb$expression,algo="sha256"),
    original_x_min=xr[1],original_x_max=xr[2],original_y_min=yr[1],original_y_max=yr[2],
    view_x_min=xlim[1],view_x_max=xlim[2],view_y_min=ylim[1],view_y_max=ylim[2],
    coordinate_ratio=1,umap_recomputed=FALSE,annotation_modified=FALSE)
  thin <- CreateSeuratObject(counts=mb$counts,min.cells=0,min.features=0,project=study)
  thin <- SetAssayData(thin,assay="RNA",layer="data",new.data=mb$expression)
  mc <- mb$umap; colnames(mc) <- c("UMAP_1","UMAP_2")
  thin[["umap"]] <- CreateDimReducObject(embeddings=mc,key="UMAP_",assay="RNA")
  stopifnot(identical(dim(Embeddings(thin,"umap")),dim(mb$umap)),
    identical(as.numeric(Embeddings(thin,"umap")),as.numeric(mb$umap)),
    identical(as.matrix(GetAssayData(thin,assay="RNA",layer="data")),as.matrix(mb$expression)))
  for (gene in markers) {
    bar <- bars[bars$study==study & bars$marker==gene,,drop=FALSE]
    stopifnot(nrow(bar)==1L)
    feature <- mb$features[[gene]]; upper <- bar$upper_limit
    tick_values <- as.numeric(strsplit(bar$tick_values,";",fixed=TRUE)[[1]])
    tick_labels <- strsplit(bar$tick_labels,";",fixed=TRUE)[[1]]
    vals <- as.numeric(mb$expression[feature,])
    p <- FeaturePlot(thin,features=feature,reduction="umap",dims=c(1,2),
      cols=c("#D3D3D3","#E31A1C"),order=TRUE,pt.size=1.3,
      min.cutoff=0,max.cutoff=upper,slot="data",raster=TRUE,raster.dpi=c(1200,1200),combine=FALSE)[[1]]+
      scale_colour_gradient(low=bar$colour_low,high=bar$colour_high,
        limits=c(0,upper),oob=scales::squish,breaks=tick_values,labels=tick_labels,
        guide=guide_colourbar(title="Log-normalized\nRNA",barheight=unit(2.6,"cm"),barwidth=unit(.30,"cm")))+
      coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
      labs(title=paste(study,"|",gene),subtitle=paste0(comma(length(vals))," cells"),x="UMAP_1",y="UMAP_2")+
      theme_uniform
    # FeaturePlot must retain the cached values after applying only the chosen cap.
    j <- match(rownames(p$data),mb$cells)
    stopifnot(!anyNA(j),nrow(p$data)==length(mb$cells),!anyDuplicated(rownames(p$data)),
      setequal(rownames(p$data),mb$cells),
      identical(as.numeric(p$data[[feature]]),pmin(pmax(vals[j],0),upper)))
    id <- paste0(study,"_",gene,"_FeaturePlot")
    plots[[id]] <- p; marker_ids <- c(marker_ids,id)
    scale_audit[[id]] <- data.frame(study=study,marker=gene,lower_limit=0,upper_limit=upper,
      visible_breaks=bar$tick_values,visible_tick_labels=bar$tick_labels,
      colour_low=bar$colour_low,colour_high=bar$colour_high,n_cells=length(vals),
      n_positive=sum(vals>0),n_saturated=sum(vals>upper),selected_pdf_sha256=bar$sha256,
      expression_identical_before_display_cap=TRUE,sorted_ascending=TRUE,
      panel_width_in=4.5,panel_height_in=4.5,coordinate_ratio=1)
  }
  rm(ct,mb,thin); gc(verbose=FALSE)
}

# Measure on the actual Cairo/Arial PDF device. Fix the panel row and column,
# then equalize every surrounding row/column so no legend can shrink a panel.
measurement <- file.path(task,"runtime/layout_measurement.pdf")
cairo_pdf(measurement,width=9,height=9,family="Arial")
grobs <- lapply(plots,ggplotGrob)
panel_mm <- 114.3 # 4.5 inches, identical in all 24 panels
shape <- vapply(grobs,function(g) paste(length(g$widths),length(g$heights)),character(1))
stopifnot(length(unique(shape))==1L)
panel_locations <- vapply(grobs,function(g) {
  loc <- g$layout[g$layout$name=="panel",]
  paste(loc$l,loc$r,loc$t,loc$b,collapse=";")
},character(1))
stopifnot(length(unique(panel_locations))==1L)
dimensions <- list()
for (id in names(grobs)) {
  g <- grobs[[id]]; loc <- g$layout[g$layout$name=="panel",]
  stopifnot(nrow(loc)==1L,loc$l==loc$r,loc$t==loc$b)
  g$widths[loc$l] <- unit(panel_mm,"mm");g$heights[loc$t] <- unit(panel_mm,"mm")
  grobs[[id]] <- g
}
all_widths <- lapply(grobs,function(g) convertWidth(g$widths,"in",valueOnly=TRUE))
all_heights <- lapply(grobs,function(g) convertHeight(g$heights,"in",valueOnly=TRUE))
common_widths <- do.call(pmax,all_widths)
common_heights <- do.call(pmax,all_heights)
for (id in names(grobs)) {
  g <- grobs[[id]];g$widths <- unit(common_widths,"in");g$heights <- unit(common_heights,"in")
  grobs[[id]] <- g
  loc <- g$layout[g$layout$name=="panel",]
  pw <- convertWidth(g$widths[loc$l],"mm",valueOnly=TRUE)
  ph <- convertHeight(g$heights[loc$t],"mm",valueOnly=TRUE)
  stopifnot(abs(pw-panel_mm)<1e-9,abs(ph-panel_mm)<1e-9)
  dimensions[[id]] <- data.frame(plot_id=id,panel_width_mm=pw,panel_height_mm=ph,
    outer_width_in=sum(common_widths),outer_height_in=sum(common_heights),font="Arial")
}
dev.off();unlink(measurement)
slot_width <- sum(common_widths);slot_height <- sum(common_heights)
draw_page <- function(ids,ncol,nrow) {
  stopifnot(length(ids)==ncol*nrow)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol,
    widths=unit(rep(slot_width,ncol),"in"),heights=unit(rep(slot_height,nrow),"in"))))
  for (i in seq_along(ids)) {
    pushViewport(viewport(layout.pos.row=(i-1L)%/%ncol+1L,layout.pos.col=(i-1L)%%ncol+1L))
    grid.draw(grobs[[ids[[i]]]])
    upViewport()
  }
  upViewport()
}
outputs <- list()
save_page <- function(ids,name,ncol=1L,nrow=1L) {
  w <- slot_width*ncol;h <- slot_height*nrow
  message(format(Sys.time(),"%F %T")," Saving ",name," (",round(w,2)," x ",round(h,2)," in)")
  cairo_pdf(file.path(task,"output/pdf",paste0(name,".pdf")),width=w,height=h,family="Arial")
  draw_page(ids,ncol,nrow);dev.off()
  ragg::agg_png(file.path(task,"output/png",paste0(name,".png")),width=w,height=h,units="in",res=240,background="white")
  draw_page(ids,ncol,nrow);dev.off()
  outputs[[name]] <<- data.frame(name=name,n_panels=length(ids),n_columns=ncol,n_rows=nrow,
    width_in=w,height_in=h,pdf=paste0("output/pdf/",name,".pdf"),png=paste0("output/png/",name,".png"))
}
for (id in c(cell_ids,marker_ids)) save_page(id,id)
save_page(cell_ids,"Six_study_celltype_UMAPs_Arial_square",3L,2L)
save_page(marker_ids[1:9],"Marker_UMAPs_3x3_page1_Chia_Fullard_Grubman",3L,3L)
save_page(marker_ids[10:18],"Marker_UMAPs_3x3_page2_Jakel_Su_Yao",3L,3L)
write.csv(do.call(rbind,coordinate_audit),file.path(task,"audit/coordinate_and_annotation_preservation.csv"),row.names=FALSE)
write.csv(do.call(rbind,scale_audit),file.path(task,"audit/marker_scale_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,label_audit),file.path(task,"audit/celltype_labels_and_anchors.csv"),row.names=FALSE)
write.csv(do.call(rbind,palette_audit),file.path(task,"audit/celltype_counts_and_colours.csv"),row.names=FALSE)
write.csv(do.call(rbind,dimensions),file.path(task,"audit/physical_panel_dimensions.csv"),row.names=FALSE)
write.csv(do.call(rbind,outputs),file.path(task,"audit/output_manifest.csv"),row.names=FALSE)
saveRDS(views,file.path(task,"audit/shared_square_view_limits.rds"))
writeLines(capture.output(sessionInfo()),file.path(task,"audit/plotting_sessionInfo.txt"))
message("UNIFORM_UMAP_PLOTTING_COMPLETE: 27 PDF and PNG files; 24 equal 114.3-mm square panels")


# Supplementary Figure 2 ------------------------------------------------------------

# S2A  Chia P2RY12 / CX3CR1 / CSF1R UMAPs
# Source: source/project/Reference/combine/results/source_celltype_umap_20260920/revisions/20261022_Arial_square_celltype_and_marker_UMAPs/code/plot_uniform_umaps.R
#!/usr/bin/env Rscript
# Display-only revision: exact selected inputs, selected marker bars, Arial,
# one shared square view per study and an identical physical data panel.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(scattermore)
  library(scales); library(grid); library(gtable); library(digest)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task <- normalizePath(args[[1]], mustWork = TRUE)
for (d in c("output/pdf", "output/png", "audit", "runtime"))
  dir.create(file.path(task, d), recursive = TRUE, showWarnings = FALSE)
studies <- c("Chia", "Fullard", "Grubman", "Jakel", "Su", "Yao")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
annotation <- c(Chia="cell.type", Fullard="celltype", Grubman="cellType",
                Jakel="CellType", Su="CellType", Yao="AuthorMajor")
expected_n <- c(Chia=87470L, Fullard=68557L, Grubman=13214L, Jakel=17799L,
                Su=224464L, Yao=73721L)
font <- systemfonts::match_fonts("Arial")
stopifnot(grepl("Arial", basename(font$path), fixed = TRUE))
writeLines(font$path, file.path(task, "audit/Arial_font_path.txt"))
base_colours <- read.csv(file.path(task,"inputs/source_tables/base_celltype_colours.csv"), stringsAsFactors=FALSE)
base_anchors <- read.csv(file.path(task,"inputs/source_tables/base_celltype_anchors.csv"), stringsAsFactors=FALSE)
yao_colours <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_colours.csv"), stringsAsFactors=FALSE)
yao_anchors <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_anchors.csv"), stringsAsFactors=FALSE)
bars <- read.csv(file.path(task,"audit/selected_marker_bar_provenance.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(bars)==18L, all(toupper(as.character(bars$exact_match))=="TRUE"),
          all(bars$upper_limit[bars$study %in% c("Grubman","Jakel")] == 1),
          all(bars$upper_limit[bars$study == "Fullard"] == 2))
theme_uniform <- theme_classic(base_size=10, base_family="Arial") + theme(
  text=element_text(family="Arial"),
  plot.title=element_text(face="bold",size=12,margin=margin(b=3)),
  plot.subtitle=element_text(size=8,colour="#333333",margin=margin(b=7)),
  axis.title=element_text(size=10), axis.text=element_text(size=9,colour="black"),
  axis.line=element_line(linewidth=.45), axis.ticks=element_line(linewidth=.35),
  legend.title=element_text(size=8),legend.text=element_text(size=8),
  legend.position="right",legend.box.margin=margin(0,0,0,5),
  plot.margin=margin(8,8,8,8),plot.title.position="panel")
plots <- list(); cell_ids <- character(); marker_ids <- character()
coordinate_audit <- list(); scale_audit <- list(); label_audit <- list()
palette_audit <- list(); views <- list()
for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Preparing ",study)
  ct <- readRDS(file.path(task,"inputs/celltype",paste0(study,".rds")))
  mb <- readRDS(file.path(task,"inputs/marker",paste0(study,".rds")))
  stopifnot(length(ct$cells)==expected_n[[study]], !anyDuplicated(ct$cells),
    identical(rownames(ct$umap),ct$cells),identical(rownames(ct$metadata),ct$cells),
    identical(rownames(mb$umap),mb$cells),identical(colnames(mb$expression),mb$cells),
    identical(colnames(mb$counts),mb$cells),setequal(ct$cells,mb$cells))
  align <- match(mb$cells,ct$cells)
  ct_aligned <- ct$umap[align,,drop=FALSE]
  # Matrix subsetting drops optional scaled:center attributes. Compare every
  # coordinate exactly after validating dimensions and cell alignment.
  stopifnot(identical(dim(ct_aligned),dim(mb$umap)),
    identical(as.numeric(ct_aligned),as.numeric(mb$umap)))
  coords <- ct$umap
  xr <- range(coords[,1]); yr <- range(coords[,2])
  span <- max(diff(xr),diff(yr))*1.16
  xlim <- mean(xr)+c(-.5,.5)*span; ylim <- mean(yr)+c(-.5,.5)*span
  views[[study]] <- list(xlim=xlim,ylim=ylim)
  values <- as.character(ct$metadata[[annotation[[study]]]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  if (study=="Yao") {
    pc <- subset(yao_colours,plot_id=="Yao_current_UMAP_author_labels")
    anchors <- subset(yao_anchors,plot_id=="Yao_current_UMAP_author_labels")
    names(pc)[names(pc)=="author_Major"] <- "celltype"
    names(anchors)[names(anchors)=="type"] <- "celltype"
  } else {
    pc <- base_colours[base_colours$plot_id==study,,drop=FALSE]
    anchors <- base_anchors[base_anchors$plot_id==study,,drop=FALSE]
  }
  types <- pc$celltype; palette <- setNames(pc$colour,types)
  counts <- table(factor(values,levels=types))
  stopifnot(setequal(unique(values),types),identical(as.integer(counts),as.integer(pc$n_cells)),
    identical(as.integer(counts[anchors$celltype]),as.integer(anchors$n_cells)),
    all(anchors$anchor_cell %in% ct$cells))
  anchor_pos <- match(anchors$anchor_cell,ct$cells)
  stopifnot(max(abs(coords[anchor_pos,1]-anchors$x))<1e-12,
            max(abs(coords[anchor_pos,2]-anchors$y))<1e-12)
  data <- data.frame(x=coords[,1],y=coords[,2],celltype=factor(values,levels=types),row.names=ct$cells)
  subtitle <- if (study=="Yao") "73,721 cells; 16 author cell types; 174 unmatched" else
    paste0(comma(length(values))," cells; ",length(types)," cell types")
  title <- if(study=="Yao") "Yao | author cell types" else paste(study,"| cell types")
  p <- ggplot(data,aes(x=x,y=y,colour=celltype))+
    geom_scattermore(pointsize=2.2,pixels=c(1400,1400),show.legend=FALSE)+
    scale_colour_manual(values=palette,drop=FALSE)+
    geom_text_repel(data=anchors,aes(x=x,y=y,label=label_text),inherit.aes=FALSE,
      family="Arial",colour="black",size=if(length(types)>20) 2.35 else 2.75,
      lineheight=.95,seed=20261016,nudge_x=anchors$nudge_x,nudge_y=anchors$nudge_y,
      box.padding=unit(.35,"lines"),point.padding=unit(.10,"lines"),
      min.segment.length=0,segment.colour="#454545",segment.size=.20,
      force=2,force_pull=.3,max.overlaps=Inf,max.time=Inf,max.iter=20000,
      show.legend=FALSE)+
    coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
    labs(title=title,subtitle=subtitle,x="UMAP_1",y="UMAP_2")+
    theme_uniform+theme(legend.position="none")
  id <- paste0(study,"_celltype_UMAP")
  plots[[id]] <- p; cell_ids <- c(cell_ids,id)
  anchors$study <- study
  label_audit[[study]] <- anchors[,c("study","celltype","n_cells","x","y","anchor_cell","label_text","nudge_x","nudge_y")]
  palette_audit[[study]] <- data.frame(study=study,celltype=types,n_cells=as.integer(counts),colour=pc$colour)
  coordinate_audit[[study]] <- data.frame(study=study,n_cells=length(values),annotation_column=annotation[[study]],
    n_celltype_labels=length(types),celltype_and_marker_barcode_sets_identical=TRUE,
    celltype_and_marker_coordinates_identical_by_barcode=TRUE,
    coordinates_sha256=digest(coords,algo="sha256"),
    celltype_values_sha256=digest(values,algo="sha256"),
    marker_expression_sha256=digest(mb$expression,algo="sha256"),
    original_x_min=xr[1],original_x_max=xr[2],original_y_min=yr[1],original_y_max=yr[2],
    view_x_min=xlim[1],view_x_max=xlim[2],view_y_min=ylim[1],view_y_max=ylim[2],
    coordinate_ratio=1,umap_recomputed=FALSE,annotation_modified=FALSE)
  thin <- CreateSeuratObject(counts=mb$counts,min.cells=0,min.features=0,project=study)
  thin <- SetAssayData(thin,assay="RNA",layer="data",new.data=mb$expression)
  mc <- mb$umap; colnames(mc) <- c("UMAP_1","UMAP_2")
  thin[["umap"]] <- CreateDimReducObject(embeddings=mc,key="UMAP_",assay="RNA")
  stopifnot(identical(dim(Embeddings(thin,"umap")),dim(mb$umap)),
    identical(as.numeric(Embeddings(thin,"umap")),as.numeric(mb$umap)),
    identical(as.matrix(GetAssayData(thin,assay="RNA",layer="data")),as.matrix(mb$expression)))
  for (gene in markers) {
    bar <- bars[bars$study==study & bars$marker==gene,,drop=FALSE]
    stopifnot(nrow(bar)==1L)
    feature <- mb$features[[gene]]; upper <- bar$upper_limit
    tick_values <- as.numeric(strsplit(bar$tick_values,";",fixed=TRUE)[[1]])
    tick_labels <- strsplit(bar$tick_labels,";",fixed=TRUE)[[1]]
    vals <- as.numeric(mb$expression[feature,])
    p <- FeaturePlot(thin,features=feature,reduction="umap",dims=c(1,2),
      cols=c("#D3D3D3","#E31A1C"),order=TRUE,pt.size=1.3,
      min.cutoff=0,max.cutoff=upper,slot="data",raster=TRUE,raster.dpi=c(1200,1200),combine=FALSE)[[1]]+
      scale_colour_gradient(low=bar$colour_low,high=bar$colour_high,
        limits=c(0,upper),oob=scales::squish,breaks=tick_values,labels=tick_labels,
        guide=guide_colourbar(title="Log-normalized\nRNA",barheight=unit(2.6,"cm"),barwidth=unit(.30,"cm")))+
      coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
      labs(title=paste(study,"|",gene),subtitle=paste0(comma(length(vals))," cells"),x="UMAP_1",y="UMAP_2")+
      theme_uniform
    # FeaturePlot must retain the cached values after applying only the chosen cap.
    j <- match(rownames(p$data),mb$cells)
    stopifnot(!anyNA(j),nrow(p$data)==length(mb$cells),!anyDuplicated(rownames(p$data)),
      setequal(rownames(p$data),mb$cells),
      identical(as.numeric(p$data[[feature]]),pmin(pmax(vals[j],0),upper)))
    id <- paste0(study,"_",gene,"_FeaturePlot")
    plots[[id]] <- p; marker_ids <- c(marker_ids,id)
    scale_audit[[id]] <- data.frame(study=study,marker=gene,lower_limit=0,upper_limit=upper,
      visible_breaks=bar$tick_values,visible_tick_labels=bar$tick_labels,
      colour_low=bar$colour_low,colour_high=bar$colour_high,n_cells=length(vals),
      n_positive=sum(vals>0),n_saturated=sum(vals>upper),selected_pdf_sha256=bar$sha256,
      expression_identical_before_display_cap=TRUE,sorted_ascending=TRUE,
      panel_width_in=4.5,panel_height_in=4.5,coordinate_ratio=1)
  }
  rm(ct,mb,thin); gc(verbose=FALSE)
}

# Measure on the actual Cairo/Arial PDF device. Fix the panel row and column,
# then equalize every surrounding row/column so no legend can shrink a panel.
measurement <- file.path(task,"runtime/layout_measurement.pdf")
cairo_pdf(measurement,width=9,height=9,family="Arial")
grobs <- lapply(plots,ggplotGrob)
panel_mm <- 114.3 # 4.5 inches, identical in all 24 panels
shape <- vapply(grobs,function(g) paste(length(g$widths),length(g$heights)),character(1))
stopifnot(length(unique(shape))==1L)
panel_locations <- vapply(grobs,function(g) {
  loc <- g$layout[g$layout$name=="panel",]
  paste(loc$l,loc$r,loc$t,loc$b,collapse=";")
},character(1))
stopifnot(length(unique(panel_locations))==1L)
dimensions <- list()
for (id in names(grobs)) {
  g <- grobs[[id]]; loc <- g$layout[g$layout$name=="panel",]
  stopifnot(nrow(loc)==1L,loc$l==loc$r,loc$t==loc$b)
  g$widths[loc$l] <- unit(panel_mm,"mm");g$heights[loc$t] <- unit(panel_mm,"mm")
  grobs[[id]] <- g
}
all_widths <- lapply(grobs,function(g) convertWidth(g$widths,"in",valueOnly=TRUE))
all_heights <- lapply(grobs,function(g) convertHeight(g$heights,"in",valueOnly=TRUE))
common_widths <- do.call(pmax,all_widths)
common_heights <- do.call(pmax,all_heights)
for (id in names(grobs)) {
  g <- grobs[[id]];g$widths <- unit(common_widths,"in");g$heights <- unit(common_heights,"in")
  grobs[[id]] <- g
  loc <- g$layout[g$layout$name=="panel",]
  pw <- convertWidth(g$widths[loc$l],"mm",valueOnly=TRUE)
  ph <- convertHeight(g$heights[loc$t],"mm",valueOnly=TRUE)
  stopifnot(abs(pw-panel_mm)<1e-9,abs(ph-panel_mm)<1e-9)
  dimensions[[id]] <- data.frame(plot_id=id,panel_width_mm=pw,panel_height_mm=ph,
    outer_width_in=sum(common_widths),outer_height_in=sum(common_heights),font="Arial")
}
dev.off();unlink(measurement)
slot_width <- sum(common_widths);slot_height <- sum(common_heights)
draw_page <- function(ids,ncol,nrow) {
  stopifnot(length(ids)==ncol*nrow)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol,
    widths=unit(rep(slot_width,ncol),"in"),heights=unit(rep(slot_height,nrow),"in"))))
  for (i in seq_along(ids)) {
    pushViewport(viewport(layout.pos.row=(i-1L)%/%ncol+1L,layout.pos.col=(i-1L)%%ncol+1L))
    grid.draw(grobs[[ids[[i]]]])
    upViewport()
  }
  upViewport()
}
outputs <- list()
save_page <- function(ids,name,ncol=1L,nrow=1L) {
  w <- slot_width*ncol;h <- slot_height*nrow
  message(format(Sys.time(),"%F %T")," Saving ",name," (",round(w,2)," x ",round(h,2)," in)")
  cairo_pdf(file.path(task,"output/pdf",paste0(name,".pdf")),width=w,height=h,family="Arial")
  draw_page(ids,ncol,nrow);dev.off()
  ragg::agg_png(file.path(task,"output/png",paste0(name,".png")),width=w,height=h,units="in",res=240,background="white")
  draw_page(ids,ncol,nrow);dev.off()
  outputs[[name]] <<- data.frame(name=name,n_panels=length(ids),n_columns=ncol,n_rows=nrow,
    width_in=w,height_in=h,pdf=paste0("output/pdf/",name,".pdf"),png=paste0("output/png/",name,".png"))
}
for (id in c(cell_ids,marker_ids)) save_page(id,id)
save_page(cell_ids,"Six_study_celltype_UMAPs_Arial_square",3L,2L)
save_page(marker_ids[1:9],"Marker_UMAPs_3x3_page1_Chia_Fullard_Grubman",3L,3L)
save_page(marker_ids[10:18],"Marker_UMAPs_3x3_page2_Jakel_Su_Yao",3L,3L)
write.csv(do.call(rbind,coordinate_audit),file.path(task,"audit/coordinate_and_annotation_preservation.csv"),row.names=FALSE)
write.csv(do.call(rbind,scale_audit),file.path(task,"audit/marker_scale_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,label_audit),file.path(task,"audit/celltype_labels_and_anchors.csv"),row.names=FALSE)
write.csv(do.call(rbind,palette_audit),file.path(task,"audit/celltype_counts_and_colours.csv"),row.names=FALSE)
write.csv(do.call(rbind,dimensions),file.path(task,"audit/physical_panel_dimensions.csv"),row.names=FALSE)
write.csv(do.call(rbind,outputs),file.path(task,"audit/output_manifest.csv"),row.names=FALSE)
saveRDS(views,file.path(task,"audit/shared_square_view_limits.rds"))
writeLines(capture.output(sessionInfo()),file.path(task,"audit/plotting_sessionInfo.txt"))
message("UNIFORM_UMAP_PLOTTING_COMPLETE: 27 PDF and PNG files; 24 equal 114.3-mm square panels")

# S2B  Fullard P2RY12 / CX3CR1 / CSF1R UMAPs
# Source: source/project/Reference/combine/results/source_celltype_umap_20260920/revisions/20261022_Arial_square_celltype_and_marker_UMAPs/code/plot_uniform_umaps.R
#!/usr/bin/env Rscript
# Display-only revision: exact selected inputs, selected marker bars, Arial,
# one shared square view per study and an identical physical data panel.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(scattermore)
  library(scales); library(grid); library(gtable); library(digest)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task <- normalizePath(args[[1]], mustWork = TRUE)
for (d in c("output/pdf", "output/png", "audit", "runtime"))
  dir.create(file.path(task, d), recursive = TRUE, showWarnings = FALSE)
studies <- c("Chia", "Fullard", "Grubman", "Jakel", "Su", "Yao")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
annotation <- c(Chia="cell.type", Fullard="celltype", Grubman="cellType",
                Jakel="CellType", Su="CellType", Yao="AuthorMajor")
expected_n <- c(Chia=87470L, Fullard=68557L, Grubman=13214L, Jakel=17799L,
                Su=224464L, Yao=73721L)
font <- systemfonts::match_fonts("Arial")
stopifnot(grepl("Arial", basename(font$path), fixed = TRUE))
writeLines(font$path, file.path(task, "audit/Arial_font_path.txt"))
base_colours <- read.csv(file.path(task,"inputs/source_tables/base_celltype_colours.csv"), stringsAsFactors=FALSE)
base_anchors <- read.csv(file.path(task,"inputs/source_tables/base_celltype_anchors.csv"), stringsAsFactors=FALSE)
yao_colours <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_colours.csv"), stringsAsFactors=FALSE)
yao_anchors <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_anchors.csv"), stringsAsFactors=FALSE)
bars <- read.csv(file.path(task,"audit/selected_marker_bar_provenance.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(bars)==18L, all(toupper(as.character(bars$exact_match))=="TRUE"),
          all(bars$upper_limit[bars$study %in% c("Grubman","Jakel")] == 1),
          all(bars$upper_limit[bars$study == "Fullard"] == 2))
theme_uniform <- theme_classic(base_size=10, base_family="Arial") + theme(
  text=element_text(family="Arial"),
  plot.title=element_text(face="bold",size=12,margin=margin(b=3)),
  plot.subtitle=element_text(size=8,colour="#333333",margin=margin(b=7)),
  axis.title=element_text(size=10), axis.text=element_text(size=9,colour="black"),
  axis.line=element_line(linewidth=.45), axis.ticks=element_line(linewidth=.35),
  legend.title=element_text(size=8),legend.text=element_text(size=8),
  legend.position="right",legend.box.margin=margin(0,0,0,5),
  plot.margin=margin(8,8,8,8),plot.title.position="panel")
plots <- list(); cell_ids <- character(); marker_ids <- character()
coordinate_audit <- list(); scale_audit <- list(); label_audit <- list()
palette_audit <- list(); views <- list()
for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Preparing ",study)
  ct <- readRDS(file.path(task,"inputs/celltype",paste0(study,".rds")))
  mb <- readRDS(file.path(task,"inputs/marker",paste0(study,".rds")))
  stopifnot(length(ct$cells)==expected_n[[study]], !anyDuplicated(ct$cells),
    identical(rownames(ct$umap),ct$cells),identical(rownames(ct$metadata),ct$cells),
    identical(rownames(mb$umap),mb$cells),identical(colnames(mb$expression),mb$cells),
    identical(colnames(mb$counts),mb$cells),setequal(ct$cells,mb$cells))
  align <- match(mb$cells,ct$cells)
  ct_aligned <- ct$umap[align,,drop=FALSE]
  # Matrix subsetting drops optional scaled:center attributes. Compare every
  # coordinate exactly after validating dimensions and cell alignment.
  stopifnot(identical(dim(ct_aligned),dim(mb$umap)),
    identical(as.numeric(ct_aligned),as.numeric(mb$umap)))
  coords <- ct$umap
  xr <- range(coords[,1]); yr <- range(coords[,2])
  span <- max(diff(xr),diff(yr))*1.16
  xlim <- mean(xr)+c(-.5,.5)*span; ylim <- mean(yr)+c(-.5,.5)*span
  views[[study]] <- list(xlim=xlim,ylim=ylim)
  values <- as.character(ct$metadata[[annotation[[study]]]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  if (study=="Yao") {
    pc <- subset(yao_colours,plot_id=="Yao_current_UMAP_author_labels")
    anchors <- subset(yao_anchors,plot_id=="Yao_current_UMAP_author_labels")
    names(pc)[names(pc)=="author_Major"] <- "celltype"
    names(anchors)[names(anchors)=="type"] <- "celltype"
  } else {
    pc <- base_colours[base_colours$plot_id==study,,drop=FALSE]
    anchors <- base_anchors[base_anchors$plot_id==study,,drop=FALSE]
  }
  types <- pc$celltype; palette <- setNames(pc$colour,types)
  counts <- table(factor(values,levels=types))
  stopifnot(setequal(unique(values),types),identical(as.integer(counts),as.integer(pc$n_cells)),
    identical(as.integer(counts[anchors$celltype]),as.integer(anchors$n_cells)),
    all(anchors$anchor_cell %in% ct$cells))
  anchor_pos <- match(anchors$anchor_cell,ct$cells)
  stopifnot(max(abs(coords[anchor_pos,1]-anchors$x))<1e-12,
            max(abs(coords[anchor_pos,2]-anchors$y))<1e-12)
  data <- data.frame(x=coords[,1],y=coords[,2],celltype=factor(values,levels=types),row.names=ct$cells)
  subtitle <- if (study=="Yao") "73,721 cells; 16 author cell types; 174 unmatched" else
    paste0(comma(length(values))," cells; ",length(types)," cell types")
  title <- if(study=="Yao") "Yao | author cell types" else paste(study,"| cell types")
  p <- ggplot(data,aes(x=x,y=y,colour=celltype))+
    geom_scattermore(pointsize=2.2,pixels=c(1400,1400),show.legend=FALSE)+
    scale_colour_manual(values=palette,drop=FALSE)+
    geom_text_repel(data=anchors,aes(x=x,y=y,label=label_text),inherit.aes=FALSE,
      family="Arial",colour="black",size=if(length(types)>20) 2.35 else 2.75,
      lineheight=.95,seed=20261016,nudge_x=anchors$nudge_x,nudge_y=anchors$nudge_y,
      box.padding=unit(.35,"lines"),point.padding=unit(.10,"lines"),
      min.segment.length=0,segment.colour="#454545",segment.size=.20,
      force=2,force_pull=.3,max.overlaps=Inf,max.time=Inf,max.iter=20000,
      show.legend=FALSE)+
    coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
    labs(title=title,subtitle=subtitle,x="UMAP_1",y="UMAP_2")+
    theme_uniform+theme(legend.position="none")
  id <- paste0(study,"_celltype_UMAP")
  plots[[id]] <- p; cell_ids <- c(cell_ids,id)
  anchors$study <- study
  label_audit[[study]] <- anchors[,c("study","celltype","n_cells","x","y","anchor_cell","label_text","nudge_x","nudge_y")]
  palette_audit[[study]] <- data.frame(study=study,celltype=types,n_cells=as.integer(counts),colour=pc$colour)
  coordinate_audit[[study]] <- data.frame(study=study,n_cells=length(values),annotation_column=annotation[[study]],
    n_celltype_labels=length(types),celltype_and_marker_barcode_sets_identical=TRUE,
    celltype_and_marker_coordinates_identical_by_barcode=TRUE,
    coordinates_sha256=digest(coords,algo="sha256"),
    celltype_values_sha256=digest(values,algo="sha256"),
    marker_expression_sha256=digest(mb$expression,algo="sha256"),
    original_x_min=xr[1],original_x_max=xr[2],original_y_min=yr[1],original_y_max=yr[2],
    view_x_min=xlim[1],view_x_max=xlim[2],view_y_min=ylim[1],view_y_max=ylim[2],
    coordinate_ratio=1,umap_recomputed=FALSE,annotation_modified=FALSE)
  thin <- CreateSeuratObject(counts=mb$counts,min.cells=0,min.features=0,project=study)
  thin <- SetAssayData(thin,assay="RNA",layer="data",new.data=mb$expression)
  mc <- mb$umap; colnames(mc) <- c("UMAP_1","UMAP_2")
  thin[["umap"]] <- CreateDimReducObject(embeddings=mc,key="UMAP_",assay="RNA")
  stopifnot(identical(dim(Embeddings(thin,"umap")),dim(mb$umap)),
    identical(as.numeric(Embeddings(thin,"umap")),as.numeric(mb$umap)),
    identical(as.matrix(GetAssayData(thin,assay="RNA",layer="data")),as.matrix(mb$expression)))
  for (gene in markers) {
    bar <- bars[bars$study==study & bars$marker==gene,,drop=FALSE]
    stopifnot(nrow(bar)==1L)
    feature <- mb$features[[gene]]; upper <- bar$upper_limit
    tick_values <- as.numeric(strsplit(bar$tick_values,";",fixed=TRUE)[[1]])
    tick_labels <- strsplit(bar$tick_labels,";",fixed=TRUE)[[1]]
    vals <- as.numeric(mb$expression[feature,])
    p <- FeaturePlot(thin,features=feature,reduction="umap",dims=c(1,2),
      cols=c("#D3D3D3","#E31A1C"),order=TRUE,pt.size=1.3,
      min.cutoff=0,max.cutoff=upper,slot="data",raster=TRUE,raster.dpi=c(1200,1200),combine=FALSE)[[1]]+
      scale_colour_gradient(low=bar$colour_low,high=bar$colour_high,
        limits=c(0,upper),oob=scales::squish,breaks=tick_values,labels=tick_labels,
        guide=guide_colourbar(title="Log-normalized\nRNA",barheight=unit(2.6,"cm"),barwidth=unit(.30,"cm")))+
      coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
      labs(title=paste(study,"|",gene),subtitle=paste0(comma(length(vals))," cells"),x="UMAP_1",y="UMAP_2")+
      theme_uniform
    # FeaturePlot must retain the cached values after applying only the chosen cap.
    j <- match(rownames(p$data),mb$cells)
    stopifnot(!anyNA(j),nrow(p$data)==length(mb$cells),!anyDuplicated(rownames(p$data)),
      setequal(rownames(p$data),mb$cells),
      identical(as.numeric(p$data[[feature]]),pmin(pmax(vals[j],0),upper)))
    id <- paste0(study,"_",gene,"_FeaturePlot")
    plots[[id]] <- p; marker_ids <- c(marker_ids,id)
    scale_audit[[id]] <- data.frame(study=study,marker=gene,lower_limit=0,upper_limit=upper,
      visible_breaks=bar$tick_values,visible_tick_labels=bar$tick_labels,
      colour_low=bar$colour_low,colour_high=bar$colour_high,n_cells=length(vals),
      n_positive=sum(vals>0),n_saturated=sum(vals>upper),selected_pdf_sha256=bar$sha256,
      expression_identical_before_display_cap=TRUE,sorted_ascending=TRUE,
      panel_width_in=4.5,panel_height_in=4.5,coordinate_ratio=1)
  }
  rm(ct,mb,thin); gc(verbose=FALSE)
}

# Measure on the actual Cairo/Arial PDF device. Fix the panel row and column,
# then equalize every surrounding row/column so no legend can shrink a panel.
measurement <- file.path(task,"runtime/layout_measurement.pdf")
cairo_pdf(measurement,width=9,height=9,family="Arial")
grobs <- lapply(plots,ggplotGrob)
panel_mm <- 114.3 # 4.5 inches, identical in all 24 panels
shape <- vapply(grobs,function(g) paste(length(g$widths),length(g$heights)),character(1))
stopifnot(length(unique(shape))==1L)
panel_locations <- vapply(grobs,function(g) {
  loc <- g$layout[g$layout$name=="panel",]
  paste(loc$l,loc$r,loc$t,loc$b,collapse=";")
},character(1))
stopifnot(length(unique(panel_locations))==1L)
dimensions <- list()
for (id in names(grobs)) {
  g <- grobs[[id]]; loc <- g$layout[g$layout$name=="panel",]
  stopifnot(nrow(loc)==1L,loc$l==loc$r,loc$t==loc$b)
  g$widths[loc$l] <- unit(panel_mm,"mm");g$heights[loc$t] <- unit(panel_mm,"mm")
  grobs[[id]] <- g
}
all_widths <- lapply(grobs,function(g) convertWidth(g$widths,"in",valueOnly=TRUE))
all_heights <- lapply(grobs,function(g) convertHeight(g$heights,"in",valueOnly=TRUE))
common_widths <- do.call(pmax,all_widths)
common_heights <- do.call(pmax,all_heights)
for (id in names(grobs)) {
  g <- grobs[[id]];g$widths <- unit(common_widths,"in");g$heights <- unit(common_heights,"in")
  grobs[[id]] <- g
  loc <- g$layout[g$layout$name=="panel",]
  pw <- convertWidth(g$widths[loc$l],"mm",valueOnly=TRUE)
  ph <- convertHeight(g$heights[loc$t],"mm",valueOnly=TRUE)
  stopifnot(abs(pw-panel_mm)<1e-9,abs(ph-panel_mm)<1e-9)
  dimensions[[id]] <- data.frame(plot_id=id,panel_width_mm=pw,panel_height_mm=ph,
    outer_width_in=sum(common_widths),outer_height_in=sum(common_heights),font="Arial")
}
dev.off();unlink(measurement)
slot_width <- sum(common_widths);slot_height <- sum(common_heights)
draw_page <- function(ids,ncol,nrow) {
  stopifnot(length(ids)==ncol*nrow)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol,
    widths=unit(rep(slot_width,ncol),"in"),heights=unit(rep(slot_height,nrow),"in"))))
  for (i in seq_along(ids)) {
    pushViewport(viewport(layout.pos.row=(i-1L)%/%ncol+1L,layout.pos.col=(i-1L)%%ncol+1L))
    grid.draw(grobs[[ids[[i]]]])
    upViewport()
  }
  upViewport()
}
outputs <- list()
save_page <- function(ids,name,ncol=1L,nrow=1L) {
  w <- slot_width*ncol;h <- slot_height*nrow
  message(format(Sys.time(),"%F %T")," Saving ",name," (",round(w,2)," x ",round(h,2)," in)")
  cairo_pdf(file.path(task,"output/pdf",paste0(name,".pdf")),width=w,height=h,family="Arial")
  draw_page(ids,ncol,nrow);dev.off()
  ragg::agg_png(file.path(task,"output/png",paste0(name,".png")),width=w,height=h,units="in",res=240,background="white")
  draw_page(ids,ncol,nrow);dev.off()
  outputs[[name]] <<- data.frame(name=name,n_panels=length(ids),n_columns=ncol,n_rows=nrow,
    width_in=w,height_in=h,pdf=paste0("output/pdf/",name,".pdf"),png=paste0("output/png/",name,".png"))
}
for (id in c(cell_ids,marker_ids)) save_page(id,id)
save_page(cell_ids,"Six_study_celltype_UMAPs_Arial_square",3L,2L)
save_page(marker_ids[1:9],"Marker_UMAPs_3x3_page1_Chia_Fullard_Grubman",3L,3L)
save_page(marker_ids[10:18],"Marker_UMAPs_3x3_page2_Jakel_Su_Yao",3L,3L)
write.csv(do.call(rbind,coordinate_audit),file.path(task,"audit/coordinate_and_annotation_preservation.csv"),row.names=FALSE)
write.csv(do.call(rbind,scale_audit),file.path(task,"audit/marker_scale_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,label_audit),file.path(task,"audit/celltype_labels_and_anchors.csv"),row.names=FALSE)
write.csv(do.call(rbind,palette_audit),file.path(task,"audit/celltype_counts_and_colours.csv"),row.names=FALSE)
write.csv(do.call(rbind,dimensions),file.path(task,"audit/physical_panel_dimensions.csv"),row.names=FALSE)
write.csv(do.call(rbind,outputs),file.path(task,"audit/output_manifest.csv"),row.names=FALSE)
saveRDS(views,file.path(task,"audit/shared_square_view_limits.rds"))
writeLines(capture.output(sessionInfo()),file.path(task,"audit/plotting_sessionInfo.txt"))
message("UNIFORM_UMAP_PLOTTING_COMPLETE: 27 PDF and PNG files; 24 equal 114.3-mm square panels")

# S2C  Grubman P2RY12 / CX3CR1 / CSF1R UMAPs
# Source: source/project/Reference/combine/results/source_celltype_umap_20260920/revisions/20261022_Arial_square_celltype_and_marker_UMAPs/code/plot_uniform_umaps.R
#!/usr/bin/env Rscript
# Display-only revision: exact selected inputs, selected marker bars, Arial,
# one shared square view per study and an identical physical data panel.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(scattermore)
  library(scales); library(grid); library(gtable); library(digest)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task <- normalizePath(args[[1]], mustWork = TRUE)
for (d in c("output/pdf", "output/png", "audit", "runtime"))
  dir.create(file.path(task, d), recursive = TRUE, showWarnings = FALSE)
studies <- c("Chia", "Fullard", "Grubman", "Jakel", "Su", "Yao")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
annotation <- c(Chia="cell.type", Fullard="celltype", Grubman="cellType",
                Jakel="CellType", Su="CellType", Yao="AuthorMajor")
expected_n <- c(Chia=87470L, Fullard=68557L, Grubman=13214L, Jakel=17799L,
                Su=224464L, Yao=73721L)
font <- systemfonts::match_fonts("Arial")
stopifnot(grepl("Arial", basename(font$path), fixed = TRUE))
writeLines(font$path, file.path(task, "audit/Arial_font_path.txt"))
base_colours <- read.csv(file.path(task,"inputs/source_tables/base_celltype_colours.csv"), stringsAsFactors=FALSE)
base_anchors <- read.csv(file.path(task,"inputs/source_tables/base_celltype_anchors.csv"), stringsAsFactors=FALSE)
yao_colours <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_colours.csv"), stringsAsFactors=FALSE)
yao_anchors <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_anchors.csv"), stringsAsFactors=FALSE)
bars <- read.csv(file.path(task,"audit/selected_marker_bar_provenance.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(bars)==18L, all(toupper(as.character(bars$exact_match))=="TRUE"),
          all(bars$upper_limit[bars$study %in% c("Grubman","Jakel")] == 1),
          all(bars$upper_limit[bars$study == "Fullard"] == 2))
theme_uniform <- theme_classic(base_size=10, base_family="Arial") + theme(
  text=element_text(family="Arial"),
  plot.title=element_text(face="bold",size=12,margin=margin(b=3)),
  plot.subtitle=element_text(size=8,colour="#333333",margin=margin(b=7)),
  axis.title=element_text(size=10), axis.text=element_text(size=9,colour="black"),
  axis.line=element_line(linewidth=.45), axis.ticks=element_line(linewidth=.35),
  legend.title=element_text(size=8),legend.text=element_text(size=8),
  legend.position="right",legend.box.margin=margin(0,0,0,5),
  plot.margin=margin(8,8,8,8),plot.title.position="panel")
plots <- list(); cell_ids <- character(); marker_ids <- character()
coordinate_audit <- list(); scale_audit <- list(); label_audit <- list()
palette_audit <- list(); views <- list()
for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Preparing ",study)
  ct <- readRDS(file.path(task,"inputs/celltype",paste0(study,".rds")))
  mb <- readRDS(file.path(task,"inputs/marker",paste0(study,".rds")))
  stopifnot(length(ct$cells)==expected_n[[study]], !anyDuplicated(ct$cells),
    identical(rownames(ct$umap),ct$cells),identical(rownames(ct$metadata),ct$cells),
    identical(rownames(mb$umap),mb$cells),identical(colnames(mb$expression),mb$cells),
    identical(colnames(mb$counts),mb$cells),setequal(ct$cells,mb$cells))
  align <- match(mb$cells,ct$cells)
  ct_aligned <- ct$umap[align,,drop=FALSE]
  # Matrix subsetting drops optional scaled:center attributes. Compare every
  # coordinate exactly after validating dimensions and cell alignment.
  stopifnot(identical(dim(ct_aligned),dim(mb$umap)),
    identical(as.numeric(ct_aligned),as.numeric(mb$umap)))
  coords <- ct$umap
  xr <- range(coords[,1]); yr <- range(coords[,2])
  span <- max(diff(xr),diff(yr))*1.16
  xlim <- mean(xr)+c(-.5,.5)*span; ylim <- mean(yr)+c(-.5,.5)*span
  views[[study]] <- list(xlim=xlim,ylim=ylim)
  values <- as.character(ct$metadata[[annotation[[study]]]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  if (study=="Yao") {
    pc <- subset(yao_colours,plot_id=="Yao_current_UMAP_author_labels")
    anchors <- subset(yao_anchors,plot_id=="Yao_current_UMAP_author_labels")
    names(pc)[names(pc)=="author_Major"] <- "celltype"
    names(anchors)[names(anchors)=="type"] <- "celltype"
  } else {
    pc <- base_colours[base_colours$plot_id==study,,drop=FALSE]
    anchors <- base_anchors[base_anchors$plot_id==study,,drop=FALSE]
  }
  types <- pc$celltype; palette <- setNames(pc$colour,types)
  counts <- table(factor(values,levels=types))
  stopifnot(setequal(unique(values),types),identical(as.integer(counts),as.integer(pc$n_cells)),
    identical(as.integer(counts[anchors$celltype]),as.integer(anchors$n_cells)),
    all(anchors$anchor_cell %in% ct$cells))
  anchor_pos <- match(anchors$anchor_cell,ct$cells)
  stopifnot(max(abs(coords[anchor_pos,1]-anchors$x))<1e-12,
            max(abs(coords[anchor_pos,2]-anchors$y))<1e-12)
  data <- data.frame(x=coords[,1],y=coords[,2],celltype=factor(values,levels=types),row.names=ct$cells)
  subtitle <- if (study=="Yao") "73,721 cells; 16 author cell types; 174 unmatched" else
    paste0(comma(length(values))," cells; ",length(types)," cell types")
  title <- if(study=="Yao") "Yao | author cell types" else paste(study,"| cell types")
  p <- ggplot(data,aes(x=x,y=y,colour=celltype))+
    geom_scattermore(pointsize=2.2,pixels=c(1400,1400),show.legend=FALSE)+
    scale_colour_manual(values=palette,drop=FALSE)+
    geom_text_repel(data=anchors,aes(x=x,y=y,label=label_text),inherit.aes=FALSE,
      family="Arial",colour="black",size=if(length(types)>20) 2.35 else 2.75,
      lineheight=.95,seed=20261016,nudge_x=anchors$nudge_x,nudge_y=anchors$nudge_y,
      box.padding=unit(.35,"lines"),point.padding=unit(.10,"lines"),
      min.segment.length=0,segment.colour="#454545",segment.size=.20,
      force=2,force_pull=.3,max.overlaps=Inf,max.time=Inf,max.iter=20000,
      show.legend=FALSE)+
    coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
    labs(title=title,subtitle=subtitle,x="UMAP_1",y="UMAP_2")+
    theme_uniform+theme(legend.position="none")
  id <- paste0(study,"_celltype_UMAP")
  plots[[id]] <- p; cell_ids <- c(cell_ids,id)
  anchors$study <- study
  label_audit[[study]] <- anchors[,c("study","celltype","n_cells","x","y","anchor_cell","label_text","nudge_x","nudge_y")]
  palette_audit[[study]] <- data.frame(study=study,celltype=types,n_cells=as.integer(counts),colour=pc$colour)
  coordinate_audit[[study]] <- data.frame(study=study,n_cells=length(values),annotation_column=annotation[[study]],
    n_celltype_labels=length(types),celltype_and_marker_barcode_sets_identical=TRUE,
    celltype_and_marker_coordinates_identical_by_barcode=TRUE,
    coordinates_sha256=digest(coords,algo="sha256"),
    celltype_values_sha256=digest(values,algo="sha256"),
    marker_expression_sha256=digest(mb$expression,algo="sha256"),
    original_x_min=xr[1],original_x_max=xr[2],original_y_min=yr[1],original_y_max=yr[2],
    view_x_min=xlim[1],view_x_max=xlim[2],view_y_min=ylim[1],view_y_max=ylim[2],
    coordinate_ratio=1,umap_recomputed=FALSE,annotation_modified=FALSE)
  thin <- CreateSeuratObject(counts=mb$counts,min.cells=0,min.features=0,project=study)
  thin <- SetAssayData(thin,assay="RNA",layer="data",new.data=mb$expression)
  mc <- mb$umap; colnames(mc) <- c("UMAP_1","UMAP_2")
  thin[["umap"]] <- CreateDimReducObject(embeddings=mc,key="UMAP_",assay="RNA")
  stopifnot(identical(dim(Embeddings(thin,"umap")),dim(mb$umap)),
    identical(as.numeric(Embeddings(thin,"umap")),as.numeric(mb$umap)),
    identical(as.matrix(GetAssayData(thin,assay="RNA",layer="data")),as.matrix(mb$expression)))
  for (gene in markers) {
    bar <- bars[bars$study==study & bars$marker==gene,,drop=FALSE]
    stopifnot(nrow(bar)==1L)
    feature <- mb$features[[gene]]; upper <- bar$upper_limit
    tick_values <- as.numeric(strsplit(bar$tick_values,";",fixed=TRUE)[[1]])
    tick_labels <- strsplit(bar$tick_labels,";",fixed=TRUE)[[1]]
    vals <- as.numeric(mb$expression[feature,])
    p <- FeaturePlot(thin,features=feature,reduction="umap",dims=c(1,2),
      cols=c("#D3D3D3","#E31A1C"),order=TRUE,pt.size=1.3,
      min.cutoff=0,max.cutoff=upper,slot="data",raster=TRUE,raster.dpi=c(1200,1200),combine=FALSE)[[1]]+
      scale_colour_gradient(low=bar$colour_low,high=bar$colour_high,
        limits=c(0,upper),oob=scales::squish,breaks=tick_values,labels=tick_labels,
        guide=guide_colourbar(title="Log-normalized\nRNA",barheight=unit(2.6,"cm"),barwidth=unit(.30,"cm")))+
      coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
      labs(title=paste(study,"|",gene),subtitle=paste0(comma(length(vals))," cells"),x="UMAP_1",y="UMAP_2")+
      theme_uniform
    # FeaturePlot must retain the cached values after applying only the chosen cap.
    j <- match(rownames(p$data),mb$cells)
    stopifnot(!anyNA(j),nrow(p$data)==length(mb$cells),!anyDuplicated(rownames(p$data)),
      setequal(rownames(p$data),mb$cells),
      identical(as.numeric(p$data[[feature]]),pmin(pmax(vals[j],0),upper)))
    id <- paste0(study,"_",gene,"_FeaturePlot")
    plots[[id]] <- p; marker_ids <- c(marker_ids,id)
    scale_audit[[id]] <- data.frame(study=study,marker=gene,lower_limit=0,upper_limit=upper,
      visible_breaks=bar$tick_values,visible_tick_labels=bar$tick_labels,
      colour_low=bar$colour_low,colour_high=bar$colour_high,n_cells=length(vals),
      n_positive=sum(vals>0),n_saturated=sum(vals>upper),selected_pdf_sha256=bar$sha256,
      expression_identical_before_display_cap=TRUE,sorted_ascending=TRUE,
      panel_width_in=4.5,panel_height_in=4.5,coordinate_ratio=1)
  }
  rm(ct,mb,thin); gc(verbose=FALSE)
}

# Measure on the actual Cairo/Arial PDF device. Fix the panel row and column,
# then equalize every surrounding row/column so no legend can shrink a panel.
measurement <- file.path(task,"runtime/layout_measurement.pdf")
cairo_pdf(measurement,width=9,height=9,family="Arial")
grobs <- lapply(plots,ggplotGrob)
panel_mm <- 114.3 # 4.5 inches, identical in all 24 panels
shape <- vapply(grobs,function(g) paste(length(g$widths),length(g$heights)),character(1))
stopifnot(length(unique(shape))==1L)
panel_locations <- vapply(grobs,function(g) {
  loc <- g$layout[g$layout$name=="panel",]
  paste(loc$l,loc$r,loc$t,loc$b,collapse=";")
},character(1))
stopifnot(length(unique(panel_locations))==1L)
dimensions <- list()
for (id in names(grobs)) {
  g <- grobs[[id]]; loc <- g$layout[g$layout$name=="panel",]
  stopifnot(nrow(loc)==1L,loc$l==loc$r,loc$t==loc$b)
  g$widths[loc$l] <- unit(panel_mm,"mm");g$heights[loc$t] <- unit(panel_mm,"mm")
  grobs[[id]] <- g
}
all_widths <- lapply(grobs,function(g) convertWidth(g$widths,"in",valueOnly=TRUE))
all_heights <- lapply(grobs,function(g) convertHeight(g$heights,"in",valueOnly=TRUE))
common_widths <- do.call(pmax,all_widths)
common_heights <- do.call(pmax,all_heights)
for (id in names(grobs)) {
  g <- grobs[[id]];g$widths <- unit(common_widths,"in");g$heights <- unit(common_heights,"in")
  grobs[[id]] <- g
  loc <- g$layout[g$layout$name=="panel",]
  pw <- convertWidth(g$widths[loc$l],"mm",valueOnly=TRUE)
  ph <- convertHeight(g$heights[loc$t],"mm",valueOnly=TRUE)
  stopifnot(abs(pw-panel_mm)<1e-9,abs(ph-panel_mm)<1e-9)
  dimensions[[id]] <- data.frame(plot_id=id,panel_width_mm=pw,panel_height_mm=ph,
    outer_width_in=sum(common_widths),outer_height_in=sum(common_heights),font="Arial")
}
dev.off();unlink(measurement)
slot_width <- sum(common_widths);slot_height <- sum(common_heights)
draw_page <- function(ids,ncol,nrow) {
  stopifnot(length(ids)==ncol*nrow)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol,
    widths=unit(rep(slot_width,ncol),"in"),heights=unit(rep(slot_height,nrow),"in"))))
  for (i in seq_along(ids)) {
    pushViewport(viewport(layout.pos.row=(i-1L)%/%ncol+1L,layout.pos.col=(i-1L)%%ncol+1L))
    grid.draw(grobs[[ids[[i]]]])
    upViewport()
  }
  upViewport()
}
outputs <- list()
save_page <- function(ids,name,ncol=1L,nrow=1L) {
  w <- slot_width*ncol;h <- slot_height*nrow
  message(format(Sys.time(),"%F %T")," Saving ",name," (",round(w,2)," x ",round(h,2)," in)")
  cairo_pdf(file.path(task,"output/pdf",paste0(name,".pdf")),width=w,height=h,family="Arial")
  draw_page(ids,ncol,nrow);dev.off()
  ragg::agg_png(file.path(task,"output/png",paste0(name,".png")),width=w,height=h,units="in",res=240,background="white")
  draw_page(ids,ncol,nrow);dev.off()
  outputs[[name]] <<- data.frame(name=name,n_panels=length(ids),n_columns=ncol,n_rows=nrow,
    width_in=w,height_in=h,pdf=paste0("output/pdf/",name,".pdf"),png=paste0("output/png/",name,".png"))
}
for (id in c(cell_ids,marker_ids)) save_page(id,id)
save_page(cell_ids,"Six_study_celltype_UMAPs_Arial_square",3L,2L)
save_page(marker_ids[1:9],"Marker_UMAPs_3x3_page1_Chia_Fullard_Grubman",3L,3L)
save_page(marker_ids[10:18],"Marker_UMAPs_3x3_page2_Jakel_Su_Yao",3L,3L)
write.csv(do.call(rbind,coordinate_audit),file.path(task,"audit/coordinate_and_annotation_preservation.csv"),row.names=FALSE)
write.csv(do.call(rbind,scale_audit),file.path(task,"audit/marker_scale_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,label_audit),file.path(task,"audit/celltype_labels_and_anchors.csv"),row.names=FALSE)
write.csv(do.call(rbind,palette_audit),file.path(task,"audit/celltype_counts_and_colours.csv"),row.names=FALSE)
write.csv(do.call(rbind,dimensions),file.path(task,"audit/physical_panel_dimensions.csv"),row.names=FALSE)
write.csv(do.call(rbind,outputs),file.path(task,"audit/output_manifest.csv"),row.names=FALSE)
saveRDS(views,file.path(task,"audit/shared_square_view_limits.rds"))
writeLines(capture.output(sessionInfo()),file.path(task,"audit/plotting_sessionInfo.txt"))
message("UNIFORM_UMAP_PLOTTING_COMPLETE: 27 PDF and PNG files; 24 equal 114.3-mm square panels")

# S2D  Jakel P2RY12 / CX3CR1 / CSF1R UMAPs
# Source: source/project/Reference/combine/results/source_celltype_umap_20260920/revisions/20261022_Arial_square_celltype_and_marker_UMAPs/code/plot_uniform_umaps.R
#!/usr/bin/env Rscript
# Display-only revision: exact selected inputs, selected marker bars, Arial,
# one shared square view per study and an identical physical data panel.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(scattermore)
  library(scales); library(grid); library(gtable); library(digest)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task <- normalizePath(args[[1]], mustWork = TRUE)
for (d in c("output/pdf", "output/png", "audit", "runtime"))
  dir.create(file.path(task, d), recursive = TRUE, showWarnings = FALSE)
studies <- c("Chia", "Fullard", "Grubman", "Jakel", "Su", "Yao")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
annotation <- c(Chia="cell.type", Fullard="celltype", Grubman="cellType",
                Jakel="CellType", Su="CellType", Yao="AuthorMajor")
expected_n <- c(Chia=87470L, Fullard=68557L, Grubman=13214L, Jakel=17799L,
                Su=224464L, Yao=73721L)
font <- systemfonts::match_fonts("Arial")
stopifnot(grepl("Arial", basename(font$path), fixed = TRUE))
writeLines(font$path, file.path(task, "audit/Arial_font_path.txt"))
base_colours <- read.csv(file.path(task,"inputs/source_tables/base_celltype_colours.csv"), stringsAsFactors=FALSE)
base_anchors <- read.csv(file.path(task,"inputs/source_tables/base_celltype_anchors.csv"), stringsAsFactors=FALSE)
yao_colours <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_colours.csv"), stringsAsFactors=FALSE)
yao_anchors <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_anchors.csv"), stringsAsFactors=FALSE)
bars <- read.csv(file.path(task,"audit/selected_marker_bar_provenance.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(bars)==18L, all(toupper(as.character(bars$exact_match))=="TRUE"),
          all(bars$upper_limit[bars$study %in% c("Grubman","Jakel")] == 1),
          all(bars$upper_limit[bars$study == "Fullard"] == 2))
theme_uniform <- theme_classic(base_size=10, base_family="Arial") + theme(
  text=element_text(family="Arial"),
  plot.title=element_text(face="bold",size=12,margin=margin(b=3)),
  plot.subtitle=element_text(size=8,colour="#333333",margin=margin(b=7)),
  axis.title=element_text(size=10), axis.text=element_text(size=9,colour="black"),
  axis.line=element_line(linewidth=.45), axis.ticks=element_line(linewidth=.35),
  legend.title=element_text(size=8),legend.text=element_text(size=8),
  legend.position="right",legend.box.margin=margin(0,0,0,5),
  plot.margin=margin(8,8,8,8),plot.title.position="panel")
plots <- list(); cell_ids <- character(); marker_ids <- character()
coordinate_audit <- list(); scale_audit <- list(); label_audit <- list()
palette_audit <- list(); views <- list()
for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Preparing ",study)
  ct <- readRDS(file.path(task,"inputs/celltype",paste0(study,".rds")))
  mb <- readRDS(file.path(task,"inputs/marker",paste0(study,".rds")))
  stopifnot(length(ct$cells)==expected_n[[study]], !anyDuplicated(ct$cells),
    identical(rownames(ct$umap),ct$cells),identical(rownames(ct$metadata),ct$cells),
    identical(rownames(mb$umap),mb$cells),identical(colnames(mb$expression),mb$cells),
    identical(colnames(mb$counts),mb$cells),setequal(ct$cells,mb$cells))
  align <- match(mb$cells,ct$cells)
  ct_aligned <- ct$umap[align,,drop=FALSE]
  # Matrix subsetting drops optional scaled:center attributes. Compare every
  # coordinate exactly after validating dimensions and cell alignment.
  stopifnot(identical(dim(ct_aligned),dim(mb$umap)),
    identical(as.numeric(ct_aligned),as.numeric(mb$umap)))
  coords <- ct$umap
  xr <- range(coords[,1]); yr <- range(coords[,2])
  span <- max(diff(xr),diff(yr))*1.16
  xlim <- mean(xr)+c(-.5,.5)*span; ylim <- mean(yr)+c(-.5,.5)*span
  views[[study]] <- list(xlim=xlim,ylim=ylim)
  values <- as.character(ct$metadata[[annotation[[study]]]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  if (study=="Yao") {
    pc <- subset(yao_colours,plot_id=="Yao_current_UMAP_author_labels")
    anchors <- subset(yao_anchors,plot_id=="Yao_current_UMAP_author_labels")
    names(pc)[names(pc)=="author_Major"] <- "celltype"
    names(anchors)[names(anchors)=="type"] <- "celltype"
  } else {
    pc <- base_colours[base_colours$plot_id==study,,drop=FALSE]
    anchors <- base_anchors[base_anchors$plot_id==study,,drop=FALSE]
  }
  types <- pc$celltype; palette <- setNames(pc$colour,types)
  counts <- table(factor(values,levels=types))
  stopifnot(setequal(unique(values),types),identical(as.integer(counts),as.integer(pc$n_cells)),
    identical(as.integer(counts[anchors$celltype]),as.integer(anchors$n_cells)),
    all(anchors$anchor_cell %in% ct$cells))
  anchor_pos <- match(anchors$anchor_cell,ct$cells)
  stopifnot(max(abs(coords[anchor_pos,1]-anchors$x))<1e-12,
            max(abs(coords[anchor_pos,2]-anchors$y))<1e-12)
  data <- data.frame(x=coords[,1],y=coords[,2],celltype=factor(values,levels=types),row.names=ct$cells)
  subtitle <- if (study=="Yao") "73,721 cells; 16 author cell types; 174 unmatched" else
    paste0(comma(length(values))," cells; ",length(types)," cell types")
  title <- if(study=="Yao") "Yao | author cell types" else paste(study,"| cell types")
  p <- ggplot(data,aes(x=x,y=y,colour=celltype))+
    geom_scattermore(pointsize=2.2,pixels=c(1400,1400),show.legend=FALSE)+
    scale_colour_manual(values=palette,drop=FALSE)+
    geom_text_repel(data=anchors,aes(x=x,y=y,label=label_text),inherit.aes=FALSE,
      family="Arial",colour="black",size=if(length(types)>20) 2.35 else 2.75,
      lineheight=.95,seed=20261016,nudge_x=anchors$nudge_x,nudge_y=anchors$nudge_y,
      box.padding=unit(.35,"lines"),point.padding=unit(.10,"lines"),
      min.segment.length=0,segment.colour="#454545",segment.size=.20,
      force=2,force_pull=.3,max.overlaps=Inf,max.time=Inf,max.iter=20000,
      show.legend=FALSE)+
    coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
    labs(title=title,subtitle=subtitle,x="UMAP_1",y="UMAP_2")+
    theme_uniform+theme(legend.position="none")
  id <- paste0(study,"_celltype_UMAP")
  plots[[id]] <- p; cell_ids <- c(cell_ids,id)
  anchors$study <- study
  label_audit[[study]] <- anchors[,c("study","celltype","n_cells","x","y","anchor_cell","label_text","nudge_x","nudge_y")]
  palette_audit[[study]] <- data.frame(study=study,celltype=types,n_cells=as.integer(counts),colour=pc$colour)
  coordinate_audit[[study]] <- data.frame(study=study,n_cells=length(values),annotation_column=annotation[[study]],
    n_celltype_labels=length(types),celltype_and_marker_barcode_sets_identical=TRUE,
    celltype_and_marker_coordinates_identical_by_barcode=TRUE,
    coordinates_sha256=digest(coords,algo="sha256"),
    celltype_values_sha256=digest(values,algo="sha256"),
    marker_expression_sha256=digest(mb$expression,algo="sha256"),
    original_x_min=xr[1],original_x_max=xr[2],original_y_min=yr[1],original_y_max=yr[2],
    view_x_min=xlim[1],view_x_max=xlim[2],view_y_min=ylim[1],view_y_max=ylim[2],
    coordinate_ratio=1,umap_recomputed=FALSE,annotation_modified=FALSE)
  thin <- CreateSeuratObject(counts=mb$counts,min.cells=0,min.features=0,project=study)
  thin <- SetAssayData(thin,assay="RNA",layer="data",new.data=mb$expression)
  mc <- mb$umap; colnames(mc) <- c("UMAP_1","UMAP_2")
  thin[["umap"]] <- CreateDimReducObject(embeddings=mc,key="UMAP_",assay="RNA")
  stopifnot(identical(dim(Embeddings(thin,"umap")),dim(mb$umap)),
    identical(as.numeric(Embeddings(thin,"umap")),as.numeric(mb$umap)),
    identical(as.matrix(GetAssayData(thin,assay="RNA",layer="data")),as.matrix(mb$expression)))
  for (gene in markers) {
    bar <- bars[bars$study==study & bars$marker==gene,,drop=FALSE]
    stopifnot(nrow(bar)==1L)
    feature <- mb$features[[gene]]; upper <- bar$upper_limit
    tick_values <- as.numeric(strsplit(bar$tick_values,";",fixed=TRUE)[[1]])
    tick_labels <- strsplit(bar$tick_labels,";",fixed=TRUE)[[1]]
    vals <- as.numeric(mb$expression[feature,])
    p <- FeaturePlot(thin,features=feature,reduction="umap",dims=c(1,2),
      cols=c("#D3D3D3","#E31A1C"),order=TRUE,pt.size=1.3,
      min.cutoff=0,max.cutoff=upper,slot="data",raster=TRUE,raster.dpi=c(1200,1200),combine=FALSE)[[1]]+
      scale_colour_gradient(low=bar$colour_low,high=bar$colour_high,
        limits=c(0,upper),oob=scales::squish,breaks=tick_values,labels=tick_labels,
        guide=guide_colourbar(title="Log-normalized\nRNA",barheight=unit(2.6,"cm"),barwidth=unit(.30,"cm")))+
      coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
      labs(title=paste(study,"|",gene),subtitle=paste0(comma(length(vals))," cells"),x="UMAP_1",y="UMAP_2")+
      theme_uniform
    # FeaturePlot must retain the cached values after applying only the chosen cap.
    j <- match(rownames(p$data),mb$cells)
    stopifnot(!anyNA(j),nrow(p$data)==length(mb$cells),!anyDuplicated(rownames(p$data)),
      setequal(rownames(p$data),mb$cells),
      identical(as.numeric(p$data[[feature]]),pmin(pmax(vals[j],0),upper)))
    id <- paste0(study,"_",gene,"_FeaturePlot")
    plots[[id]] <- p; marker_ids <- c(marker_ids,id)
    scale_audit[[id]] <- data.frame(study=study,marker=gene,lower_limit=0,upper_limit=upper,
      visible_breaks=bar$tick_values,visible_tick_labels=bar$tick_labels,
      colour_low=bar$colour_low,colour_high=bar$colour_high,n_cells=length(vals),
      n_positive=sum(vals>0),n_saturated=sum(vals>upper),selected_pdf_sha256=bar$sha256,
      expression_identical_before_display_cap=TRUE,sorted_ascending=TRUE,
      panel_width_in=4.5,panel_height_in=4.5,coordinate_ratio=1)
  }
  rm(ct,mb,thin); gc(verbose=FALSE)
}

# Measure on the actual Cairo/Arial PDF device. Fix the panel row and column,
# then equalize every surrounding row/column so no legend can shrink a panel.
measurement <- file.path(task,"runtime/layout_measurement.pdf")
cairo_pdf(measurement,width=9,height=9,family="Arial")
grobs <- lapply(plots,ggplotGrob)
panel_mm <- 114.3 # 4.5 inches, identical in all 24 panels
shape <- vapply(grobs,function(g) paste(length(g$widths),length(g$heights)),character(1))
stopifnot(length(unique(shape))==1L)
panel_locations <- vapply(grobs,function(g) {
  loc <- g$layout[g$layout$name=="panel",]
  paste(loc$l,loc$r,loc$t,loc$b,collapse=";")
},character(1))
stopifnot(length(unique(panel_locations))==1L)
dimensions <- list()
for (id in names(grobs)) {
  g <- grobs[[id]]; loc <- g$layout[g$layout$name=="panel",]
  stopifnot(nrow(loc)==1L,loc$l==loc$r,loc$t==loc$b)
  g$widths[loc$l] <- unit(panel_mm,"mm");g$heights[loc$t] <- unit(panel_mm,"mm")
  grobs[[id]] <- g
}
all_widths <- lapply(grobs,function(g) convertWidth(g$widths,"in",valueOnly=TRUE))
all_heights <- lapply(grobs,function(g) convertHeight(g$heights,"in",valueOnly=TRUE))
common_widths <- do.call(pmax,all_widths)
common_heights <- do.call(pmax,all_heights)
for (id in names(grobs)) {
  g <- grobs[[id]];g$widths <- unit(common_widths,"in");g$heights <- unit(common_heights,"in")
  grobs[[id]] <- g
  loc <- g$layout[g$layout$name=="panel",]
  pw <- convertWidth(g$widths[loc$l],"mm",valueOnly=TRUE)
  ph <- convertHeight(g$heights[loc$t],"mm",valueOnly=TRUE)
  stopifnot(abs(pw-panel_mm)<1e-9,abs(ph-panel_mm)<1e-9)
  dimensions[[id]] <- data.frame(plot_id=id,panel_width_mm=pw,panel_height_mm=ph,
    outer_width_in=sum(common_widths),outer_height_in=sum(common_heights),font="Arial")
}
dev.off();unlink(measurement)
slot_width <- sum(common_widths);slot_height <- sum(common_heights)
draw_page <- function(ids,ncol,nrow) {
  stopifnot(length(ids)==ncol*nrow)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol,
    widths=unit(rep(slot_width,ncol),"in"),heights=unit(rep(slot_height,nrow),"in"))))
  for (i in seq_along(ids)) {
    pushViewport(viewport(layout.pos.row=(i-1L)%/%ncol+1L,layout.pos.col=(i-1L)%%ncol+1L))
    grid.draw(grobs[[ids[[i]]]])
    upViewport()
  }
  upViewport()
}
outputs <- list()
save_page <- function(ids,name,ncol=1L,nrow=1L) {
  w <- slot_width*ncol;h <- slot_height*nrow
  message(format(Sys.time(),"%F %T")," Saving ",name," (",round(w,2)," x ",round(h,2)," in)")
  cairo_pdf(file.path(task,"output/pdf",paste0(name,".pdf")),width=w,height=h,family="Arial")
  draw_page(ids,ncol,nrow);dev.off()
  ragg::agg_png(file.path(task,"output/png",paste0(name,".png")),width=w,height=h,units="in",res=240,background="white")
  draw_page(ids,ncol,nrow);dev.off()
  outputs[[name]] <<- data.frame(name=name,n_panels=length(ids),n_columns=ncol,n_rows=nrow,
    width_in=w,height_in=h,pdf=paste0("output/pdf/",name,".pdf"),png=paste0("output/png/",name,".png"))
}
for (id in c(cell_ids,marker_ids)) save_page(id,id)
save_page(cell_ids,"Six_study_celltype_UMAPs_Arial_square",3L,2L)
save_page(marker_ids[1:9],"Marker_UMAPs_3x3_page1_Chia_Fullard_Grubman",3L,3L)
save_page(marker_ids[10:18],"Marker_UMAPs_3x3_page2_Jakel_Su_Yao",3L,3L)
write.csv(do.call(rbind,coordinate_audit),file.path(task,"audit/coordinate_and_annotation_preservation.csv"),row.names=FALSE)
write.csv(do.call(rbind,scale_audit),file.path(task,"audit/marker_scale_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,label_audit),file.path(task,"audit/celltype_labels_and_anchors.csv"),row.names=FALSE)
write.csv(do.call(rbind,palette_audit),file.path(task,"audit/celltype_counts_and_colours.csv"),row.names=FALSE)
write.csv(do.call(rbind,dimensions),file.path(task,"audit/physical_panel_dimensions.csv"),row.names=FALSE)
write.csv(do.call(rbind,outputs),file.path(task,"audit/output_manifest.csv"),row.names=FALSE)
saveRDS(views,file.path(task,"audit/shared_square_view_limits.rds"))
writeLines(capture.output(sessionInfo()),file.path(task,"audit/plotting_sessionInfo.txt"))
message("UNIFORM_UMAP_PLOTTING_COMPLETE: 27 PDF and PNG files; 24 equal 114.3-mm square panels")

# S2E  Su P2RY12 / CX3CR1 / CSF1R UMAPs
# Source: source/project/Reference/combine/results/source_celltype_umap_20260920/revisions/20261022_Arial_square_celltype_and_marker_UMAPs/code/plot_uniform_umaps.R
#!/usr/bin/env Rscript
# Display-only revision: exact selected inputs, selected marker bars, Arial,
# one shared square view per study and an identical physical data panel.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(scattermore)
  library(scales); library(grid); library(gtable); library(digest)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task <- normalizePath(args[[1]], mustWork = TRUE)
for (d in c("output/pdf", "output/png", "audit", "runtime"))
  dir.create(file.path(task, d), recursive = TRUE, showWarnings = FALSE)
studies <- c("Chia", "Fullard", "Grubman", "Jakel", "Su", "Yao")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
annotation <- c(Chia="cell.type", Fullard="celltype", Grubman="cellType",
                Jakel="CellType", Su="CellType", Yao="AuthorMajor")
expected_n <- c(Chia=87470L, Fullard=68557L, Grubman=13214L, Jakel=17799L,
                Su=224464L, Yao=73721L)
font <- systemfonts::match_fonts("Arial")
stopifnot(grepl("Arial", basename(font$path), fixed = TRUE))
writeLines(font$path, file.path(task, "audit/Arial_font_path.txt"))
base_colours <- read.csv(file.path(task,"inputs/source_tables/base_celltype_colours.csv"), stringsAsFactors=FALSE)
base_anchors <- read.csv(file.path(task,"inputs/source_tables/base_celltype_anchors.csv"), stringsAsFactors=FALSE)
yao_colours <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_colours.csv"), stringsAsFactors=FALSE)
yao_anchors <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_anchors.csv"), stringsAsFactors=FALSE)
bars <- read.csv(file.path(task,"audit/selected_marker_bar_provenance.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(bars)==18L, all(toupper(as.character(bars$exact_match))=="TRUE"),
          all(bars$upper_limit[bars$study %in% c("Grubman","Jakel")] == 1),
          all(bars$upper_limit[bars$study == "Fullard"] == 2))
theme_uniform <- theme_classic(base_size=10, base_family="Arial") + theme(
  text=element_text(family="Arial"),
  plot.title=element_text(face="bold",size=12,margin=margin(b=3)),
  plot.subtitle=element_text(size=8,colour="#333333",margin=margin(b=7)),
  axis.title=element_text(size=10), axis.text=element_text(size=9,colour="black"),
  axis.line=element_line(linewidth=.45), axis.ticks=element_line(linewidth=.35),
  legend.title=element_text(size=8),legend.text=element_text(size=8),
  legend.position="right",legend.box.margin=margin(0,0,0,5),
  plot.margin=margin(8,8,8,8),plot.title.position="panel")
plots <- list(); cell_ids <- character(); marker_ids <- character()
coordinate_audit <- list(); scale_audit <- list(); label_audit <- list()
palette_audit <- list(); views <- list()
for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Preparing ",study)
  ct <- readRDS(file.path(task,"inputs/celltype",paste0(study,".rds")))
  mb <- readRDS(file.path(task,"inputs/marker",paste0(study,".rds")))
  stopifnot(length(ct$cells)==expected_n[[study]], !anyDuplicated(ct$cells),
    identical(rownames(ct$umap),ct$cells),identical(rownames(ct$metadata),ct$cells),
    identical(rownames(mb$umap),mb$cells),identical(colnames(mb$expression),mb$cells),
    identical(colnames(mb$counts),mb$cells),setequal(ct$cells,mb$cells))
  align <- match(mb$cells,ct$cells)
  ct_aligned <- ct$umap[align,,drop=FALSE]
  # Matrix subsetting drops optional scaled:center attributes. Compare every
  # coordinate exactly after validating dimensions and cell alignment.
  stopifnot(identical(dim(ct_aligned),dim(mb$umap)),
    identical(as.numeric(ct_aligned),as.numeric(mb$umap)))
  coords <- ct$umap
  xr <- range(coords[,1]); yr <- range(coords[,2])
  span <- max(diff(xr),diff(yr))*1.16
  xlim <- mean(xr)+c(-.5,.5)*span; ylim <- mean(yr)+c(-.5,.5)*span
  views[[study]] <- list(xlim=xlim,ylim=ylim)
  values <- as.character(ct$metadata[[annotation[[study]]]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  if (study=="Yao") {
    pc <- subset(yao_colours,plot_id=="Yao_current_UMAP_author_labels")
    anchors <- subset(yao_anchors,plot_id=="Yao_current_UMAP_author_labels")
    names(pc)[names(pc)=="author_Major"] <- "celltype"
    names(anchors)[names(anchors)=="type"] <- "celltype"
  } else {
    pc <- base_colours[base_colours$plot_id==study,,drop=FALSE]
    anchors <- base_anchors[base_anchors$plot_id==study,,drop=FALSE]
  }
  types <- pc$celltype; palette <- setNames(pc$colour,types)
  counts <- table(factor(values,levels=types))
  stopifnot(setequal(unique(values),types),identical(as.integer(counts),as.integer(pc$n_cells)),
    identical(as.integer(counts[anchors$celltype]),as.integer(anchors$n_cells)),
    all(anchors$anchor_cell %in% ct$cells))
  anchor_pos <- match(anchors$anchor_cell,ct$cells)
  stopifnot(max(abs(coords[anchor_pos,1]-anchors$x))<1e-12,
            max(abs(coords[anchor_pos,2]-anchors$y))<1e-12)
  data <- data.frame(x=coords[,1],y=coords[,2],celltype=factor(values,levels=types),row.names=ct$cells)
  subtitle <- if (study=="Yao") "73,721 cells; 16 author cell types; 174 unmatched" else
    paste0(comma(length(values))," cells; ",length(types)," cell types")
  title <- if(study=="Yao") "Yao | author cell types" else paste(study,"| cell types")
  p <- ggplot(data,aes(x=x,y=y,colour=celltype))+
    geom_scattermore(pointsize=2.2,pixels=c(1400,1400),show.legend=FALSE)+
    scale_colour_manual(values=palette,drop=FALSE)+
    geom_text_repel(data=anchors,aes(x=x,y=y,label=label_text),inherit.aes=FALSE,
      family="Arial",colour="black",size=if(length(types)>20) 2.35 else 2.75,
      lineheight=.95,seed=20261016,nudge_x=anchors$nudge_x,nudge_y=anchors$nudge_y,
      box.padding=unit(.35,"lines"),point.padding=unit(.10,"lines"),
      min.segment.length=0,segment.colour="#454545",segment.size=.20,
      force=2,force_pull=.3,max.overlaps=Inf,max.time=Inf,max.iter=20000,
      show.legend=FALSE)+
    coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
    labs(title=title,subtitle=subtitle,x="UMAP_1",y="UMAP_2")+
    theme_uniform+theme(legend.position="none")
  id <- paste0(study,"_celltype_UMAP")
  plots[[id]] <- p; cell_ids <- c(cell_ids,id)
  anchors$study <- study
  label_audit[[study]] <- anchors[,c("study","celltype","n_cells","x","y","anchor_cell","label_text","nudge_x","nudge_y")]
  palette_audit[[study]] <- data.frame(study=study,celltype=types,n_cells=as.integer(counts),colour=pc$colour)
  coordinate_audit[[study]] <- data.frame(study=study,n_cells=length(values),annotation_column=annotation[[study]],
    n_celltype_labels=length(types),celltype_and_marker_barcode_sets_identical=TRUE,
    celltype_and_marker_coordinates_identical_by_barcode=TRUE,
    coordinates_sha256=digest(coords,algo="sha256"),
    celltype_values_sha256=digest(values,algo="sha256"),
    marker_expression_sha256=digest(mb$expression,algo="sha256"),
    original_x_min=xr[1],original_x_max=xr[2],original_y_min=yr[1],original_y_max=yr[2],
    view_x_min=xlim[1],view_x_max=xlim[2],view_y_min=ylim[1],view_y_max=ylim[2],
    coordinate_ratio=1,umap_recomputed=FALSE,annotation_modified=FALSE)
  thin <- CreateSeuratObject(counts=mb$counts,min.cells=0,min.features=0,project=study)
  thin <- SetAssayData(thin,assay="RNA",layer="data",new.data=mb$expression)
  mc <- mb$umap; colnames(mc) <- c("UMAP_1","UMAP_2")
  thin[["umap"]] <- CreateDimReducObject(embeddings=mc,key="UMAP_",assay="RNA")
  stopifnot(identical(dim(Embeddings(thin,"umap")),dim(mb$umap)),
    identical(as.numeric(Embeddings(thin,"umap")),as.numeric(mb$umap)),
    identical(as.matrix(GetAssayData(thin,assay="RNA",layer="data")),as.matrix(mb$expression)))
  for (gene in markers) {
    bar <- bars[bars$study==study & bars$marker==gene,,drop=FALSE]
    stopifnot(nrow(bar)==1L)
    feature <- mb$features[[gene]]; upper <- bar$upper_limit
    tick_values <- as.numeric(strsplit(bar$tick_values,";",fixed=TRUE)[[1]])
    tick_labels <- strsplit(bar$tick_labels,";",fixed=TRUE)[[1]]
    vals <- as.numeric(mb$expression[feature,])
    p <- FeaturePlot(thin,features=feature,reduction="umap",dims=c(1,2),
      cols=c("#D3D3D3","#E31A1C"),order=TRUE,pt.size=1.3,
      min.cutoff=0,max.cutoff=upper,slot="data",raster=TRUE,raster.dpi=c(1200,1200),combine=FALSE)[[1]]+
      scale_colour_gradient(low=bar$colour_low,high=bar$colour_high,
        limits=c(0,upper),oob=scales::squish,breaks=tick_values,labels=tick_labels,
        guide=guide_colourbar(title="Log-normalized\nRNA",barheight=unit(2.6,"cm"),barwidth=unit(.30,"cm")))+
      coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
      labs(title=paste(study,"|",gene),subtitle=paste0(comma(length(vals))," cells"),x="UMAP_1",y="UMAP_2")+
      theme_uniform
    # FeaturePlot must retain the cached values after applying only the chosen cap.
    j <- match(rownames(p$data),mb$cells)
    stopifnot(!anyNA(j),nrow(p$data)==length(mb$cells),!anyDuplicated(rownames(p$data)),
      setequal(rownames(p$data),mb$cells),
      identical(as.numeric(p$data[[feature]]),pmin(pmax(vals[j],0),upper)))
    id <- paste0(study,"_",gene,"_FeaturePlot")
    plots[[id]] <- p; marker_ids <- c(marker_ids,id)
    scale_audit[[id]] <- data.frame(study=study,marker=gene,lower_limit=0,upper_limit=upper,
      visible_breaks=bar$tick_values,visible_tick_labels=bar$tick_labels,
      colour_low=bar$colour_low,colour_high=bar$colour_high,n_cells=length(vals),
      n_positive=sum(vals>0),n_saturated=sum(vals>upper),selected_pdf_sha256=bar$sha256,
      expression_identical_before_display_cap=TRUE,sorted_ascending=TRUE,
      panel_width_in=4.5,panel_height_in=4.5,coordinate_ratio=1)
  }
  rm(ct,mb,thin); gc(verbose=FALSE)
}

# Measure on the actual Cairo/Arial PDF device. Fix the panel row and column,
# then equalize every surrounding row/column so no legend can shrink a panel.
measurement <- file.path(task,"runtime/layout_measurement.pdf")
cairo_pdf(measurement,width=9,height=9,family="Arial")
grobs <- lapply(plots,ggplotGrob)
panel_mm <- 114.3 # 4.5 inches, identical in all 24 panels
shape <- vapply(grobs,function(g) paste(length(g$widths),length(g$heights)),character(1))
stopifnot(length(unique(shape))==1L)
panel_locations <- vapply(grobs,function(g) {
  loc <- g$layout[g$layout$name=="panel",]
  paste(loc$l,loc$r,loc$t,loc$b,collapse=";")
},character(1))
stopifnot(length(unique(panel_locations))==1L)
dimensions <- list()
for (id in names(grobs)) {
  g <- grobs[[id]]; loc <- g$layout[g$layout$name=="panel",]
  stopifnot(nrow(loc)==1L,loc$l==loc$r,loc$t==loc$b)
  g$widths[loc$l] <- unit(panel_mm,"mm");g$heights[loc$t] <- unit(panel_mm,"mm")
  grobs[[id]] <- g
}
all_widths <- lapply(grobs,function(g) convertWidth(g$widths,"in",valueOnly=TRUE))
all_heights <- lapply(grobs,function(g) convertHeight(g$heights,"in",valueOnly=TRUE))
common_widths <- do.call(pmax,all_widths)
common_heights <- do.call(pmax,all_heights)
for (id in names(grobs)) {
  g <- grobs[[id]];g$widths <- unit(common_widths,"in");g$heights <- unit(common_heights,"in")
  grobs[[id]] <- g
  loc <- g$layout[g$layout$name=="panel",]
  pw <- convertWidth(g$widths[loc$l],"mm",valueOnly=TRUE)
  ph <- convertHeight(g$heights[loc$t],"mm",valueOnly=TRUE)
  stopifnot(abs(pw-panel_mm)<1e-9,abs(ph-panel_mm)<1e-9)
  dimensions[[id]] <- data.frame(plot_id=id,panel_width_mm=pw,panel_height_mm=ph,
    outer_width_in=sum(common_widths),outer_height_in=sum(common_heights),font="Arial")
}
dev.off();unlink(measurement)
slot_width <- sum(common_widths);slot_height <- sum(common_heights)
draw_page <- function(ids,ncol,nrow) {
  stopifnot(length(ids)==ncol*nrow)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol,
    widths=unit(rep(slot_width,ncol),"in"),heights=unit(rep(slot_height,nrow),"in"))))
  for (i in seq_along(ids)) {
    pushViewport(viewport(layout.pos.row=(i-1L)%/%ncol+1L,layout.pos.col=(i-1L)%%ncol+1L))
    grid.draw(grobs[[ids[[i]]]])
    upViewport()
  }
  upViewport()
}
outputs <- list()
save_page <- function(ids,name,ncol=1L,nrow=1L) {
  w <- slot_width*ncol;h <- slot_height*nrow
  message(format(Sys.time(),"%F %T")," Saving ",name," (",round(w,2)," x ",round(h,2)," in)")
  cairo_pdf(file.path(task,"output/pdf",paste0(name,".pdf")),width=w,height=h,family="Arial")
  draw_page(ids,ncol,nrow);dev.off()
  ragg::agg_png(file.path(task,"output/png",paste0(name,".png")),width=w,height=h,units="in",res=240,background="white")
  draw_page(ids,ncol,nrow);dev.off()
  outputs[[name]] <<- data.frame(name=name,n_panels=length(ids),n_columns=ncol,n_rows=nrow,
    width_in=w,height_in=h,pdf=paste0("output/pdf/",name,".pdf"),png=paste0("output/png/",name,".png"))
}
for (id in c(cell_ids,marker_ids)) save_page(id,id)
save_page(cell_ids,"Six_study_celltype_UMAPs_Arial_square",3L,2L)
save_page(marker_ids[1:9],"Marker_UMAPs_3x3_page1_Chia_Fullard_Grubman",3L,3L)
save_page(marker_ids[10:18],"Marker_UMAPs_3x3_page2_Jakel_Su_Yao",3L,3L)
write.csv(do.call(rbind,coordinate_audit),file.path(task,"audit/coordinate_and_annotation_preservation.csv"),row.names=FALSE)
write.csv(do.call(rbind,scale_audit),file.path(task,"audit/marker_scale_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,label_audit),file.path(task,"audit/celltype_labels_and_anchors.csv"),row.names=FALSE)
write.csv(do.call(rbind,palette_audit),file.path(task,"audit/celltype_counts_and_colours.csv"),row.names=FALSE)
write.csv(do.call(rbind,dimensions),file.path(task,"audit/physical_panel_dimensions.csv"),row.names=FALSE)
write.csv(do.call(rbind,outputs),file.path(task,"audit/output_manifest.csv"),row.names=FALSE)
saveRDS(views,file.path(task,"audit/shared_square_view_limits.rds"))
writeLines(capture.output(sessionInfo()),file.path(task,"audit/plotting_sessionInfo.txt"))
message("UNIFORM_UMAP_PLOTTING_COMPLETE: 27 PDF and PNG files; 24 equal 114.3-mm square panels")

# S2F  Yao P2RY12 / CX3CR1 / CSF1R UMAPs
# Source: source/project/Reference/combine/results/source_celltype_umap_20260920/revisions/20261022_Arial_square_celltype_and_marker_UMAPs/code/plot_uniform_umaps.R
#!/usr/bin/env Rscript
# Display-only revision: exact selected inputs, selected marker bars, Arial,
# one shared square view per study and an identical physical data panel.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggrepel); library(scattermore)
  library(scales); library(grid); library(gtable); library(digest)
})
options(Seurat.object.assay.version = "v5")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
task <- normalizePath(args[[1]], mustWork = TRUE)
for (d in c("output/pdf", "output/png", "audit", "runtime"))
  dir.create(file.path(task, d), recursive = TRUE, showWarnings = FALSE)
studies <- c("Chia", "Fullard", "Grubman", "Jakel", "Su", "Yao")
markers <- c("P2RY12", "CX3CR1", "CSF1R")
annotation <- c(Chia="cell.type", Fullard="celltype", Grubman="cellType",
                Jakel="CellType", Su="CellType", Yao="AuthorMajor")
expected_n <- c(Chia=87470L, Fullard=68557L, Grubman=13214L, Jakel=17799L,
                Su=224464L, Yao=73721L)
font <- systemfonts::match_fonts("Arial")
stopifnot(grepl("Arial", basename(font$path), fixed = TRUE))
writeLines(font$path, file.path(task, "audit/Arial_font_path.txt"))
base_colours <- read.csv(file.path(task,"inputs/source_tables/base_celltype_colours.csv"), stringsAsFactors=FALSE)
base_anchors <- read.csv(file.path(task,"inputs/source_tables/base_celltype_anchors.csv"), stringsAsFactors=FALSE)
yao_colours <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_colours.csv"), stringsAsFactors=FALSE)
yao_anchors <- read.csv(file.path(task,"inputs/source_tables/yao_celltype_anchors.csv"), stringsAsFactors=FALSE)
bars <- read.csv(file.path(task,"audit/selected_marker_bar_provenance.csv"), stringsAsFactors=FALSE)
stopifnot(nrow(bars)==18L, all(toupper(as.character(bars$exact_match))=="TRUE"),
          all(bars$upper_limit[bars$study %in% c("Grubman","Jakel")] == 1),
          all(bars$upper_limit[bars$study == "Fullard"] == 2))
theme_uniform <- theme_classic(base_size=10, base_family="Arial") + theme(
  text=element_text(family="Arial"),
  plot.title=element_text(face="bold",size=12,margin=margin(b=3)),
  plot.subtitle=element_text(size=8,colour="#333333",margin=margin(b=7)),
  axis.title=element_text(size=10), axis.text=element_text(size=9,colour="black"),
  axis.line=element_line(linewidth=.45), axis.ticks=element_line(linewidth=.35),
  legend.title=element_text(size=8),legend.text=element_text(size=8),
  legend.position="right",legend.box.margin=margin(0,0,0,5),
  plot.margin=margin(8,8,8,8),plot.title.position="panel")
plots <- list(); cell_ids <- character(); marker_ids <- character()
coordinate_audit <- list(); scale_audit <- list(); label_audit <- list()
palette_audit <- list(); views <- list()
for (study in studies) {
  message(format(Sys.time(), "%F %T"), " Preparing ",study)
  ct <- readRDS(file.path(task,"inputs/celltype",paste0(study,".rds")))
  mb <- readRDS(file.path(task,"inputs/marker",paste0(study,".rds")))
  stopifnot(length(ct$cells)==expected_n[[study]], !anyDuplicated(ct$cells),
    identical(rownames(ct$umap),ct$cells),identical(rownames(ct$metadata),ct$cells),
    identical(rownames(mb$umap),mb$cells),identical(colnames(mb$expression),mb$cells),
    identical(colnames(mb$counts),mb$cells),setequal(ct$cells,mb$cells))
  align <- match(mb$cells,ct$cells)
  ct_aligned <- ct$umap[align,,drop=FALSE]
  # Matrix subsetting drops optional scaled:center attributes. Compare every
  # coordinate exactly after validating dimensions and cell alignment.
  stopifnot(identical(dim(ct_aligned),dim(mb$umap)),
    identical(as.numeric(ct_aligned),as.numeric(mb$umap)))
  coords <- ct$umap
  xr <- range(coords[,1]); yr <- range(coords[,2])
  span <- max(diff(xr),diff(yr))*1.16
  xlim <- mean(xr)+c(-.5,.5)*span; ylim <- mean(yr)+c(-.5,.5)*span
  views[[study]] <- list(xlim=xlim,ylim=ylim)
  values <- as.character(ct$metadata[[annotation[[study]]]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  if (study=="Yao") {
    pc <- subset(yao_colours,plot_id=="Yao_current_UMAP_author_labels")
    anchors <- subset(yao_anchors,plot_id=="Yao_current_UMAP_author_labels")
    names(pc)[names(pc)=="author_Major"] <- "celltype"
    names(anchors)[names(anchors)=="type"] <- "celltype"
  } else {
    pc <- base_colours[base_colours$plot_id==study,,drop=FALSE]
    anchors <- base_anchors[base_anchors$plot_id==study,,drop=FALSE]
  }
  types <- pc$celltype; palette <- setNames(pc$colour,types)
  counts <- table(factor(values,levels=types))
  stopifnot(setequal(unique(values),types),identical(as.integer(counts),as.integer(pc$n_cells)),
    identical(as.integer(counts[anchors$celltype]),as.integer(anchors$n_cells)),
    all(anchors$anchor_cell %in% ct$cells))
  anchor_pos <- match(anchors$anchor_cell,ct$cells)
  stopifnot(max(abs(coords[anchor_pos,1]-anchors$x))<1e-12,
            max(abs(coords[anchor_pos,2]-anchors$y))<1e-12)
  data <- data.frame(x=coords[,1],y=coords[,2],celltype=factor(values,levels=types),row.names=ct$cells)
  subtitle <- if (study=="Yao") "73,721 cells; 16 author cell types; 174 unmatched" else
    paste0(comma(length(values))," cells; ",length(types)," cell types")
  title <- if(study=="Yao") "Yao | author cell types" else paste(study,"| cell types")
  p <- ggplot(data,aes(x=x,y=y,colour=celltype))+
    geom_scattermore(pointsize=2.2,pixels=c(1400,1400),show.legend=FALSE)+
    scale_colour_manual(values=palette,drop=FALSE)+
    geom_text_repel(data=anchors,aes(x=x,y=y,label=label_text),inherit.aes=FALSE,
      family="Arial",colour="black",size=if(length(types)>20) 2.35 else 2.75,
      lineheight=.95,seed=20261016,nudge_x=anchors$nudge_x,nudge_y=anchors$nudge_y,
      box.padding=unit(.35,"lines"),point.padding=unit(.10,"lines"),
      min.segment.length=0,segment.colour="#454545",segment.size=.20,
      force=2,force_pull=.3,max.overlaps=Inf,max.time=Inf,max.iter=20000,
      show.legend=FALSE)+
    coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
    labs(title=title,subtitle=subtitle,x="UMAP_1",y="UMAP_2")+
    theme_uniform+theme(legend.position="none")
  id <- paste0(study,"_celltype_UMAP")
  plots[[id]] <- p; cell_ids <- c(cell_ids,id)
  anchors$study <- study
  label_audit[[study]] <- anchors[,c("study","celltype","n_cells","x","y","anchor_cell","label_text","nudge_x","nudge_y")]
  palette_audit[[study]] <- data.frame(study=study,celltype=types,n_cells=as.integer(counts),colour=pc$colour)
  coordinate_audit[[study]] <- data.frame(study=study,n_cells=length(values),annotation_column=annotation[[study]],
    n_celltype_labels=length(types),celltype_and_marker_barcode_sets_identical=TRUE,
    celltype_and_marker_coordinates_identical_by_barcode=TRUE,
    coordinates_sha256=digest(coords,algo="sha256"),
    celltype_values_sha256=digest(values,algo="sha256"),
    marker_expression_sha256=digest(mb$expression,algo="sha256"),
    original_x_min=xr[1],original_x_max=xr[2],original_y_min=yr[1],original_y_max=yr[2],
    view_x_min=xlim[1],view_x_max=xlim[2],view_y_min=ylim[1],view_y_max=ylim[2],
    coordinate_ratio=1,umap_recomputed=FALSE,annotation_modified=FALSE)
  thin <- CreateSeuratObject(counts=mb$counts,min.cells=0,min.features=0,project=study)
  thin <- SetAssayData(thin,assay="RNA",layer="data",new.data=mb$expression)
  mc <- mb$umap; colnames(mc) <- c("UMAP_1","UMAP_2")
  thin[["umap"]] <- CreateDimReducObject(embeddings=mc,key="UMAP_",assay="RNA")
  stopifnot(identical(dim(Embeddings(thin,"umap")),dim(mb$umap)),
    identical(as.numeric(Embeddings(thin,"umap")),as.numeric(mb$umap)),
    identical(as.matrix(GetAssayData(thin,assay="RNA",layer="data")),as.matrix(mb$expression)))
  for (gene in markers) {
    bar <- bars[bars$study==study & bars$marker==gene,,drop=FALSE]
    stopifnot(nrow(bar)==1L)
    feature <- mb$features[[gene]]; upper <- bar$upper_limit
    tick_values <- as.numeric(strsplit(bar$tick_values,";",fixed=TRUE)[[1]])
    tick_labels <- strsplit(bar$tick_labels,";",fixed=TRUE)[[1]]
    vals <- as.numeric(mb$expression[feature,])
    p <- FeaturePlot(thin,features=feature,reduction="umap",dims=c(1,2),
      cols=c("#D3D3D3","#E31A1C"),order=TRUE,pt.size=1.3,
      min.cutoff=0,max.cutoff=upper,slot="data",raster=TRUE,raster.dpi=c(1200,1200),combine=FALSE)[[1]]+
      scale_colour_gradient(low=bar$colour_low,high=bar$colour_high,
        limits=c(0,upper),oob=scales::squish,breaks=tick_values,labels=tick_labels,
        guide=guide_colourbar(title="Log-normalized\nRNA",barheight=unit(2.6,"cm"),barwidth=unit(.30,"cm")))+
      coord_fixed(ratio=1,xlim=xlim,ylim=ylim,expand=FALSE,clip="off")+
      labs(title=paste(study,"|",gene),subtitle=paste0(comma(length(vals))," cells"),x="UMAP_1",y="UMAP_2")+
      theme_uniform
    # FeaturePlot must retain the cached values after applying only the chosen cap.
    j <- match(rownames(p$data),mb$cells)
    stopifnot(!anyNA(j),nrow(p$data)==length(mb$cells),!anyDuplicated(rownames(p$data)),
      setequal(rownames(p$data),mb$cells),
      identical(as.numeric(p$data[[feature]]),pmin(pmax(vals[j],0),upper)))
    id <- paste0(study,"_",gene,"_FeaturePlot")
    plots[[id]] <- p; marker_ids <- c(marker_ids,id)
    scale_audit[[id]] <- data.frame(study=study,marker=gene,lower_limit=0,upper_limit=upper,
      visible_breaks=bar$tick_values,visible_tick_labels=bar$tick_labels,
      colour_low=bar$colour_low,colour_high=bar$colour_high,n_cells=length(vals),
      n_positive=sum(vals>0),n_saturated=sum(vals>upper),selected_pdf_sha256=bar$sha256,
      expression_identical_before_display_cap=TRUE,sorted_ascending=TRUE,
      panel_width_in=4.5,panel_height_in=4.5,coordinate_ratio=1)
  }
  rm(ct,mb,thin); gc(verbose=FALSE)
}

# Measure on the actual Cairo/Arial PDF device. Fix the panel row and column,
# then equalize every surrounding row/column so no legend can shrink a panel.
measurement <- file.path(task,"runtime/layout_measurement.pdf")
cairo_pdf(measurement,width=9,height=9,family="Arial")
grobs <- lapply(plots,ggplotGrob)
panel_mm <- 114.3 # 4.5 inches, identical in all 24 panels
shape <- vapply(grobs,function(g) paste(length(g$widths),length(g$heights)),character(1))
stopifnot(length(unique(shape))==1L)
panel_locations <- vapply(grobs,function(g) {
  loc <- g$layout[g$layout$name=="panel",]
  paste(loc$l,loc$r,loc$t,loc$b,collapse=";")
},character(1))
stopifnot(length(unique(panel_locations))==1L)
dimensions <- list()
for (id in names(grobs)) {
  g <- grobs[[id]]; loc <- g$layout[g$layout$name=="panel",]
  stopifnot(nrow(loc)==1L,loc$l==loc$r,loc$t==loc$b)
  g$widths[loc$l] <- unit(panel_mm,"mm");g$heights[loc$t] <- unit(panel_mm,"mm")
  grobs[[id]] <- g
}
all_widths <- lapply(grobs,function(g) convertWidth(g$widths,"in",valueOnly=TRUE))
all_heights <- lapply(grobs,function(g) convertHeight(g$heights,"in",valueOnly=TRUE))
common_widths <- do.call(pmax,all_widths)
common_heights <- do.call(pmax,all_heights)
for (id in names(grobs)) {
  g <- grobs[[id]];g$widths <- unit(common_widths,"in");g$heights <- unit(common_heights,"in")
  grobs[[id]] <- g
  loc <- g$layout[g$layout$name=="panel",]
  pw <- convertWidth(g$widths[loc$l],"mm",valueOnly=TRUE)
  ph <- convertHeight(g$heights[loc$t],"mm",valueOnly=TRUE)
  stopifnot(abs(pw-panel_mm)<1e-9,abs(ph-panel_mm)<1e-9)
  dimensions[[id]] <- data.frame(plot_id=id,panel_width_mm=pw,panel_height_mm=ph,
    outer_width_in=sum(common_widths),outer_height_in=sum(common_heights),font="Arial")
}
dev.off();unlink(measurement)
slot_width <- sum(common_widths);slot_height <- sum(common_heights)
draw_page <- function(ids,ncol,nrow) {
  stopifnot(length(ids)==ncol*nrow)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol,
    widths=unit(rep(slot_width,ncol),"in"),heights=unit(rep(slot_height,nrow),"in"))))
  for (i in seq_along(ids)) {
    pushViewport(viewport(layout.pos.row=(i-1L)%/%ncol+1L,layout.pos.col=(i-1L)%%ncol+1L))
    grid.draw(grobs[[ids[[i]]]])
    upViewport()
  }
  upViewport()
}
outputs <- list()
save_page <- function(ids,name,ncol=1L,nrow=1L) {
  w <- slot_width*ncol;h <- slot_height*nrow
  message(format(Sys.time(),"%F %T")," Saving ",name," (",round(w,2)," x ",round(h,2)," in)")
  cairo_pdf(file.path(task,"output/pdf",paste0(name,".pdf")),width=w,height=h,family="Arial")
  draw_page(ids,ncol,nrow);dev.off()
  ragg::agg_png(file.path(task,"output/png",paste0(name,".png")),width=w,height=h,units="in",res=240,background="white")
  draw_page(ids,ncol,nrow);dev.off()
  outputs[[name]] <<- data.frame(name=name,n_panels=length(ids),n_columns=ncol,n_rows=nrow,
    width_in=w,height_in=h,pdf=paste0("output/pdf/",name,".pdf"),png=paste0("output/png/",name,".png"))
}
for (id in c(cell_ids,marker_ids)) save_page(id,id)
save_page(cell_ids,"Six_study_celltype_UMAPs_Arial_square",3L,2L)
save_page(marker_ids[1:9],"Marker_UMAPs_3x3_page1_Chia_Fullard_Grubman",3L,3L)
save_page(marker_ids[10:18],"Marker_UMAPs_3x3_page2_Jakel_Su_Yao",3L,3L)
write.csv(do.call(rbind,coordinate_audit),file.path(task,"audit/coordinate_and_annotation_preservation.csv"),row.names=FALSE)
write.csv(do.call(rbind,scale_audit),file.path(task,"audit/marker_scale_parameters.csv"),row.names=FALSE)
write.csv(do.call(rbind,label_audit),file.path(task,"audit/celltype_labels_and_anchors.csv"),row.names=FALSE)
write.csv(do.call(rbind,palette_audit),file.path(task,"audit/celltype_counts_and_colours.csv"),row.names=FALSE)
write.csv(do.call(rbind,dimensions),file.path(task,"audit/physical_panel_dimensions.csv"),row.names=FALSE)
write.csv(do.call(rbind,outputs),file.path(task,"audit/output_manifest.csv"),row.names=FALSE)
saveRDS(views,file.path(task,"audit/shared_square_view_limits.rds"))
writeLines(capture.output(sessionInfo()),file.path(task,"audit/plotting_sessionInfo.txt"))
message("UNIFORM_UMAP_PLOTTING_COMPLETE: 27 PDF and PNG files; 24 equal 114.3-mm square panels")


# Supplementary Figure 3 ------------------------------------------------------------

# S3A  Human MG UMAP split by Condition
# Source: source/figure_settings/S3A_Condition.R
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


# S3B  Human MG UMAP split by Region
# Source: source/figure_settings/S3B_Region.R
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


# S3C  Nuclei and mean genes by Study
# Source: source/figure_settings/S3C_Study.R
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


# S3D  Nuclei and mean genes by Region
# Source: source/figure_settings/S3D_Region.R
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


# S3E  Condition x MG0-MG5 ROE heatmap
# Source: source/figure_settings/S3E_Condition.R
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


# S3F  Region x MG0-MG5 ROE heatmap
# Source: source/figure_settings/S3F_Region.R
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


# S3G  NFIX SOX2 SOX10 SREBF1 expression and AUCell scores
# Source: source/project/code_upload/SCENIC.r
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


# Figure 2 ------------------------------------------------------------

# F2A  Condition and healthy-age MG0-MG5 100% composition
# Source: source/project/Reference/combine/scripts/age_MG0_MG5_composition_disease_and_Health_age_from_original_module4.R
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

# F2B  23 myelin/OLM-associated genes across MG0-MG5
# Source: source/project/Reference/combine/scripts/20260520_dotplot.r
## 2.1 Parameters --------------------------------------------------------------
cluster_col <- 'CellType1'      
top_n       <- 5                            

# Auto-filenames
res_suffix  <- gsub("integrated_snn_", "", cluster_col)
png_path    <- paste0("results/figure/new_dotplot.png")
pdf_path    <- paste0("results/figure/new_dotplot.pdf")

## 2.2 Process Data ------------------------------------------------------------
tmp_obj <- qread('data/processed_data/mg_cca_meta_v2.qs')

genes_list <- c("TRIM2","RNF220","MBP","MOG",
                "ST18","SLC44A1","PHLPP1","PTPRD","NCAM2","PLP1","NCAM1","KCNMB4",
                "GNAO1","UGT8","MOBP","MYRF","OPALIN","QKI","CLDN11","FA2H","DLG1","SH3TC2","PTN")

## 2.3 Plotting ----------------------------------------------------------------
DefaultAssay(tmp_obj) <- 'RNA'
p_base <- DotPlot(tmp_obj, features = genes_list, group.by = cluster_col)

p_dot <- p_base$data %>%
  dplyr::rename(pct_exp = pct.exp, avg_exp = avg.exp.scaled, Gene = features.plot, Cluster = id) %>%
  # 1. 交换 X 和 Y 的映射
  ggplot(aes(x = Gene, y = Cluster, size = pct_exp, color = avg_exp)) + 
  geom_point() +
  scale_color_gradient2(low = "#000000", mid = "#007965", high = "#ffcc29", midpoint = 0) +
  # 2. 将 scale_y 换成 scale_x，并去掉 rev()，让基因按列表顺序从左到右排列
  scale_x_discrete(limits = genes_list) + 
  theme_classic() +
  labs(title = paste0("Genes Dotplot"), x = NULL, y = NULL) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    # 3. 将 vjust = 1 改为 vjust = 0.5，这样旋转 90 度后的基因名能更好地对齐刻度线居中
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.margin = margin(t = 50, r = 10, b = 10, l = 10)
  )

## 2.4 Save --------------------------------------------------------------------
# 4. 对调保存时的 width 和 height
ggsave(png_path, p_dot, width = 6, height = 3.5, limitsize = FALSE)
ggsave(pdf_path, p_dot, width = 6, height = 3.5, limitsize = FALSE)

# F2C  TRIM2 RNF220 MBP MOG FeaturePlots
# Source: source/project/Reference/combine/scripts/20260520_featureplot.r
input_file      <- 'data/processed_data/mg_cca_meta_v2.qs'
output_fig_dir  <- 'results/figure/'
save_path   <- paste0(output_fig_dir, '20260520_featurePlot.pdf')

dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

gene_list <- c('TRIM2','RNF220','MBP','MOG')

plot_pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(10)
tmp_obj <- qread(input_file)

plot <- FeaturePlot(tmp_obj,features = gene_list,cols = plot_pal, raster = FALSE)
ggsave(save_path, plot, width = 7.5, height = 6, limitsize = FALSE)

# F2D  Nine-term OLM GSEA
# Source: source/project/code_upload/GSEA.r
# GSEA ----------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(qs)
library(openxlsx)
library(future)
library(pheatmap)
library(ggplotify)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(patchwork)
library(stringr)
library(GseaVis)

input_file      <- 'data/processed_data/mg_cca_meta_v2.qs'
dir_data        <- "scripts/GSEA/data/"
dir_res         <- "scripts/GSEA/results/"
custom_set_file <- paste0(dir_data, 'my_gene_set.xlsx')
go_gmt_file     <- paste0(dir_data, "c5.go.v2023.2.Hs.symbols.gmt")

target_group_by <- 'CellType1'
target_cluster  <- 'MG0'  

output_qs       <- paste0(dir_data, 'GSEA_kk_', target_cluster, '.qs')
output_txt      <- paste0(dir_res, target_cluster, "_combined_genesets_GSEA_result.txt")
output_pdf_p    <- paste0(dir_res, target_cluster, "_combined_GSEA_select_plot.v2.pdf")
output_pdf_nop  <- paste0(dir_res, target_cluster, "_combined_GSEA_select_plot_no_pvalue.v2.pdf")

dir.create(dir_res, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_data, recursive = TRUE, showWarnings = FALSE)

target_go_terms <- c(
  'GOBP_ENSHEATHMENT_OF_NEURONS',
  'GOBP_OLIGODENDROCYTE_DEVELOPMENT_AND_MYELINATION',
  'GOCC_PLASMA_MEMBRANE_REGION',
  'GOCC_MYELIN_SHEATH',
  'GOBP_OLIGODENDROCYTE_DIFFERENTIATION',
  'GOBP_AXON_DEVELOPMENT',
  'GOCC_NEURON_TO_NEURON_SYNAPSE',
  'GOCC_CELL_CELL_JUNCTION',
  'GOCC_MICROTUBULE'
)

selected_plot_ids <- c(
  'ensheathment of neurons',
  'oligodendrocyte development and myelination',
  'gocc plasma membrane region',
  'gocc myelin sheath',
  'oligodendrocyte differentiation',
  'axon development',
  'gocc neuron to neuron synapse',
  'gocc cell cell junction',
  'gocc microtubule'
)

plot_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
  "#FF7F00", "#A65628", "#F781BF", "#1B9E77", "#666666"
)

if(!file.exists(input_file)) stop("Input file not found.")
tmp_obj <- qread(input_file)

gene_sets_raw <- read.xlsx(custom_set_file, colNames = FALSE, na.strings = "")
temp_gmt_path <- paste0(dir_data, 'temp_my_genesets.gmt')
write.table(gene_sets_raw, temp_gmt_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

custom_gmt <- read.gmt(temp_gmt_path)
custom_gmt$gene <- toupper(custom_gmt$gene)

gmt_go <- read.gmt(go_gmt_file)
filtered_go <- gmt_go %>%
  filter(term %in% target_go_terms) %>%
  mutate(term = gsub("GOBP_", "", term)) %>%
  mutate(term = tolower(gsub("_", " ", term)))

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
rt_sorted <- deg_genes %>%
  arrange(desc(avg_log2FC)) %>%
  distinct(gene, .keep_all = TRUE)

logFC_vector <- rt_sorted$avg_log2FC
names(logFC_vector) <- toupper(rt_sorted$gene)

kk <- GSEA(
  geneList      = logFC_vector, 
  TERM2GENE     = combined_gmt, 
  pvalueCutoff  = 1,
  pAdjustMethod = 'BH',
  minGSSize     = 1,
  maxGSSize     = 500,
  scoreType     = "pos",
  eps           = 0,
  verbose       = FALSE
)

# Fig.01 ---------------
p_gsea_pval <- gseaNb(
  object    = kk, 
  geneSetID = selected_plot_ids,
  subPlot   = 3, 
  curveCol  = plot_colors,
  addPval   = TRUE,
  pvalSize  = 5,
  base_size = 15
)

# Fig.02 ---------------
p_gsea_clean <- gseaNb(
  object    = kk, 
  geneSetID = selected_plot_ids,
  subPlot   = 3, 
  curveCol  = plot_colors,
  addPval   = FALSE,
  pvalSize  = 5,
  base_size = 15
)

qsave(kk, file = output_qs)

write.table(as.data.frame(kk), 
            file = output_txt, 
            sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file = output_pdf_p, width = 12, height = 8)
print(p_gsea_pval)
dev.off()

pdf(file = output_pdf_nop, width = 12, height = 8)
print(p_gsea_clean)
dev.off()

if(file.exists(temp_gmt_path)) unlink(temp_gmt_path)


# Supplementary Figure 4 ------------------------------------------------------------

# S4A  Human OLM top-50 module score over all Human/Monkey/Mouse microglia, positive-only yellow/red UMAP
# Source: source/project/Reference/combine/scripts/MultiSpecies/OLM_MG0_top50_cross_species_20260915/scripts/01_human_MG0_top50_score_and_cross_species_OLM_violin.R
#!/usr/bin/env Rscript

# Cross-species OLM analysis using the integrated Human/Monkey/Mouse object.
#
# Biological definitions, retained exactly from the source combine_3.r script:
#   Human OLM       = Species == "Human"  & CellType1 == "MG0"
#   Monkey OLM-like = Species == "Monkey" & integrated_snn_res.0.1 == "1"
#   Mouse OLM-like  = Species == "Mouse"  & integrated_snn_res.0.1 == "1"
#
# The Human MG0 top-50 signature is calculated in the RNA assay from
# Human MG0 versus Human MG1--MG5. It replaces the earlier externally supplied
# reference gene set. Every cell is scored using AddModuleScore(ctrl = 10), as
# in the supplied spatial_module.r reference. UMAP coordinates are read only.

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/OLM_MG0_top50_cross_species_20260915")
input_qs <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/02.cross_species_integrated_noref_Human_Mouse_Monkey.qs")
result_dir <- file.path(project_dir, "results", "01_human_MG0_top50_score_cross_species_OLM")
figure_dir <- file.path(result_dir, "figures")
table_dir <- file.path(result_dir, "tables")
object_dir <- file.path(result_dir, "objects")
log_dir <- file.path(project_dir, "logs")
tmp_dir <- file.path(project_dir, "tmp")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

# All temporary/cache locations remain in the project rather than the home directory.
Sys.setenv(
  TMPDIR = tmp_dir,
  TMP = tmp_dir,
  TEMP = tmp_dir,
  XDG_CACHE_HOME = file.path(project_dir, "cache"),
  R_USER_CACHE_DIR = file.path(project_dir, "cache", "R")
)
dir.create(Sys.getenv("XDG_CACHE_HOME"), recursive = TRUE, showWarnings = FALSE)
dir.create(Sys.getenv("R_USER_CACHE_DIR"), recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library"), file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(SCP)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260915)

if (!file.exists(input_qs)) stop("Input object does not exist: ", input_qs)

obj <- qread(input_qs)
required_metadata <- c("Species", "CellType1", "integrated_snn_res.0.1")
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata) > 0L) {
  stop("Input object lacks required metadata: ", paste(missing_metadata, collapse = ", "))
}
if (!"RNA" %in% Assays(obj)) stop("Input object lacks an RNA assay.")
if (!"umap" %in% Reductions(obj)) stop("Input object lacks the stored UMAP reduction.")
if (!all(c("Human", "Monkey", "Mouse") %in% unique(as.character(obj$Species)))) {
  stop("All Human, Monkey and Mouse cells must be present.")
}
DefaultAssay(obj) <- "RNA"

# 1. Human MG0 FindMarkers: only Human cells enter the reference contrast.
human_cells <- colnames(obj)[as.character(obj$Species) == "Human"]
human_obj <- subset(obj, cells = human_cells)
DefaultAssay(human_obj) <- "RNA"
Idents(human_obj) <- "CellType1"

human_levels <- levels(Idents(human_obj))
required_human_types <- c("MG0", "MG1", "MG2", "MG3", "MG4", "MG5")
if (!all(required_human_types %in% human_levels)) {
  stop("Human CellType1 is missing one or more MG0--MG5 labels.")
}

human_mg0_markers <- FindMarkers(
  object = human_obj,
  ident.1 = "MG0",
  ident.2 = c("MG1", "MG2", "MG3", "MG4", "MG5"),
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0
)

human_mg0_markers$gene <- rownames(human_mg0_markers)
human_mg0_markers <- human_mg0_markers[
  order(human_mg0_markers$p_val_adj, -human_mg0_markers$avg_log2FC, human_mg0_markers$gene),
  c("gene", setdiff(colnames(human_mg0_markers), "gene")),
  drop = FALSE
]
if (nrow(human_mg0_markers) < 50L) {
  stop("Fewer than 50 positive Human MG0 markers were returned.")
}

human_mg0_top50 <- human_mg0_markers[seq_len(50L), , drop = FALSE]
top50_genes <- human_mg0_top50$gene
if (length(unique(top50_genes)) != 50L) stop("The top-50 signature contains duplicate genes.")
if (length(setdiff(top50_genes, rownames(obj[["RNA"]]))) > 0L) {
  stop("At least one Human MG0 top-50 gene is absent from the integrated RNA assay.")
}

write.table(
  human_mg0_markers,
  file = file.path(table_dir, "Human_MG0_vs_Human_MG1toMG5_FindMarkers_all_positive.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.csv(
  human_mg0_markers,
  file = file.path(table_dir, "Human_MG0_vs_Human_MG1toMG5_FindMarkers_all_positive.csv"),
  row.names = FALSE
)
write.table(
  human_mg0_top50,
  file = file.path(table_dir, "Human_MG0_FindMarkers_top50_reference_signature.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.csv(
  human_mg0_top50,
  file = file.path(table_dir, "Human_MG0_FindMarkers_top50_reference_signature.csv"),
  row.names = FALSE
)

# 2. Score every cell from the identical Human MG0 top-50 signature.
# This is intentionally the same AddModuleScore control-gene design as spatial_module.r.
obj <- AddModuleScore(
  object = obj,
  assay = "RNA",
  features = list(top50_genes),
  ctrl = 10,
  name = "Human_MG0_top50_score",
  seed = 20260915
)
score_column_from_AddModuleScore <- "Human_MG0_top50_score1"
if (!score_column_from_AddModuleScore %in% colnames(obj[[]])) stop("AddModuleScore did not produce ", score_column_from_AddModuleScore)
# A descriptive alias prevents the automatic trailing "1" from appearing in the final figures.
score_column <- "Human_MG0_top50_module_score"
obj[[score_column]] <- obj[[score_column_from_AddModuleScore]][, 1]

# 3. Apply the author-defined cross-species OLM/MG0-like assignments.
obj$OLM_cross_species_definition <- "Other MG"
obj$OLM_cross_species_definition[
  as.character(obj$Species) == "Human" & as.character(obj$CellType1) == "MG0"
] <- "OLM"
obj$OLM_cross_species_definition[
  as.character(obj$Species) == "Monkey" & as.character(obj$integrated_snn_res.0.1) == "1"
] <- "OLM"
obj$OLM_cross_species_definition[
  as.character(obj$Species) == "Mouse" & as.character(obj$integrated_snn_res.0.1) == "1"
] <- "OLM"
obj$OLM_cross_species_definition <- factor(obj$OLM_cross_species_definition, levels = c("OLM", "Other MG"))

obj$OLM_definition_rule <- "Other MG"
obj$OLM_definition_rule[
  as.character(obj$Species) == "Human" & as.character(obj$CellType1) == "MG0"
] <- "Human: CellType1 == MG0"
obj$OLM_definition_rule[
  as.character(obj$Species) == "Monkey" & as.character(obj$integrated_snn_res.0.1) == "1"
] <- "Monkey: integrated_snn_res.0.1 == 1"
obj$OLM_definition_rule[
  as.character(obj$Species) == "Mouse" & as.character(obj$integrated_snn_res.0.1) == "1"
] <- "Mouse: integrated_snn_res.0.1 == 1"

# 4. Stored integrated UMAP.  This is the exact positive-expression drawing
# scheme used for the prior LPC/CPZ panels: all cells are grey first, then only
# score-positive cells are drawn in ascending score order from yellow to dark red.
score_umap_coordinates <- Embeddings(obj, reduction = "umap")
score_plot_frame <- data.frame(
  UMAP_1 = score_umap_coordinates[, 1],
  UMAP_2 = score_umap_coordinates[, 2],
  module_score = as.numeric(obj[[score_column]][, 1]),
  stringsAsFactors = FALSE
)
score_positive_frame <- score_plot_frame[score_plot_frame$module_score > 0, , drop = FALSE]
score_positive_frame <- score_positive_frame[order(score_positive_frame$module_score), , drop = FALSE]
score_upper_limit <- if (nrow(score_positive_frame) == 0L) {
  1
} else {
  as.numeric(quantile(score_positive_frame$module_score, 0.99, names = FALSE))
}
score_umap_all <- ggplot(score_plot_frame, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    color = "#B8B8B8",
    size = 0.20,
    alpha = 0.60,
    shape = 16,
    stroke = 0
  ) +
  geom_point(
    data = score_positive_frame,
    aes(color = module_score),
    size = 0.20 * 1.55,
    alpha = 1,
    shape = 16,
    stroke = 0
  ) +
  scale_color_gradient(
    low = "#FFD166",
    high = "#B2182B",
    limits = c(0, score_upper_limit),
    oob = scales::squish
  ) +
  labs(
    title = "Human MG0 top-50 module score",
    subtitle = "Human MG0 signature scored across all Human, Monkey and Mouse cells",
    color = "Module\nscore",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  coord_equal() +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )
ggsave(file.path(figure_dir, "01_Human_MG0_top50_module_score_all_species_UMAP.png"), score_umap_all, width = 7.2, height = 6.4, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "01_Human_MG0_top50_module_score_all_species_UMAP.pdf"), score_umap_all, width = 7.2, height = 6.4, device = grDevices::cairo_pdf, bg = "white")

# S4B  Within each species: OLM top-50 module-score violin, OLM versus Other MG
# Source: source/project/Reference/combine/scripts/MultiSpecies/OLM_MG0_top50_cross_species_20260915/scripts/01_human_MG0_top50_score_and_cross_species_OLM_violin.R
#!/usr/bin/env Rscript

# Cross-species OLM analysis using the integrated Human/Monkey/Mouse object.
#
# Biological definitions, retained exactly from the source combine_3.r script:
#   Human OLM       = Species == "Human"  & CellType1 == "MG0"
#   Monkey OLM-like = Species == "Monkey" & integrated_snn_res.0.1 == "1"
#   Mouse OLM-like  = Species == "Mouse"  & integrated_snn_res.0.1 == "1"
#
# The Human MG0 top-50 signature is calculated in the RNA assay from
# Human MG0 versus Human MG1--MG5. It replaces the earlier externally supplied
# reference gene set. Every cell is scored using AddModuleScore(ctrl = 10), as
# in the supplied spatial_module.r reference. UMAP coordinates are read only.

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/OLM_MG0_top50_cross_species_20260915")
input_qs <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/02.cross_species_integrated_noref_Human_Mouse_Monkey.qs")
result_dir <- file.path(project_dir, "results", "01_human_MG0_top50_score_cross_species_OLM")
figure_dir <- file.path(result_dir, "figures")
table_dir <- file.path(result_dir, "tables")
object_dir <- file.path(result_dir, "objects")
log_dir <- file.path(project_dir, "logs")
tmp_dir <- file.path(project_dir, "tmp")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

# All temporary/cache locations remain in the project rather than the home directory.
Sys.setenv(
  TMPDIR = tmp_dir,
  TMP = tmp_dir,
  TEMP = tmp_dir,
  XDG_CACHE_HOME = file.path(project_dir, "cache"),
  R_USER_CACHE_DIR = file.path(project_dir, "cache", "R")
)
dir.create(Sys.getenv("XDG_CACHE_HOME"), recursive = TRUE, showWarnings = FALSE)
dir.create(Sys.getenv("R_USER_CACHE_DIR"), recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library"), file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(SCP)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260915)

if (!file.exists(input_qs)) stop("Input object does not exist: ", input_qs)

obj <- qread(input_qs)
required_metadata <- c("Species", "CellType1", "integrated_snn_res.0.1")
missing_metadata <- setdiff(required_metadata, colnames(obj[[]]))
if (length(missing_metadata) > 0L) {
  stop("Input object lacks required metadata: ", paste(missing_metadata, collapse = ", "))
}
if (!"RNA" %in% Assays(obj)) stop("Input object lacks an RNA assay.")
if (!"umap" %in% Reductions(obj)) stop("Input object lacks the stored UMAP reduction.")
if (!all(c("Human", "Monkey", "Mouse") %in% unique(as.character(obj$Species)))) {
  stop("All Human, Monkey and Mouse cells must be present.")
}
DefaultAssay(obj) <- "RNA"

# 1. Human MG0 FindMarkers: only Human cells enter the reference contrast.
human_cells <- colnames(obj)[as.character(obj$Species) == "Human"]
human_obj <- subset(obj, cells = human_cells)
DefaultAssay(human_obj) <- "RNA"
Idents(human_obj) <- "CellType1"

human_levels <- levels(Idents(human_obj))
required_human_types <- c("MG0", "MG1", "MG2", "MG3", "MG4", "MG5")
if (!all(required_human_types %in% human_levels)) {
  stop("Human CellType1 is missing one or more MG0--MG5 labels.")
}

human_mg0_markers <- FindMarkers(
  object = human_obj,
  ident.1 = "MG0",
  ident.2 = c("MG1", "MG2", "MG3", "MG4", "MG5"),
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0
)

human_mg0_markers$gene <- rownames(human_mg0_markers)
human_mg0_markers <- human_mg0_markers[
  order(human_mg0_markers$p_val_adj, -human_mg0_markers$avg_log2FC, human_mg0_markers$gene),
  c("gene", setdiff(colnames(human_mg0_markers), "gene")),
  drop = FALSE
]
if (nrow(human_mg0_markers) < 50L) {
  stop("Fewer than 50 positive Human MG0 markers were returned.")
}

human_mg0_top50 <- human_mg0_markers[seq_len(50L), , drop = FALSE]
top50_genes <- human_mg0_top50$gene
if (length(unique(top50_genes)) != 50L) stop("The top-50 signature contains duplicate genes.")
if (length(setdiff(top50_genes, rownames(obj[["RNA"]]))) > 0L) {
  stop("At least one Human MG0 top-50 gene is absent from the integrated RNA assay.")
}

write.table(
  human_mg0_markers,
  file = file.path(table_dir, "Human_MG0_vs_Human_MG1toMG5_FindMarkers_all_positive.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.csv(
  human_mg0_markers,
  file = file.path(table_dir, "Human_MG0_vs_Human_MG1toMG5_FindMarkers_all_positive.csv"),
  row.names = FALSE
)
write.table(
  human_mg0_top50,
  file = file.path(table_dir, "Human_MG0_FindMarkers_top50_reference_signature.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.csv(
  human_mg0_top50,
  file = file.path(table_dir, "Human_MG0_FindMarkers_top50_reference_signature.csv"),
  row.names = FALSE
)

# 2. Score every cell from the identical Human MG0 top-50 signature.
# This is intentionally the same AddModuleScore control-gene design as spatial_module.r.
obj <- AddModuleScore(
  object = obj,
  assay = "RNA",
  features = list(top50_genes),
  ctrl = 10,
  name = "Human_MG0_top50_score",
  seed = 20260915
)
score_column_from_AddModuleScore <- "Human_MG0_top50_score1"
if (!score_column_from_AddModuleScore %in% colnames(obj[[]])) stop("AddModuleScore did not produce ", score_column_from_AddModuleScore)
# A descriptive alias prevents the automatic trailing "1" from appearing in the final figures.
score_column <- "Human_MG0_top50_module_score"
obj[[score_column]] <- obj[[score_column_from_AddModuleScore]][, 1]

# 3. Apply the author-defined cross-species OLM/MG0-like assignments.
obj$OLM_cross_species_definition <- "Other MG"
obj$OLM_cross_species_definition[
  as.character(obj$Species) == "Human" & as.character(obj$CellType1) == "MG0"
] <- "OLM"
obj$OLM_cross_species_definition[
  as.character(obj$Species) == "Monkey" & as.character(obj$integrated_snn_res.0.1) == "1"
] <- "OLM"
obj$OLM_cross_species_definition[
  as.character(obj$Species) == "Mouse" & as.character(obj$integrated_snn_res.0.1) == "1"
] <- "OLM"
obj$OLM_cross_species_definition <- factor(obj$OLM_cross_species_definition, levels = c("OLM", "Other MG"))

obj$OLM_definition_rule <- "Other MG"
obj$OLM_definition_rule[
  as.character(obj$Species) == "Human" & as.character(obj$CellType1) == "MG0"
] <- "Human: CellType1 == MG0"
obj$OLM_definition_rule[
  as.character(obj$Species) == "Monkey" & as.character(obj$integrated_snn_res.0.1) == "1"
] <- "Monkey: integrated_snn_res.0.1 == 1"
obj$OLM_definition_rule[
  as.character(obj$Species) == "Mouse" & as.character(obj$integrated_snn_res.0.1) == "1"
] <- "Mouse: integrated_snn_res.0.1 == 1"

# 4. Stored integrated UMAP.  This is the exact positive-expression drawing
# scheme used for the prior LPC/CPZ panels: all cells are grey first, then only
# score-positive cells are drawn in ascending score order from yellow to dark red.
score_umap_coordinates <- Embeddings(obj, reduction = "umap")
score_plot_frame <- data.frame(
  UMAP_1 = score_umap_coordinates[, 1],
  UMAP_2 = score_umap_coordinates[, 2],
  module_score = as.numeric(obj[[score_column]][, 1]),
  stringsAsFactors = FALSE
)
score_positive_frame <- score_plot_frame[score_plot_frame$module_score > 0, , drop = FALSE]
score_positive_frame <- score_positive_frame[order(score_positive_frame$module_score), , drop = FALSE]
score_upper_limit <- if (nrow(score_positive_frame) == 0L) {
  1
} else {
  as.numeric(quantile(score_positive_frame$module_score, 0.99, names = FALSE))
}
score_umap_all <- ggplot(score_plot_frame, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    color = "#B8B8B8",
    size = 0.20,
    alpha = 0.60,
    shape = 16,
    stroke = 0
  ) +
  geom_point(
    data = score_positive_frame,
    aes(color = module_score),
    size = 0.20 * 1.55,
    alpha = 1,
    shape = 16,
    stroke = 0
  ) +
  scale_color_gradient(
    low = "#FFD166",
    high = "#B2182B",
    limits = c(0, score_upper_limit),
    oob = scales::squish
  ) +
  labs(
    title = "Human MG0 top-50 module score",
    subtitle = "Human MG0 signature scored across all Human, Monkey and Mouse cells",
    color = "Module\nscore",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  coord_equal() +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )
ggsave(file.path(figure_dir, "01_Human_MG0_top50_module_score_all_species_UMAP.png"), score_umap_all, width = 7.2, height = 6.4, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "01_Human_MG0_top50_module_score_all_species_UMAP.pdf"), score_umap_all, width = 7.2, height = 6.4, device = grDevices::cairo_pdf, bg = "white")

# 5. Per-species OLM versus Other MG module-score violins.
# The FeatureStatPlot settings and red/blue palette duplicate the prior M0-vs-Other reference plots.
olm_other_colours <- c("OLM" = "#E64B35", "Other MG" = "#4DBBD5")

human_score_object <- subset(obj, cells = colnames(obj)[as.character(obj$Species) == "Human"])
human_score_object$OLM_vs_Other_MG <- factor(as.character(human_score_object$OLM_cross_species_definition), levels = c("OLM", "Other MG"))
human_violin <- SCP::FeatureStatPlot(
  srt = human_score_object,
  assay = "RNA",
  slot = "data",
  group.by = "OLM_vs_Other_MG",
  stat.by = score_column,
  add_box = FALSE,
  comparisons = list(c("Other MG", "OLM")),
  bg.by = "OLM_vs_Other_MG",
  palcolor = olm_other_colours,
  bg_palcolor = olm_other_colours,
  legend.position = "none",
  xlab = "",
  theme_args = list(base_size = 10),
  ncol = 1
)
human_violin <- human_violin & labs(title = "Human") & theme(
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 11),
  plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
  plot.margin = ggplot2::margin(t = 24, r = 4, b = 4, l = 4, unit = "pt")
)
ggsave(file.path(figure_dir, "02_Human_OLM_vs_Other_MG_top50_score_violin.png"), human_violin, width = 4.6, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "02_Human_OLM_vs_Other_MG_top50_score_violin.pdf"), human_violin, width = 4.6, height = 6.2, device = grDevices::cairo_pdf, bg = "white")

monkey_score_object <- subset(obj, cells = colnames(obj)[as.character(obj$Species) == "Monkey"])
monkey_score_object$OLM_vs_Other_MG <- factor(as.character(monkey_score_object$OLM_cross_species_definition), levels = c("OLM", "Other MG"))
monkey_violin <- SCP::FeatureStatPlot(
  srt = monkey_score_object,
  assay = "RNA",
  slot = "data",
  group.by = "OLM_vs_Other_MG",
  stat.by = score_column,
  add_box = FALSE,
  comparisons = list(c("Other MG", "OLM")),
  bg.by = "OLM_vs_Other_MG",
  palcolor = olm_other_colours,
  bg_palcolor = olm_other_colours,
  legend.position = "none",
  xlab = "",
  theme_args = list(base_size = 10),
  ncol = 1
)
monkey_violin <- monkey_violin & labs(title = "Monkey") & theme(
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 11),
  plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
  plot.margin = ggplot2::margin(t = 24, r = 4, b = 4, l = 4, unit = "pt")
)
ggsave(file.path(figure_dir, "03_Monkey_OLM_vs_Other_MG_top50_score_violin.png"), monkey_violin, width = 4.6, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "03_Monkey_OLM_vs_Other_MG_top50_score_violin.pdf"), monkey_violin, width = 4.6, height = 6.2, device = grDevices::cairo_pdf, bg = "white")

mouse_score_object <- subset(obj, cells = colnames(obj)[as.character(obj$Species) == "Mouse"])
mouse_score_object$OLM_vs_Other_MG <- factor(as.character(mouse_score_object$OLM_cross_species_definition), levels = c("OLM", "Other MG"))
mouse_violin <- SCP::FeatureStatPlot(
  srt = mouse_score_object,
  assay = "RNA",
  slot = "data",
  group.by = "OLM_vs_Other_MG",
  stat.by = score_column,
  add_box = FALSE,
  comparisons = list(c("Other MG", "OLM")),
  bg.by = "OLM_vs_Other_MG",
  palcolor = olm_other_colours,
  bg_palcolor = olm_other_colours,
  legend.position = "none",
  xlab = "",
  theme_args = list(base_size = 10),
  ncol = 1
)
mouse_violin <- mouse_violin & labs(title = "Mouse") & theme(
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 11),
  plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
  plot.margin = ggplot2::margin(t = 24, r = 4, b = 4, l = 4, unit = "pt")
)
ggsave(file.path(figure_dir, "04_Mouse_OLM_vs_Other_MG_top50_score_violin.png"), mouse_violin, width = 4.6, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "04_Mouse_OLM_vs_Other_MG_top50_score_violin.pdf"), mouse_violin, width = 4.6, height = 6.2, device = grDevices::cairo_pdf, bg = "white")

# 6. Per-species statistics and auditable cell/group assignments.
human_score_values <- as.numeric(human_score_object[[score_column]][, 1])
human_groups <- as.character(human_score_object$OLM_vs_Other_MG)
monkey_score_values <- as.numeric(monkey_score_object[[score_column]][, 1])
monkey_groups <- as.character(monkey_score_object$OLM_vs_Other_MG)
mouse_score_values <- as.numeric(mouse_score_object[[score_column]][, 1])
mouse_groups <- as.character(mouse_score_object$OLM_vs_Other_MG)

human_test <- wilcox.test(human_score_values[human_groups == "OLM"], human_score_values[human_groups == "Other MG"], alternative = "two.sided", exact = FALSE)
monkey_test <- wilcox.test(monkey_score_values[monkey_groups == "OLM"], monkey_score_values[monkey_groups == "Other MG"], alternative = "two.sided", exact = FALSE)
mouse_test <- wilcox.test(mouse_score_values[mouse_groups == "OLM"], mouse_score_values[mouse_groups == "Other MG"], alternative = "two.sided", exact = FALSE)

species_statistics <- data.frame(
  Species = c("Human", "Monkey", "Mouse"),
  OLM_definition = c("CellType1 == MG0", "integrated_snn_res.0.1 == 1", "integrated_snn_res.0.1 == 1"),
  OLM_cells = c(sum(human_groups == "OLM"), sum(monkey_groups == "OLM"), sum(mouse_groups == "OLM")),
  Other_MG_cells = c(sum(human_groups == "Other MG"), sum(monkey_groups == "Other MG"), sum(mouse_groups == "Other MG")),
  OLM_mean_module_score = c(mean(human_score_values[human_groups == "OLM"]), mean(monkey_score_values[monkey_groups == "OLM"]), mean(mouse_score_values[mouse_groups == "OLM"])),
  OLM_median_module_score = c(median(human_score_values[human_groups == "OLM"]), median(monkey_score_values[monkey_groups == "OLM"]), median(mouse_score_values[mouse_groups == "OLM"])),
  Other_MG_mean_module_score = c(mean(human_score_values[human_groups == "Other MG"]), mean(monkey_score_values[monkey_groups == "Other MG"]), mean(mouse_score_values[mouse_groups == "Other MG"])),
  Other_MG_median_module_score = c(median(human_score_values[human_groups == "Other MG"]), median(monkey_score_values[monkey_groups == "Other MG"]), median(mouse_score_values[mouse_groups == "Other MG"])),
  Wilcoxon_W = c(unname(human_test$statistic), unname(monkey_test$statistic), unname(mouse_test$statistic)),
  p_value_unadjusted = c(human_test$p.value, monkey_test$p.value, mouse_test$p.value),
  stringsAsFactors = FALSE
)
species_statistics$p_value_BH_three_species <- p.adjust(species_statistics$p_value_unadjusted, method = "BH")
species_statistics$plot_significance <- ifelse(species_statistics$p_value_unadjusted <= 1e-4, "****", ifelse(species_statistics$p_value_unadjusted <= 1e-3, "***", ifelse(species_statistics$p_value_unadjusted <= 1e-2, "**", ifelse(species_statistics$p_value_unadjusted <= 0.05, "*", "ns"))))
write.table(species_statistics, file.path(table_dir, "Human_Monkey_Mouse_OLM_vs_Other_MG_top50_score_Wilcoxon.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# S4C  Monkey_MG0_like top10 GO Biological Process barplot
# Source: source/project/Reference/combine/scripts/MultiSpecies/scripts/Monkey_MG0_like_GO.r
library(qs)
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(stringr)
library(ggplot2)

# MODULE : GO Enrichment and Barplot for a Single Cluster -------------------

## 1 Parameters --------------------------------------------------------------
# -- 数据输入参数 --
input_filepath  <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/monkey_mg_markers.qs")
target_cluster  <- "Monkey_MG0_like"       # 目标提取的 cluster 名称

# -- 过滤与分析阈值 --
p_val_adj_cut   <- 0.05                   # 显著性阈值
log2FC_cut      <- 0.5                    # Log2FC 阈值
db_org          <- "org.Mmu.eg.db"         # 物种数据库: "org.Mm.eg.db" (鼠) / "org.Hs.eg.db" (人) / "org.Mmu.eg.db" (恒河猴)
ontology_type   <- "BP"                   # GO 类型: "BP", "CC", "MF" 或 "ALL"
top_n_go        <- 10                     # 图中展示的 Top GO 通路数量

# -- 绘图与美化参数 --
cluster_color   <- "#B2DF8A"              # 柱状图的填充颜色
wrap_width      <- 60                     # 文本自动换行的字符数限制

# -- 输出参数 --
out_dir         <- "results/figure/cluster_go/"  # 图片保存的文件夹路径
out_width       <- 8                      # 输出图片宽度
out_height      <- 5                      # 输出图片高度


## 2 Process Data ------------------------------------------------------------
# 动态加载物种数据库
require(db_org, character.only = TRUE)

# 1. 读取 Marker 数据框 (修复：使用 qread 真正载入数据)
marker_df <- qread(input_filepath)

# 2. 过滤目标群的 Marker
cluster_markers <- marker_df %>%
  dplyr::filter(cluster == target_cluster) %>%
  dplyr::filter(p_val_adj < p_val_adj_cut) %>% 
  dplyr::filter(avg_log2FC > log2FC_cut) %>%
  dplyr::arrange(desc(avg_log2FC))

# 3. ID 转换 (SYMBOL -> ENTREZID)
gene_conv <- bitr(cluster_markers$gene, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = db_org)

# 4. 运行 GO 富集分析
ego_cluster <- enrichGO(
  gene          = gene_conv$ENTREZID,
  OrgDb         = db_org,
  ont           = ontology_type,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# 5. 提取富集结果并准备作图数据
plot_data <- as.data.frame(ego_cluster)
dir.create(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/table"), recursive = TRUE)
write.xlsx(plot_data,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/table/monkey_go.xlsx"))

if(nrow(plot_data) > 0) {
  plot_data <- plot_data %>%
    dplyr::arrange(p.adjust) %>%
    head(top_n_go) %>%
    dplyr::mutate(
      logP = -log10(p.adjust),
      # 对过长的 Description 自动换行
      Description = stringr::str_wrap(Description, width = wrap_width) 
    )
  
  # 锁定因子顺序，确保显著性最高的排在最上面
  plot_data$Description <- factor(plot_data$Description, levels = rev(plot_data$Description))
} else {
  warning(paste("Warning: No significant GO terms found for", target_cluster))
}


## 3 Generate Plots ----------------------------------------------------------
# 初始化画图对象为空
p_barplot <- NULL

if(nrow(plot_data) > 0) {
  
  # 计算文字的 X 轴起始坐标 (放在坐标轴起点偏右一点点，确保整齐对齐)
  text_x_pos <- max(plot_data$logP) * 0.02 
  
  p_barplot <- ggplot(plot_data, aes(x = logP, y = Description)) +
    geom_col(width = 0.8, alpha = 0.8, fill = cluster_color, color = "black", linewidth = 0.3) +
    # 文字标签：固定在图表左侧对齐
    geom_text(aes(x = text_x_pos, label = Description), 
              size = 4.5, hjust = 0, color = "black") + 
    theme_minimal() +
    labs(
      title = paste0("Cluster ", target_cluster, " Top ", top_n_go, " GO Terms (", ontology_type, ")"), 
      x = "-log10(p.adjust)", 
      y = NULL
    ) +
    theme(
      axis.text.y  = element_blank(), # 隐藏原生 Y 轴文字
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(), # 隐藏网格线，画面更干净
      axis.line    = element_line(color = "black"),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    # 向右延展 X 轴空间，防止特别长的通路名字越界被裁切
    scale_x_continuous(expand = expansion(mult = c(0, 0.4))) 
}


## 4 Save Results ------------------------------------------------------------
if(!is.null(p_barplot)) {
  # 如果文件夹不存在，自动创建
  if(!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # 替换 cluster 名字里的非法字符（如果有），防止保存报错
  safe_name <- gsub("[^a-zA-Z0-9_]", "_", target_cluster) 
  file_base <- file.path(out_dir, paste0("Cluster_", safe_name, "_GO_", ontology_type))
  
  # 保存 PDF
  ggsave(filename = paste0(file_base, ".pdf"), 
         plot     = p_barplot, 
         width    = out_width, 
         height   = out_height)
  
  # 保存 PNG
  ggsave(filename = paste0(file_base, ".png"), 
         plot     = p_barplot, 
         width    = out_width, 
         height   = out_height, 
         dpi      = 300)
  
  print(paste("✅ GO barplot successfully generated and saved for:", target_cluster))
}

# Source: source/project/Reference/combine/scripts/MultiSpecies/scripts/finallmarkers.r
# Monkey findallmarkers----------------
monkey_wm_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/monkey_hip_mg_v2.qs"))
DefaultAssay(monkey_wm_mg) <- "RNA"
Idents(monkey_wm_mg) <- "new_celltype"
monkey_mg_markers <- FindAllMarkers(monkey_wm_mg, 
                                   only.pos = TRUE, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.25)
qsave(monkey_mg_markers,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/monkey_mg_markers.qs"))

# S4D  Mouse_MG0_like top10 GO Biological Process barplot
# Source: source/project/Reference/combine/scripts/MultiSpecies/scripts/Mouse_M0_like_GO.r
library(qs)
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(stringr)
library(ggplot2)

# MODULE : GO Enrichment and Barplot for a Single Cluster -------------------

## 1 Parameters --------------------------------------------------------------
# -- 数据输入参数 --
input_filepath  <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/mouse_mg_markers.qs")
target_cluster  <- "Mouse_MG0_like"       # 目标提取的 cluster 名称

# -- 过滤与分析阈值 --
p_val_adj_cut   <- 0.05                   # 显著性阈值
log2FC_cut      <- 0.5                    # Log2FC 阈值
db_org          <- "org.Mm.eg.db"         # 物种数据库: "org.Mm.eg.db" (鼠) / "org.Hs.eg.db" (人)
ontology_type   <- "BP"                   # GO 类型: "BP", "CC", "MF" 或 "ALL"
top_n_go        <- 10                     # 图中展示的 Top GO 通路数量

# -- 绘图与美化参数 --
cluster_color   <- "#33A02C"              # 柱状图的填充颜色
wrap_width      <- 60                     # 文本自动换行的字符数限制

# -- 输出参数 --
out_dir         <- "results/figure/cluster_go/"  # 图片保存的文件夹路径
out_width       <- 8                      # 输出图片宽度
out_height      <- 5                      # 输出图片高度


## 2 Process Data ------------------------------------------------------------
# 动态加载物种数据库
require(db_org, character.only = TRUE)

# 1. 读取 Marker 数据框 (修复：使用 qread 真正载入数据)
marker_df <- qread(input_filepath)

# 2. 过滤目标群的 Marker
cluster_markers <- marker_df %>%
  dplyr::filter(cluster == target_cluster) %>%
  dplyr::filter(p_val_adj < p_val_adj_cut) %>% 
  dplyr::filter(avg_log2FC > log2FC_cut) %>%
  dplyr::arrange(desc(avg_log2FC))

# 3. ID 转换 (SYMBOL -> ENTREZID)
gene_conv <- bitr(cluster_markers$gene, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = db_org)

# 4. 运行 GO 富集分析
ego_cluster <- enrichGO(
  gene          = gene_conv$ENTREZID,
  OrgDb         = db_org,
  ont           = ontology_type,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# 5. 提取富集结果并准备作图数据
plot_data <- as.data.frame(ego_cluster)
dir.create(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/table"), recursive = TRUE)
write.xlsx(plot_data,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/table/mouse_go.xlsx"))

if(nrow(plot_data) > 0) {
  plot_data <- plot_data %>%
    dplyr::arrange(p.adjust) %>%
    head(top_n_go) %>%
    dplyr::mutate(
      logP = -log10(p.adjust),
      # 对过长的 Description 自动换行
      Description = stringr::str_wrap(Description, width = wrap_width) 
    )
  
  # 锁定因子顺序，确保显著性最高的排在最上面
  plot_data$Description <- factor(plot_data$Description, levels = rev(plot_data$Description))
} else {
  warning(paste("Warning: No significant GO terms found for", target_cluster))
}


## 3 Generate Plots ----------------------------------------------------------
# 初始化画图对象为空
p_barplot <- NULL

if(nrow(plot_data) > 0) {
  
  # 计算文字的 X 轴起始坐标 (放在坐标轴起点偏右一点点，确保整齐对齐)
  text_x_pos <- max(plot_data$logP) * 0.02 
  
  p_barplot <- ggplot(plot_data, aes(x = logP, y = Description)) +
    geom_col(width = 0.8, alpha = 0.8, fill = cluster_color, color = "black", linewidth = 0.3) +
    # 文字标签：固定在图表左侧对齐
    geom_text(aes(x = text_x_pos, label = Description), 
              size = 4.5, hjust = 0, color = "black") + 
    theme_minimal() +
    labs(
      title = paste0("Cluster ", target_cluster, " Top ", top_n_go, " GO Terms (", ontology_type, ")"), 
      x = "-log10(p.adjust)", 
      y = NULL
    ) +
    theme(
      axis.text.y  = element_blank(), # 隐藏原生 Y 轴文字
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(), # 隐藏网格线，画面更干净
      axis.line    = element_line(color = "black"),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    # 向右延展 X 轴空间，防止特别长的通路名字越界被裁切
    scale_x_continuous(expand = expansion(mult = c(0, 0.4))) 
}


## 4 Save Results ------------------------------------------------------------
if(!is.null(p_barplot)) {
  # 如果文件夹不存在，自动创建
  if(!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # 替换 cluster 名字里的非法字符（如果有），防止保存报错
  safe_name <- gsub("[^a-zA-Z0-9_]", "_", target_cluster) 
  file_base <- file.path(out_dir, paste0("Cluster_", safe_name, "_GO_", ontology_type))
  
  # 保存 PDF
  ggsave(filename = paste0(file_base, ".pdf"), 
         plot     = p_barplot, 
         width    = out_width, 
         height   = out_height)
  
  # 保存 PNG
  ggsave(filename = paste0(file_base, ".png"), 
         plot     = p_barplot, 
         width    = out_width, 
         height   = out_height, 
         dpi      = 300)
  
  print(paste("✅ GO barplot successfully generated and saved for:", target_cluster))
}

# Source: source/project/Reference/combine/scripts/MultiSpecies/scripts/finallmarkers.r
# Mouse findallmarkers----------------
mouse_wm_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/mouse_wm_mg_v2.qs"))
DefaultAssay(mouse_wm_mg) <- "RNA"
Idents(mouse_wm_mg) <- "new_celltype"
mouse_mg_markers <- FindAllMarkers(mouse_wm_mg, 
                                   only.pos = TRUE, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.25)
qsave(mouse_mg_markers,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/mouse_mg_markers.qs"))

# S4E  LPC GSE148676 stored UMAP with merged cluster0-cluster10 labels
# Source: source/external/LPC/scripts/34_gse148676_reference_spatial_module.R
# MODULE : Cell Type Annotation and Module Scoring -------------------
#
# Direct adaptation of spatial_module.r for merged GSE148676 cells.
# The original package calls and plotting structure are retained. Only the
# input/output paths and the GSE148676 cluster mapping differ. The saveRDS()
# directly after AddModuleScore makes the original reload-for-plot step work
# when this script is sourced from a fresh R session.

## 1 Parameters --------------------------------------------------------------
# Load dependencies
library(Seurat)
library(qs)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SCP)

# File paths
project <- Sys.getenv("MG_LPC_ROOT")
file_in_merfish <- file.path(
  project,
  "objects/simple_microglia/mouse_GSE148676_merged_groups_seurat_v5.rds"
)
file_in_markers <- file.path(
  project,
  "raw_data/reference/mouse_mg_markers.qs"
)
file_out_merfish <- file.path(
  project,
  "objects/simple_microglia/mouse_GSE148676_merged_groups_reference_module_score.rds"
)
file_out_umap_pdf <- file.path(
  project,
  "figures/simple_microglia/mouse_GSE148676_reference_module_score_UMAP.pdf"
)
file_out_umap_png <- sub("\\.pdf$", ".png", file_out_umap_pdf)
file_out_violin_pdf <- file.path(
  project,
  "figures/simple_microglia/mouse_GSE148676_reference_module_score_violin.pdf"
)
file_out_violin_png <- sub("\\.pdf$", ".png", file_out_violin_pdf)
dir.create(dirname(file_out_merfish), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(file_out_umap_pdf), recursive = TRUE, showWarnings = FALSE)

# Cell-type mapping dictionary (old merged-cluster label -> plot label)
col_old_celltype <- "merged_cluster"
col_new_celltype <- "plot_celltype"
celltype_dict <- c(
  "0" = "cluster0",
  "1" = "cluster1",
  "2" = "cluster2",
  "3" = "cluster3",
  "4" = "cluster4",
  "5" = "cluster5",
  "6" = "cluster6",
  "7" = "cluster7",
  "8" = "cluster8",
  "9" = "cluster9",
  "10" = "cluster10"
)

# Marker extraction and module-score parameters
marker_cluster_col <- "cluster"
target_cluster <- "cluster2"
marker_gene_col <- "gene"

score_module_name <- "myelin"
score_feature_name <- paste0(score_module_name, "1")
score_ctrl_size <- 10

# Plotting parameter retained from the reference script
plot_pt_size <- 0


## 2 Process Data ------------------------------------------------------------
# 2.1 Read data
merfish <- readRDS(file_in_merfish)
mouse_mg_markers <- qread(file_in_markers)

# 2.2 Replace cell-type names
merfish@meta.data[[col_new_celltype]] <- as.character(
  merfish@meta.data[[col_old_celltype]]
)
for (old_name in names(celltype_dict)) {
  new_name <- celltype_dict[[old_name]]
  merfish@meta.data[[col_new_celltype]][
    merfish@meta.data[[col_new_celltype]] %in% old_name
  ] <- new_name
}

# 2.3 Extract target-cluster genes and intersect with GSE148676 genes
filtered_data <- mouse_mg_markers %>%
  filter(.data[[marker_cluster_col]] == target_cluster)
int_gene <- intersect(filtered_data[[marker_gene_col]], rownames(merfish))

# 2.4 Calculate module score
merfish <- AddModuleScore(
  object = merfish,
  features = list(int_gene),
  ctrl = score_ctrl_size,
  name = score_module_name
)

# The reference SCP function accesses the Seurat v4 RNA@data slot. This
# conversion preserves the existing normalized RNA layer for plotting only.
merfish[["RNA"]] <- as(merfish[["RNA"]], "Assay")

# The reference script re-opens a previously saved scored object for plotting.
# Save here so the same operation is valid in a new session.
saveRDS(merfish, file = file_out_merfish)


## 3 Generate Plots ----------------------------------------------------------
merfish <- readRDS(file_out_merfish)

st_mg_colors_npg <- c(
  "cluster1" = "#E64B35",
  "cluster0" = "#4DBBD5",
  "cluster2" = "#00A087",
  "cluster3" = "#3C5488",
  "cluster4" = "#F39B7F",
  "cluster5" = "#7E6148",
  "cluster6" = "#8E6BBE",
  "cluster7" = "#7F7F7F",
  "cluster8" = "#66C2A5",
  "cluster9" = "#E6AB02",
  "cluster10" = "#A65628"
)

(p1 <- DimPlot(merfish, group.by = col_new_celltype, cols = st_mg_colors_npg))
## 4 Save Results ------------------------------------------------------------
# Save updated Seurat object (with the new cell-type and module score)
saveRDS(merfish, file = file_out_merfish)

# Save standalone UMAP and violin panels
ggsave(filename = file_out_umap_pdf, plot = p1, width = 5, height = 5)
ggsave(filename = file_out_umap_png, plot = p1, width = 5, height = 5, dpi = 300)

# Source: source/external/LPC/scripts/31_merge_gse148676_clusters_export.R
#!/usr/bin/env Rscript

# GSE148676 microglia: merge selected Seurat clusters without changing the
# original seurat_clusters column, then export group-level expression and UMAP.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(Seurat)
})

project <- Sys.getenv("MG_LPC_ROOT")
input_object <- file.path(
  project,
  "objects",
  "simple_microglia",
  "mouse_GSE148676_simple_microglia_seurat_v5.rds"
)
output_object <- file.path(
  project,
  "objects",
  "simple_microglia",
  "mouse_GSE148676_merged_groups_seurat_v5.rds"
)
report_dir <- file.path(project, "reports", "simple_microglia")
figure_dir <- file.path(project, "figures", "simple_microglia")

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

merged_3_9_10 <- c("3", "9", "10")
merged_0_1_2_4_5_8 <- c("0", "1", "2", "4", "5", "8")
all_expected_clusters <- as.character(0:17)
remaining_clusters <- setdiff(
  all_expected_clusters,
  c(merged_3_9_10, merged_0_1_2_4_5_8)
)
new_cluster_levels <- as.character(0:10)
old_to_new_cluster <- c(
  # Final IDs are ordered by merged-group cell count, from largest to smallest.
  setNames(rep("1", length(merged_3_9_10)), merged_3_9_10),
  setNames(rep("0", length(merged_0_1_2_4_5_8)), merged_0_1_2_4_5_8),
  setNames(as.character(2:10), remaining_clusters)
)

message("Reading ", input_object)
microglia <- readRDS(input_object)
DefaultAssay(microglia) <- "RNA"

original_cluster <- as.character(microglia$seurat_clusters)
unexpected_clusters <- setdiff(unique(original_cluster), all_expected_clusters)
if (length(unexpected_clusters) > 0) {
  stop(
    "Unexpected original cluster labels: ",
    paste(unexpected_clusters, collapse = ", ")
  )
}

merged_cluster_source <- original_cluster
merged_cluster_source[original_cluster %in% merged_3_9_10] <- "3_9_10"
merged_cluster_source[original_cluster %in% merged_0_1_2_4_5_8] <- "0_1_2_4_5_8"
microglia$merged_cluster <- factor(
  unname(old_to_new_cluster[original_cluster]),
  levels = new_cluster_levels
)
microglia$merged_cluster_source <- merged_cluster_source
Idents(microglia) <- "merged_cluster"

# S4F  LPC myelin1 signature score by merged clusters; cluster1 compared to each other cluster
# Source: source/external/LPC/scripts/34_gse148676_reference_spatial_module.R
# MODULE : Cell Type Annotation and Module Scoring -------------------
#
# Direct adaptation of spatial_module.r for merged GSE148676 cells.
# The original package calls and plotting structure are retained. Only the
# input/output paths and the GSE148676 cluster mapping differ. The saveRDS()
# directly after AddModuleScore makes the original reload-for-plot step work
# when this script is sourced from a fresh R session.

## 1 Parameters --------------------------------------------------------------
# Load dependencies
library(Seurat)
library(qs)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SCP)

# File paths
project <- Sys.getenv("MG_LPC_ROOT")
file_in_merfish <- file.path(
  project,
  "objects/simple_microglia/mouse_GSE148676_merged_groups_seurat_v5.rds"
)
file_in_markers <- file.path(
  project,
  "raw_data/reference/mouse_mg_markers.qs"
)
file_out_merfish <- file.path(
  project,
  "objects/simple_microglia/mouse_GSE148676_merged_groups_reference_module_score.rds"
)
file_out_umap_pdf <- file.path(
  project,
  "figures/simple_microglia/mouse_GSE148676_reference_module_score_UMAP.pdf"
)
file_out_umap_png <- sub("\\.pdf$", ".png", file_out_umap_pdf)
file_out_violin_pdf <- file.path(
  project,
  "figures/simple_microglia/mouse_GSE148676_reference_module_score_violin.pdf"
)
file_out_violin_png <- sub("\\.pdf$", ".png", file_out_violin_pdf)
dir.create(dirname(file_out_merfish), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(file_out_umap_pdf), recursive = TRUE, showWarnings = FALSE)

# Cell-type mapping dictionary (old merged-cluster label -> plot label)
col_old_celltype <- "merged_cluster"
col_new_celltype <- "plot_celltype"
celltype_dict <- c(
  "0" = "cluster0",
  "1" = "cluster1",
  "2" = "cluster2",
  "3" = "cluster3",
  "4" = "cluster4",
  "5" = "cluster5",
  "6" = "cluster6",
  "7" = "cluster7",
  "8" = "cluster8",
  "9" = "cluster9",
  "10" = "cluster10"
)

# Marker extraction and module-score parameters
marker_cluster_col <- "cluster"
target_cluster <- "cluster2"
marker_gene_col <- "gene"

score_module_name <- "myelin"
score_feature_name <- paste0(score_module_name, "1")
score_ctrl_size <- 10

# Plotting parameter retained from the reference script
plot_pt_size <- 0


## 2 Process Data ------------------------------------------------------------
# 2.1 Read data
merfish <- readRDS(file_in_merfish)
mouse_mg_markers <- qread(file_in_markers)

# 2.2 Replace cell-type names
merfish@meta.data[[col_new_celltype]] <- as.character(
  merfish@meta.data[[col_old_celltype]]
)
for (old_name in names(celltype_dict)) {
  new_name <- celltype_dict[[old_name]]
  merfish@meta.data[[col_new_celltype]][
    merfish@meta.data[[col_new_celltype]] %in% old_name
  ] <- new_name
}

# 2.3 Extract target-cluster genes and intersect with GSE148676 genes
filtered_data <- mouse_mg_markers %>%
  filter(.data[[marker_cluster_col]] == target_cluster)
int_gene <- intersect(filtered_data[[marker_gene_col]], rownames(merfish))

# 2.4 Calculate module score
merfish <- AddModuleScore(
  object = merfish,
  features = list(int_gene),
  ctrl = score_ctrl_size,
  name = score_module_name
)

# The reference SCP function accesses the Seurat v4 RNA@data slot. This
# conversion preserves the existing normalized RNA layer for plotting only.
merfish[["RNA"]] <- as(merfish[["RNA"]], "Assay")

# The reference script re-opens a previously saved scored object for plotting.
# Save here so the same operation is valid in a new session.
saveRDS(merfish, file = file_out_merfish)


## 3 Generate Plots ----------------------------------------------------------
merfish <- readRDS(file_out_merfish)

st_mg_colors_npg <- c(
  "cluster1" = "#E64B35",
  "cluster0" = "#4DBBD5",
  "cluster2" = "#00A087",
  "cluster3" = "#3C5488",
  "cluster4" = "#F39B7F",
  "cluster5" = "#7E6148",
  "cluster6" = "#8E6BBE",
  "cluster7" = "#7F7F7F",
  "cluster8" = "#66C2A5",
  "cluster9" = "#E6AB02",
  "cluster10" = "#A65628"
)
# The first group is the merged candidate group (old clusters 3, 9, and 10),
# preserving the reference comparison pattern of first group versus all others.
target_order <- c(
  "cluster1", "cluster0", "cluster2", "cluster3", "cluster4", "cluster5",
  "cluster6", "cluster7", "cluster8", "cluster9", "cluster10"
)
merfish$plot_celltype <- factor(merfish$plot_celltype, levels = target_order)

celltype <- levels(merfish$plot_celltype)
my_comparisons <- lapply(celltype[-1], function(x) c(x, celltype[1]))

p2 <- FeatureStatPlot(
  srt = merfish, group.by = "plot_celltype",
  stat.by = c(score_feature_name), add_box = FALSE,
  comparisons = my_comparisons,
  bg.by = "plot_celltype",
  palcolor = st_mg_colors_npg,
  bg_palcolor = st_mg_colors_npg
)
ggsave(filename = file_out_violin_pdf, plot = p2, width = 7, height = 5)
ggsave(filename = file_out_violin_png, plot = p2, width = 7, height = 5, dpi = 300)

# S4G  LPC Trim2/Rnf220/Mbp/Mog/Olig2/St18 normalized-expression UMAP panels
# Source: source/external/LPC/scripts/28_gse148676_extended_marker_featureplots.R
project <- Sys.getenv("MG_LPC_ROOT")
object_path <- file.path(
  project,
  "objects",
  "simple_microglia",
  "mouse_GSE148676_simple_microglia_seurat_v5.rds"
)
figure_dir <- file.path(project, "figures", "simple_microglia")
individual_figure_dir <- file.path(
  figure_dir,
  "mouse_GSE148676_individual_featureplots"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_figure_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading ", object_path)
microglia <- readRDS(object_path)
DefaultAssay(microglia) <- "RNA"

requested_genes <- c("Trim2", "Rnf220", "St18", "Olig2", "Mbp", "Mog")
available_genes <- rownames(microglia)
matched_genes <- vapply(
  requested_genes,
  function(gene) {
    hits <- available_genes[toupper(available_genes) == toupper(gene)]
    if (length(hits) == 0) {
      stop("Gene is absent from the RNA assay: ", gene)
    }
    hits[[1]]
  },
  character(1)
)
names(matched_genes) <- requested_genes
message(
  "Matched features: ",
  paste(names(matched_genes), matched_genes, sep = "=", collapse = ", ")
)

umap_frame <- as.data.frame(Embeddings(microglia, "umap"))
colnames(umap_frame)[1:2] <- c("UMAP_1", "UMAP_2")
cell_order <- rownames(umap_frame)

normalized_matrix <- LayerData(
  microglia[["RNA"]],
  layer = "data"
)[matched_genes, cell_order, drop = FALSE]
for (gene_label in requested_genes) {
  feature_name <- matched_genes[[gene_label]]
  expression <- as.numeric(normalized_matrix[feature_name, ])
  expression_values[[gene_label]] <- expression

  plot_frame <- umap_frame
  plot_frame$expression <- expression
  positive_frame <- plot_frame[plot_frame$expression > 0, , drop = FALSE]
  positive_frame <- positive_frame[
    order(positive_frame$expression),
    ,
    drop = FALSE
  ]

  if (nrow(positive_frame) == 0) {
    upper_limit <- 1
  } else {
    upper_limit <- as.numeric(
      quantile(
        positive_frame$expression,
        0.99,
        names = FALSE
      )
    )
  }

  gene_plots[[gene_label]] <- ggplot(
    plot_frame,
    aes(x = UMAP_1, y = UMAP_2)
  ) +
    geom_point(
      color = "#B8B8B8",
      size = point_size,
      alpha = 0.60,
      shape = 16,
      stroke = 0
    ) +
    geom_point(
      data = positive_frame,
      aes(color = expression),
      size = point_size * 1.55,
      alpha = 1,
      shape = 16,
      stroke = 0
    ) +
    scale_color_gradient(
      low = "#FFD166",
      high = "#B2182B",
      limits = c(0, upper_limit),
      oob = scales::squish
    ) +
    labs(
      title = gene_label,
      color = "Log-normalized\nexpression",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    coord_equal() +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
}
# Export every marker panel as a standalone figure. These files retain the
# same UMAP coordinates and normalized RNA data used in the combined panel.
for (gene_label in requested_genes) {
  individual_path <- file.path(
    individual_figure_dir,
    paste0(
      "mouse_GSE148676_",
      gene_label,
      "_log_normalized_FeaturePlot.png"
    )
  )
  ggsave(
    individual_path,
    gene_plots[[gene_label]],
    width = 7.2,
    height = 6.4,
    dpi = 300,
    bg = "white"
  )
  message("Wrote ", individual_path)

  individual_pdf_path <- sub("\\.png$", ".pdf", individual_path)
  ggsave(
    individual_pdf_path,
    gene_plots[[gene_label]],
    width = 7.2,
    height = 6.4,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  message("Wrote ", individual_pdf_path)
}

# S4H  100% retained-MG composition, LPC OLM and Monkey/Mouse MG0-like
# Source: source/project/Reference/combine/scripts/MultiSpecies/results/MG0_like_proportion_LPC_GSE148676_Monkey_Mouse_20260921/100pct_composition_and_umap_source_audit/plot_LPC_MG0_Monkey_Mouse_100pct_composition_umap_audit_20260921.R
lpc_metadata_path <- file.path(
  lpc_project_root,
  "reports/simple_microglia/mouse_GSE148676_merged_groups_cell_metadata.csv"
)
multi_species_object_path <- file.path(
  project_root,
  "combine/scripts/MultiSpecies/data/processed_data/02.cross_species_integrated_noref_Human_Mouse_Monkey.qs"
)
if (!file.exists(lpc_metadata_path)) stop("LPC metadata/UMAP export is absent: ", lpc_metadata_path)
if (!file.exists(multi_species_object_path)) stop("Multi-species object is absent: ", multi_species_object_path)

# LPC (GSE148676): use the existing per-cell export from the merged-cluster
# microglia object. Its UMAP is the stored UMAP and cluster 1 is the project's
# designated M0/OLM group.
lpc_metadata <- read.csv(lpc_metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
required_lpc_columns <- c("cell", "merged_cluster", "UMAP_1", "UMAP_2")
if (!all(required_lpc_columns %in% colnames(lpc_metadata))) {
  stop("LPC metadata is missing: ", paste(setdiff(required_lpc_columns, colnames(lpc_metadata)), collapse = ", "))
}
if (anyDuplicated(lpc_metadata$cell) > 0L) stop("LPC metadata contains duplicated cell IDs.")
lpc_frame <- data.frame(
  Dataset = "LPC (GSE148676, mouse)",
  UMAP_1 = as.numeric(lpc_metadata$UMAP_1),
  UMAP_2 = as.numeric(lpc_metadata$UMAP_2),
  Target_group = ifelse(as.character(lpc_metadata$merged_cluster) == "1", "OLM", "Other MG"),
  stringsAsFactors = FALSE
)
lpc_total <- nrow(lpc_frame)
lpc_target <- sum(lpc_frame$Target_group == "OLM")
lpc_cluster_values <- paste(sort(unique(as.character(lpc_metadata$merged_cluster))), collapse = ", ")
if (lpc_total == 0L || lpc_target == 0L) stop("LPC has no cells and/or no merged cluster 1 cells.")
rm(lpc_metadata)
gc(verbose = FALSE)

# Monkey and Mouse: retain the stored integrated UMAP and use the original
# author-defined integrated cluster 1 as MG0-like for each species.
multi_species <- qread(multi_species_object_path)
multi_metadata <- multi_species[[]]
required_multi_columns <- c("Species", "integrated_snn_res.0.1")
if (!all(required_multi_columns %in% colnames(multi_metadata))) {
  stop("Multi-species metadata is missing: ", paste(setdiff(required_multi_columns, colnames(multi_metadata)), collapse = ", "))
}
if (!"umap" %in% Reductions(multi_species)) stop("The multi-species object has no stored UMAP reduction.")
multi_umap <- Embeddings(multi_species, reduction = "umap")
if (!identical(rownames(multi_metadata), rownames(multi_umap))) {
  stop("Multi-species metadata and stored UMAP have non-identical cell order.")
}

make_multispecies_frame <- function(species_name) {
  keep <- as.character(multi_metadata$Species) == species_name
  frame <- data.frame(
    Dataset = species_name,
    UMAP_1 = as.numeric(multi_umap[keep, 1]),
    UMAP_2 = as.numeric(multi_umap[keep, 2]),
    Target_group = ifelse(
      as.character(multi_metadata$integrated_snn_res.0.1[keep]) == "1",
      "MG0-like", "Other MG"
    ),
    stringsAsFactors = FALSE
  )
  if (nrow(frame) == 0L || sum(frame$Target_group == "MG0-like") == 0L) {
    stop(species_name, " has no cells and/or no integrated cluster 1 cells.")
  }
  frame
}
monkey_frame <- make_multispecies_frame("Monkey")
mouse_frame <- make_multispecies_frame("Mouse")
monkey_cluster_values <- paste(sort(unique(as.character(
  multi_metadata$integrated_snn_res.0.1[as.character(multi_metadata$Species) == "Monkey"]
))), collapse = ", ")
mouse_cluster_values <- paste(sort(unique(as.character(
  multi_metadata$integrated_snn_res.0.1[as.character(multi_metadata$Species) == "Mouse"]
))), collapse = ", ")
composition_table <- data.frame(
  Dataset = c("LPC (GSE148676, mouse)", "Monkey", "Mouse"),
  Target_label = c("OLM", "MG0-like", "MG0-like"),
  Target_definition = c(
    "merged_cluster == 1",
    "integrated_snn_res.0.1 == 1",
    "integrated_snn_res.0.1 == 1"
  ),
  Target_cells = c(
    lpc_target,
    sum(monkey_frame$Target_group == "MG0-like"),
    sum(mouse_frame$Target_group == "MG0-like")
  ),
  Total_microglia_cells = c(nrow(lpc_frame), nrow(monkey_frame), nrow(mouse_frame)),
  stringsAsFactors = FALSE
)
composition_table$Target_percent <- 100 * composition_table$Target_cells / composition_table$Total_microglia_cells
composition_table$Other_MG_cells <- composition_table$Total_microglia_cells - composition_table$Target_cells
composition_table$Other_MG_percent <- 100 - composition_table$Target_percent
composition_table$Dataset <- factor(
  composition_table$Dataset,
  # coord_flip displays the last factor level at the top. This ordering gives
  # the requested LPC, Monkey, Mouse top-to-bottom reading sequence.
  levels = c("Mouse", "Monkey", "LPC (GSE148676, mouse)")
)
write.csv(composition_table, file.path(result_dir, "MG0_like_100pct_composition_counts.csv"), row.names = FALSE)
composition_plot <- ggplot(composition_long, aes(x = Dataset, y = Percent, fill = Segment)) +
  # The reference composition plot places the target subtype at the right end
  # of the 100% bar, so Other MG is stacked first and the target follows it.
  geom_col(
    width = 0.72,
    color = "#1A1A1A",
    linewidth = 0.35,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    data = composition_table,
    aes(x = Dataset, y = Other_MG_percent / 2, label = sprintf("Other MG %.1f%%", Other_MG_percent)),
    inherit.aes = FALSE,
    size = 4.0,
    color = "#303030"
  ) +
  geom_text(
    data = composition_table,
    aes(x = Dataset, y = 103.0, label = target_annotation),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.9,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("Other MG" = "#D0D0D0", "OLM" = "#4DAF4A", "MG0-like" = "#8FCB6A"),
    breaks = c("OLM", "MG0-like", "Other MG"),
    name = "Microglial group"
  ) +
  scale_y_continuous(
    limits = c(0, 115),
    breaks = seq(0, 100, by = 25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_flip(clip = "off") +
  labs(
    title = "MG0/MG0-like composition among all microglia",
    subtitle = "Every bar is normalized to 100% of the retained microglial cells in its input dataset",
    x = NULL,
    y = "Proportion"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.y = element_text(face = "bold", size = 11),
    axis.title.x = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    plot.margin = margin(t = 10, r = 165, b = 8, l = 8, unit = "pt")
  )
ggsave(
  file.path(result_dir, "MG0_like_100pct_composition.pdf"),
  composition_plot,
  width = 10.8, height = 4.9, device = grDevices::cairo_pdf, bg = "white"
)
ggsave(
  file.path(result_dir, "MG0_like_100pct_composition.png"),
  composition_plot,
  width = 10.8, height = 4.9, dpi = 300, bg = "white"
)


# Figure 3 ------------------------------------------------------------

# F3A  Human/monkey/mouse integrated microglial UMAP, species-colored UMAP and OLM-like highlights; decorative species icons and manually drawn red circles
# Source: source/project/code_upload/MultiSpecies.r
manual_colors <- c("MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A", 
                   "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628")

custom_colors <- c("MG0" = "#4DAF4A", "Mouse_MG0_like" = "#33A02C",  
                   "Monkey_MG0_like" = "#B2DF8A", "Other_MG" = "#D9D9D9")

species_colors <- c("Human" = "#E69F00", "Mouse" = "#56B4E9", "Monkey" = "#009E73")

integrated_mg <- qread('data/processed_data/02.cross_species_integrated_noref_Human_Mouse_Monkey.qs')

integrated_mg$new_celltype <- as.character(integrated_mg$CellType1)

target_cells_mouse <- integrated_mg$integrated_snn_res.0.1 == '1' & integrated_mg$Species == 'Mouse'
integrated_mg$new_celltype[target_cells_mouse] <- "Mouse_MG0_like"

target_cells_monkey <- integrated_mg$integrated_snn_res.0.1 == '1' & integrated_mg$Species == 'Monkey'
integrated_mg$new_celltype[target_cells_monkey] <- "Monkey_MG0_like"

keep_types <- c("MG0", "Mouse_MG0_like", "Monkey_MG0_like")
integrated_mg$new_celltype <- ifelse(
  integrated_mg$new_celltype %in% keep_types,
  integrated_mg$new_celltype, 
  "Other_MG"
)
Idents(integrated_mg) <- "new_celltype"

# Fig.01 ---------------
p1 <- DimPlot(integrated_mg, 
              reduction = "umap", 
              group.by = "new_celltype", 
              cols = custom_colors,
              split.by = 'Species',
              pt.size = 0.3,       
              raster = FALSE) + 
  ggtitle("MG0 Cross-Species Subsets") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoAxes()

ggsave('results/figure/species_dimplot.pdf', p1, width = 9, height = 3)
ggsave('results/figure/species_dimplot.png', p1, width = 9, height = 3)

# Fig.02 ---------------
p2 <- DimPlot(integrated_mg, 
              reduction = "umap", 
              group.by = "Species",     
              cols = species_colors,    
              pt.size = 0.2,            
              raster = FALSE) + 
  ggtitle("Species") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoAxes()                             

ggsave('results/figure/species_celltype.pdf', p2, width = 4, height = 3)
ggsave('results/figure/species_celltype.png', p2, width = 4, height = 3)

# F3B  Monkey OLM-like versus other MG expression of MOG, MBP, RNF220, TRIM2
# Source: source/project/code_upload/MultiSpecies.r
monkey_hip_mg <- qread('data/processed_data/monkey_hip_mg.qs')
monkey_integrated_barcodes <- colnames(integrated_mg)[integrated_mg$Species == "Monkey"]
monkey_celltypes <- integrated_mg$new_celltype[integrated_mg$Species == "Monkey"]
orig_monkey_barcodes <- gsub("^Monkey_", "", monkey_integrated_barcodes)

monkey_anno_df <- data.frame(new_celltype = monkey_celltypes, row.names = orig_monkey_barcodes)
monkey_hip_mg <- AddMetaData(monkey_hip_mg, metadata = monkey_anno_df)
qsave(monkey_hip_mg, 'data/processed_data/monkey_hip_mg_v2.qs')
# Fig.04 ---------------
DefaultAssay(monkey_hip_mg) <- 'RNA'
monkey_hip_mg$new_celltype <- factor(monkey_hip_mg$new_celltype[, 1], levels = c("Monkey_MG0_like", "Other_MG"))

compare_groups_monkey <- setdiff(levels(monkey_hip_mg$new_celltype), "Monkey_MG0_like")
my_comparisons_monkey <- lapply(compare_groups_monkey, function(x) c(x, "Monkey_MG0_like"))

p_violin_monkey <- FeatureStatPlot(
  srt         = monkey_hip_mg, 
  group.by    = "new_celltype",
  stat.by     = c("MOG", "MBP", "RNF220", "TRIM2"), 
  add_box     = FALSE, 
  comparisons = my_comparisons_monkey,
  bg.by       = "new_celltype",
  palcolor    = custom_colors,
  bg_palcolor = custom_colors,
  ncol        = 4
)

ggsave(filename = 'results/figure/Monkey_Violin_Stat_Markers.pdf', plot = p_violin_monkey, width = 16, height = 4)

# F3C  Mouse OLM-like versus other MG expression of Mog, Mbp, Rnf220, Trim2
# Source: source/project/code_upload/MultiSpecies.r
mouse_wm_mg <- qread('data/processed_data/mouse_wm_mg.qs')
mouse_integrated_barcodes <- colnames(integrated_mg)[integrated_mg$Species == "Mouse"]
mouse_celltypes <- integrated_mg$new_celltype[integrated_mg$Species == "Mouse"]
orig_mouse_barcodes <- gsub("^Mouse_", "", mouse_integrated_barcodes)

mouse_anno_df <- data.frame(new_celltype = mouse_celltypes, row.names = orig_mouse_barcodes)
mouse_wm_mg <- AddMetaData(mouse_wm_mg, metadata = mouse_anno_df)
qsave(mouse_wm_mg, 'data/processed_data/mouse_wm_mg_v2.qs')

# Fig.03 ---------------
DefaultAssay(mouse_wm_mg) <- 'RNA'
mouse_wm_mg$new_celltype <- factor(mouse_wm_mg$new_celltype[, 1], levels = c("Mouse_MG0_like", "Other_MG"))

compare_groups_mouse <- setdiff(levels(mouse_wm_mg$new_celltype), "Mouse_MG0_like")
my_comparisons_mouse <- lapply(compare_groups_mouse, function(x) c(x, "Mouse_MG0_like"))

p_violin_mouse <- FeatureStatPlot(
  srt         = mouse_wm_mg, 
  group.by    = "new_celltype",
  stat.by     = c("Mog", "Mbp", "Rnf220", "Trim2"), 
  add_box     = FALSE, 
  comparisons = my_comparisons_mouse,
  bg.by       = "new_celltype",
  palcolor    = custom_colors,
  bg_palcolor = custom_colors,
  ncol        = 4
)

ggsave(filename = 'results/figure/Mouse_Violin_Stat_Markers.pdf', plot = p_violin_mouse, width = 16, height = 4)

# F3D  Monkey-3M Iba1/Trim2/DAPI immunofluorescence, 10um

# F3E  C57-CC-LPS mouse corpus callosum Iba1/Trim2/DAPI immunofluorescence, 10um

# F3F  C57-CC-LPS mouse MBP/Iba1/DAPI immunofluorescence, 2um

# F3G  C57-CC-LPS mouse RNF220/Iba1/DAPI immunofluorescence, 10um

# F3H  C57-LPS dual immunogold TEM, Trim2 10nm and Iba1 4nm, overview and zooms


# Supplementary Figure 5 ------------------------------------------------------------

# S5A  Iba1/Trim2/DAPI immunofluorescence, needle-like cells and manually marked arrows,20um

# S5B  Iba1/Trim2/DAPI overview and zoom with 500um/200um scalebars

# S5C  OLM proportion in five observations
# Source: source/recovered/experimental/S3C_OLM_proportion.R
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[[1L]] else file.path(Sys.getenv("MG_PLOT_INPUT_ROOT", "data/plot_inputs"), "2026-7-29-统计图.xlsx")
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "S3C")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
raw_v <- read_excel(input_file, sheet = "Sheet2", range = "A1:F11", col_names = FALSE)
value <- suppressWarnings(as.numeric(unlist(raw_v[, 6])))
olm <- data.frame(group = "OLM", value = value[!is.na(value)])
stopifnot(nrow(olm) == 5L)
p_olm <- ggplot(olm, aes(x = group, y = value)) +
  geom_boxplot(
    width = 0.38, outlier.shape = NA,
    fill = "#1786C0", alpha = 0.35,
    color = "#1786C0", linewidth = 0.9,
    staplewidth = 0.55
  ) +
  geom_point(
    shape = 21, size = 2.6, stroke = 0.9,
    fill = "#1786C0", color = "black",
    position = position_jitter(width = 0.075, height = 0, seed = 20260729)
  ) +
  scale_y_continuous(
    name = expression(Trim2^"+" * Iba1^"+" * OLM/Iba1^"+" ~ MG),
    breaks = seq(0, 0.12, 0.02),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL) +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(linewidth = 0.9),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black"),
    plot.margin = margin(8, 8, 7, 7)
  ) +
  coord_cartesian(ylim = c(0, 0.12))
ggsave(file.path(output_dir, "S3C_OLM_proportion.png"), p_olm, width = 3.0, height = 4.3, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "S3C_OLM_proportion.pdf"), p_olm, width = 3.0, height = 4.3, device = cairo_pdf)
write.csv(olm, file.path(output_dir, "S3C_plot_values.csv"), row.names = FALSE)
cat("S3C values:\n")
print(olm, digits = 17)
cat("\nSession:\n")
print(sessionInfo())

# S5D  DLB CellChat circle, major Chia source annotations and MG0-MG5
# Source: source/project/Reference/combine/scripts/CellChat/scripts/03.con_DLB_cellchat.r
# CellChat流程

# 全局设置----------------
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/CellChat"))
rm(list = ls())

# library------------------
library(CellChat)
library(Seurat)
library(qs)
library(ggplot2)
library(patchwork)
library(igraph)
library(ComplexHeatmap)

# 加载数据---------------------
seu <- qread('data/processed_data/Chia_DLB.qs')

unique(seu$Sample)
# [1] DLB_2 DLB_4 DLB_5 DLB_6 DLB_7 DLB_8
# Levels: ADD_4 ADD_5 ADD_6 ADD_7 CTRL_4 CTRL_5 CTRL_6 CTRL_7 DLB_2 DLB_4 DLB_5 DLB_6 DLB_7 DLB_8 PDD_1 PDD_2 PDD_4 PDD_5 PDD_7 PDD_9
# seu <- subset(seu,subset = Sample %in% c('DLB_2'))

seu.data.input <- GetAssayData(seu, assay = "SCT", slot = "data")
seu.identity <- subset(seu@meta.data, select = "cellchat_type")
cellchat <- createCellChat(object = seu.data.input, meta = seu.identity, group.by = "cellchat_type")

cellchat <- addMeta(cellchat, meta = seu.identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") 
groupSize <- as.numeric(table(cellchat@idents))

# 设置配体-受体相互作用数据库
CellChatDB <- CellChatDB.human  
CellChatDB.use <- subsetDB(CellChatDB, search = c("Cell-Cell Contact"))
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = FALSE
                              # ,population.size = T
                              )
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# weight<-as.data.frame(cellchat@net$weight)
# write.table(weight,file="results/table/weight_con_DLB.xls",sep="\t")
# count<-as.data.frame(cellchat@net$count)
# write.table(count,file="results/table/count_con_DLB.xls",sep="\t")
# netP<-as.data.frame(cellchat@netP[["prob"]])
# write.table(netP,file="results/table/netP_con_DLB.xls",sep="\t")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# 保存-----------
qsave(cellchat,'data/processed_data/con_DLB_cellchat.qs')

cellchat <- qread('data/processed_data/con_DLB_cellchat.qs')
# netVisual_circle(cellchat@net$weight, vertex.weight = 13, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cell_colors <- c(
  "MG0" = "#4DAF4A", 
  "MG1" = "#A65628",
  "MG2" = "#377EB8", 
  "MG3" = "#984EA3", 
  "MG4" = "#E6AB02", 
  "MG5" = "#8B0000", 
  "Excitatory neuron"         = "#E31A1C",
  "Inhibitory neuron"         = "#E32636",
  "Astrocyte"      = "#FF7F00",
  "Oligodendrocyte"           = "#1F78B4",
  "OPC"            = "#A6CEE3",
  "Ependymal cell"      = "#6A3D9A",
  "Endothelial cell"       = "#FB9A99"
)


# 1. 获取当前通讯矩阵中所有细胞类型的准确名称和顺序
group.names <- rownames(cellchat@net$weight)

# 2. 根据该顺序，从你的颜色字典中提取并重新排列颜色
# 确保所有 group.names 都在你的 cell_colors 字典中有定义，否则会产生 NA
ordered_colors <- cell_colors[group.names]

pdf('results/figure/DLB_cellchat.pdf',width = 5,height = 5)
netVisual_circle(
  cellchat@net$weight, 
  vertex.weight = 13, 
  weight.scale = T, 
  label.edge = F, 
  title.name = "Interaction weights/strength",
  sources.use = c("MG0","MG1","MG2","MG3","MG4","MG5"),
  color.use = ordered_colors  # 注意这里使用的是重新对齐顺序的 ordered_colors
)
dev.off()

# S5E  DLB CellChat selected ligand-receptor bubbles, OLM-major cells/MG-OLIG/OLIG-MG
# Source: source/project/Reference/combine/scripts/CellChat/scripts/multi_DLB_cellchat_pop.r
# 全局设置----------------
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/CellChat"))
rm(list = ls())

# library------------------
library(CellChat)
library(Seurat)
library(qs)
library(ggplot2)
library(patchwork)
library(igraph)
library(ComplexHeatmap)

# 加载全局数据---------------------
seu <- qread('data/processed_data/Chia_DLB.qs')

# 获取所有以 "DLB" 开头的样本名
all_samples <- unique(seu$Sample)
dlb_samples <- grep("^DLB", all_samples, value = TRUE)
print(paste("找到以下 DLB 样本将进行分析:", paste(dlb_samples, collapse = ", ")))

# 预设颜色字典 (移到循环外部，避免重复定义)
cell_colors <- c(
  "MG0" = "#4DAF4A", 
  "MG1" = "#A65628",
  "MG2" = "#377EB8", 
  "MG3" = "#984EA3", 
  "MG4" = "#E6AB02", 
  "MG5" = "#8B0000", 
  "Excitatory neuron"         = "#E31A1C",
  "Inhibitory neuron"         = "#E32636",
  "Astrocyte"      = "#FF7F00",
  "Oligodendrocyte"            = "#1F78B4",
  "OPC"             = "#A6CEE3",
  "Ependymal cell"      = "#6A3D9A",
  "Endothelial cell"       = "#FB9A99"
)

# 提前设置好受体-配体相互作用数据库
CellChatDB <- CellChatDB.human  
CellChatDB.use <- subsetDB(CellChatDB, search = c("Cell-Cell Contact"))

# 开始循环---------------------
for (sample_id in dlb_samples) {
  
  cat("\n======================================================\n")
  cat("正在处理样本:", sample_id, "\n")
  cat("======================================================\n")
  
  # 使用 tryCatch 捕获潜在错误，防止循环中断
  tryCatch({
    # 1. 提取当前样本的数据
    seu_sub <- subset(seu, subset = Sample == sample_id)
    
    seu.data.input <- GetAssayData(seu_sub, assay = "SCT", slot = "data")
    seu.identity <- subset(seu_sub@meta.data, select = "cellchat_type")
    
    # 2. 创建并初始化 CellChat 对象
    cellchat <- createCellChat(object = seu.data.input, meta = seu.identity, group.by = "cellchat_type")
    cellchat <- addMeta(cellchat, meta = seu.identity, meta.name = "labels")
    cellchat <- setIdent(cellchat, ident.use = "labels") 
    
    cellchat@DB <- CellChatDB.use
    
    # 3. 核心计算流程
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- smoothData(cellchat, adj = PPI.human)
    
    cellchat <- computeCommunProb(cellchat, raw.use = FALSE,population.size = T)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    
    # 4. 保存 CellChat QS 文件 (动态命名)
    qs_file_path <- paste0('data/processed_data/con_', sample_id, '_pop_cellchat.qs')
    qsave(cellchat, qs_file_path)
    cat("已保存 QS 数据至:", qs_file_path, "\n")
    
    # 5. 获取当前通讯矩阵中的细胞名称和对应的颜色
    group.names <- rownames(cellchat@net$weight)
    ordered_colors <- cell_colors[group.names]
    
    # 6. 绘制并保存 Circle Plot (动态命名)
    circle_pdf_path <- paste0('results/figure/', sample_id, '_pop_cellchat_circle.pdf')
    pdf(circle_pdf_path, width = 5, height = 5)
    netVisual_circle(
      cellchat@net$weight, 
      vertex.weight = 13, 
      weight.scale = T, 
      label.edge = F, 
      title.name = paste("Interaction weights/strength -", sample_id),
      sources.use = c("MG0","MG1","MG2","MG3","MG4","MG5"),
      color.use = ordered_colors 
    )
    dev.off()
    cat("已保存 Circle 图至:", circle_pdf_path, "\n")
    
    # 7. 绘制并保存 Bubble Plot (动态命名)
    bubble_pdf_path <- paste0('results/figure/', sample_id, '_pop_cellchat_bubble.pdf')
    # 确保图赋值给 p1，这样 ggsave 才能正确保存它
    p1 <- netVisual_bubble(cellchat, sources.use = c("MG1"), remove.isolate = TRUE, angle.x = 45) +
      ggtitle(paste("MG1 Interactions -", sample_id))
    ggsave(bubble_pdf_path, plot = p1, width = 5, height = 5)
    cat("已保存 Bubble 图至:", bubble_pdf_path, "\n")
    
  }, error = function(e) {
    # 如果发生错误，打印报错信息并继续下一个样本
    cat("处理样本", sample_id, "时发生错误:\n")
    print(e$message)
    cat("跳过此样本，继续处理下一个...\n")
  })
}

cat("\n所有样本处理完毕！\n")





# MODULE : CellChat Interaction Visualization -------------------

## 1 Parameters --------------------------------------------------------------
# 加载所需包 (请确保已安装)
library(qs)
library(CellChat)
library(ggplot2)
library(patchwork)

# 文件路径设置
input_file <- 'data/processed_data/con_DLB_2_pop_cellchat.qs'
out_dir <- 'results/figure/'
out_file_circle_1  <- paste0(out_dir, "DLB_2_cellchat_circle_multi.pdf")
out_file_circle_2  <- paste0(out_dir, "DLB_2_cellchat_circle_olig.pdf")
out_file_circle_3  <- paste0(out_dir, "DLB_2_cellchat_circle_subset.pdf")
out_file_bubble    <- paste0(out_dir, "DLB_2_bubble_cellchat.pdf")

# 定义细胞群组 (方便后续随意更改)
clusters_source_multi <- c("MG0", "MG1", "MG2", "MG3", "MG4", "MG5")
cluster_olig          <- c("Oligodendrocyte")
cluster_source_single <- c("MG0")
clusters_target_main  <- c("Excitatory neuron","Inhibitory neuron","Astrocyte","Oligodendrocyte","OPC","Ependymal cell","Endothelial cell")

# 子集分析群组 (将上述群组结合)
subset_celltypes <- c(cluster_olig, clusters_source_multi)

# 气泡图关注的受体-配体对
lr_pairs_list <- c(
  "NCAM1_NCAM2", "NCAM1_NCAM1", "JAM3_JAM3", "CD6_ALCAM", "ALCAM_CD6",
  "CNTN2_CNTN2_CNTNAP2", "CNTN2_CNTN2", "MAG_MAG",
  "SEMA4D_PLXNB1", "CDH1_ITGA1_ITGB1", "CD86_CD28",
  "CD80_CD28", "HLA-B_CD8A","HLA-A_CD8A", "ICAM1_ITGAX_ITGB2"
)


# 统一的颜色配置 (通过名称匹配，支持自动提取子集颜色)
cell_colors <- c(
  "MG0" = "#4DAF4A", 
  "MG1" = "#A65628",
  "MG2" = "#377EB8", 
  "MG3" = "#984EA3", 
  "MG4" = "#E6AB02", 
  "MG5" = "#8B0000", 
  "Excitatory neuron"         = "#E31A1C",
  "Inhibitory neuron"         = "#E32636",
  "Astrocyte"      = "#FF7F00",
  "Oligodendrocyte"            = "#1F78B4",
  "OPC"             = "#A6CEE3",
  "Ependymal cell"      = "#6A3D9A",
  "Endothelial cell"       = "#FB9A99"
)

## 2 Process Data ------------------------------------------------------------
# 检查并创建输出目录
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 加载数据
cellchat <- qread(input_file)

# 提取并处理权重矩阵 (用于子集圆圈图)
cell_weight <- cellchat@net$weight

net_weight_no_self <- cellchat@net$weight
diag(net_weight_no_self) <- 0

mat_subset_rows <- cell_weight[rownames(cell_weight) %in% subset_celltypes, ]
mat_subset <- mat_subset_rows[, colnames(mat_subset_rows) %in% subset_celltypes]
diag(mat_subset) <- 0

# 处理气泡图所需的 L-R pairs 数据框
pairLR.use_df <- data.frame(interaction_name = lr_pairs_list, stringsAsFactors = FALSE)

# 提取子集作图专属颜色字典
new_col_subset <- cell_colors[names(cell_colors) %in% subset_celltypes]


## 3 Generate Plots ----------------------------------------------------------
# --- Base R Plots (由于Base R直接写入设备，这里用函数封装作图逻辑) ---
plot_circle_1 <- function() {
  netVisual_circle(net_weight_no_self,
                   arrow.size = 0.5, arrow.width = 0.8, vertex.weight = 13,
                   weight.scale = TRUE, label.edge = FALSE,
                   color.use = cell_colors, title.name = "Interaction weights/strength",
                   sources.use = clusters_source_multi)
}

plot_circle_2 <- function() {
  netVisual_circle(net_weight_no_self,
                   arrow.size = 0.5, arrow.width = 0.8, vertex.weight = 13,
                   weight.scale = TRUE, label.edge = FALSE,
                   color.use = cell_colors, title.name = "Interaction weights/strength",
                   sources.use = cluster_olig)
}

plot_circle_3 <- function() {
  netVisual_circle(mat_subset,
                   arrow.size = 1, arrow.width = 0.8, vertex.weight = 13,
                   weight.scale = TRUE, label.edge = FALSE,
                   color.use = new_col_subset, title.name = "Interaction weights/strength",
                   edge.width.max = 15)
}

# --- ggplot2 Plots (气泡图) ---
p1 <- netVisual_bubble(cellchat,
                       sources.use = cluster_source_single,
                       targets.use = clusters_target_main,
                       pairLR.use = pairLR.use_df,
                       remove.isolate = FALSE, show.legend = FALSE) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

p2 <- netVisual_bubble(cellchat,
                       sources.use = clusters_source_multi,
                       targets.use = cluster_olig,
                       pairLR.use = pairLR.use_df,
                       remove.isolate = FALSE, show.legend = FALSE) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p3 <- netVisual_bubble(cellchat,
                       sources.use = cluster_olig,
                       targets.use = clusters_source_multi,
                       pairLR.use = pairLR.use_df,
                       remove.isolate = FALSE) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

# 拼接气泡图
combined_bubble <- p1 | p2 | p3


## 4 Save Results ------------------------------------------------------------
# 保存圆圈图 1
pdf(out_file_circle_1, width = 5, height = 5)
plot_circle_1()
invisible(dev.off())

# 保存圆圈图 2 (修复了原代码中覆盖同名文件的问题)
pdf(out_file_circle_2, width = 5, height = 5)
plot_circle_2()
invisible(dev.off())

# 保存圆圈图 3 (子集)
pdf(out_file_circle_3, width = 5, height = 5)
plot_circle_3()
invisible(dev.off())

# 保存拼接后的气泡图
pdf(out_file_bubble, height = 3.5, width = 6.5)
print(combined_bubble)
invisible(dev.off())



# Figure 4 ------------------------------------------------------------

# F4A  Cell-cell-contact circle network, MG0-MG5 sources and major brain-cell targets
# Source: source/project/code_upload/CellChat.r
# CellChat ----------------
library(CellChat)
library(Seurat)
library(qs)
library(ggplot2)
library(patchwork)
library(igraph)
library(ComplexHeatmap)

# Data Loading & Initialization
seu <- qread('data/processed_data/stroke_cellchat.qs')
seu.data.input <- GetAssayData(seu, assay = "RNA", slot = "data")
seu.identity <- subset(seu@meta.data, select = "cellchat_type")

cellchat <- createCellChat(object = seu.data.input, meta = seu.identity, group.by = "cellchat_type")
cellchat <- addMeta(cellchat, meta = seu.identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") 
groupSize <- as.numeric(table(cellchat@idents))

# CellChat Database Setup
CellChatDB <- CellChatDB.human  
CellChatDB.use <- subsetDB(CellChatDB, search = c("Cell-Cell Contact"))
cellchat@DB <- CellChatDB.use

# Core CellChat Processing
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Save Processed Object
qsave(cellchat, 'data/processed_data/con_stroke_cellchat.qs')

# Visualization Setup
out_dir <- 'results/figure/'
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

clusters_source_multi <- c("MG0", "MG1", "MG2", "MG3", "MG4", "MG5")
cluster_olig          <- c("OLIG")
cluster_source_single <- c("MG0")
clusters_target_main  <- c("Neuron", "Astrocyte", "OLIG", "OPC",  "Ependymal", "Vascular", "T cell")
subset_celltypes      <- c(cluster_olig, clusters_source_multi)

lr_pairs_list <- c(
  "NCAM1_NCAM2", "NCAM1_NCAM1", "JAM3_JAM3", "CD6_ALCAM", "ALCAM_CD6",
  "CNTN2_CNTN2_CNTNAP2", "CNTN2_CNTN2", "MAG_MAG",
  "SEMA4D_PLXNB1", "CDH1_ITGA1_ITGB1", "CD86_CD28",
  "CD80_CD28", "HLA-B_CD8A","HLA-A_CD8A", "ICAM1_ITGAX_ITGB2"
)

cell_colors <- c(
  "MG0" = "#4DAF4A", 
  "MG1" = "#A65628",
  "MG2" = "#377EB8", 
  "MG3" = "#984EA3", 
  "MG4" = "#E6AB02", 
  "MG5" = "#8B0000", 
  "Neuron"         = "#E31A1C",
  "Astrocyte"      = "#FF7F00",
  "OLIG"           = "#1F78B4",
  "OPC"            = "#A6CEE3",
  "Ependymal"      = "#6A3D9A",
  "Vascular"       = "#FB9A99",
  "T cell"         = "#FDBF6F"
)

# Matrix Preparation
cell_weight <- cellchat@net$weight
net_weight_no_self <- cell_weight
diag(net_weight_no_self) <- 0

mat_subset_rows <- cell_weight[rownames(cell_weight) %in% subset_celltypes, ]
mat_subset <- mat_subset_rows[, colnames(mat_subset_rows) %in% subset_celltypes]
diag(mat_subset) <- 0

pairLR.use_df <- data.frame(interaction_name = lr_pairs_list, stringsAsFactors = FALSE)
new_col_subset <- cell_colors[names(cell_colors) %in% subset_celltypes]


# Fig.01 ---------------
pdf(paste0(out_dir, "stroke_cellchat_circle_multi.pdf"), width = 5, height = 5)
netVisual_circle(net_weight_no_self,
                 arrow.size = 0.5, arrow.width = 0.8, vertex.weight = 13,
                 weight.scale = TRUE, label.edge = FALSE,
                 color.use = cell_colors, title.name = "Interaction weights/strength",
                 sources.use = clusters_source_multi)
invisible(dev.off())

# F4B  Three-block selected ligand-receptor bubbles: OLM->major cells, MG->OLIG, OLIG->MG
# Source: source/project/code_upload/CellChat.r
# CellChat ----------------
library(CellChat)
library(Seurat)
library(qs)
library(ggplot2)
library(patchwork)
library(igraph)
library(ComplexHeatmap)
qsave(cellchat, 'data/processed_data/con_stroke_cellchat.qs')

# Visualization Setup
out_dir <- 'results/figure/'
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

clusters_source_multi <- c("MG0", "MG1", "MG2", "MG3", "MG4", "MG5")
cluster_olig          <- c("OLIG")
cluster_source_single <- c("MG0")
clusters_target_main  <- c("Neuron", "Astrocyte", "OLIG", "OPC",  "Ependymal", "Vascular", "T cell")
subset_celltypes      <- c(cluster_olig, clusters_source_multi)

lr_pairs_list <- c(
  "NCAM1_NCAM2", "NCAM1_NCAM1", "JAM3_JAM3", "CD6_ALCAM", "ALCAM_CD6",
  "CNTN2_CNTN2_CNTNAP2", "CNTN2_CNTN2", "MAG_MAG",
  "SEMA4D_PLXNB1", "CDH1_ITGA1_ITGB1", "CD86_CD28",
  "CD80_CD28", "HLA-B_CD8A","HLA-A_CD8A", "ICAM1_ITGAX_ITGB2"
)

cell_colors <- c(
  "MG0" = "#4DAF4A", 
  "MG1" = "#A65628",
  "MG2" = "#377EB8", 
  "MG3" = "#984EA3", 
  "MG4" = "#E6AB02", 
  "MG5" = "#8B0000", 
  "Neuron"         = "#E31A1C",
  "Astrocyte"      = "#FF7F00",
  "OLIG"           = "#1F78B4",
  "OPC"            = "#A6CEE3",
  "Ependymal"      = "#6A3D9A",
  "Vascular"       = "#FB9A99",
  "T cell"         = "#FDBF6F"
)

# Matrix Preparation
cell_weight <- cellchat@net$weight
net_weight_no_self <- cell_weight
diag(net_weight_no_self) <- 0

mat_subset_rows <- cell_weight[rownames(cell_weight) %in% subset_celltypes, ]
mat_subset <- mat_subset_rows[, colnames(mat_subset_rows) %in% subset_celltypes]
diag(mat_subset) <- 0

pairLR.use_df <- data.frame(interaction_name = lr_pairs_list, stringsAsFactors = FALSE)
new_col_subset <- cell_colors[names(cell_colors) %in% subset_celltypes]
# Fig.04 ---------------
p1 <- netVisual_bubble(cellchat,
                       sources.use = cluster_source_single,
                       targets.use = clusters_target_main,
                       pairLR.use = pairLR.use_df,
                       remove.isolate = FALSE, show.legend = FALSE) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

p2 <- netVisual_bubble(cellchat,
                       sources.use = clusters_source_multi,
                       targets.use = cluster_olig,
                       pairLR.use = pairLR.use_df,
                       remove.isolate = FALSE, show.legend = FALSE) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p3 <- netVisual_bubble(cellchat,
                       sources.use = cluster_olig,
                       targets.use = clusters_source_multi,
                       pairLR.use = pairLR.use_df,
                       remove.isolate = FALSE) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

combined_bubble <- p1 | p2 | p3

pdf(paste0(out_dir, "stroke_bubble_cellchat.pdf"), height = 3.5, width = 6.5)
print(combined_bubble)
invisible(dev.off())

# F4C  Integrated three-slice MERFISH microglial UMAP ST-OLM/ST_MG1-ST_MG5
# Source: source/project/Reference/combine/scripts/Spatial/scripts/spatial_module_verify_single.r
# Verification copy of spatial_module.r: single-plot reproduction only.
# The original merged object retained in this migration already contains
# plot_celltype and myelin1. No label, score, or source object is changed.

scp_lib <- file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library")
.libPaths(c(scp_lib, .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(SCP)
})

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
file_in_merfish <- file.path(project_dir, "data/processed_data/merfish_immune.combined_with4.rds")
dir_out_fig <- file.path(project_dir, "results/figure/code_verification_module")
dir.create(dir_out_fig, recursive = TRUE, showWarnings = FALSE)

# Original script settings retained verbatim.
col_new_celltype <- "plot_celltype"
score_feature_name <- "myelin1"
st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35",
  "ST_MG1" = "#4DBBD5",
  "ST_MG2" = "#00A087",
  "ST_MG3" = "#3C5488",
  "ST_MG4" = "#F39B7F",
  "ST_MG5" = "#7E6148"
)
target_order <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")

merfish <- readRDS(file_in_merfish)
stopifnot(all(c(col_new_celltype, score_feature_name) %in% colnames(merfish@meta.data)))
merfish$plot_celltype <- factor(merfish$plot_celltype, levels = target_order)

# Blue-box panel D: original DimPlot call, saved as a standalone figure.
p_umap <- DimPlot(merfish, group.by = col_new_celltype, cols = st_mg_colors_npg) +
  ggtitle("Spatial MG subclusters")

# Blue-box panel E: original FeatureStatPlot call, saved as a standalone figure.
celltype <- levels(merfish$plot_celltype)
my_comparisons <- lapply(celltype[-1], function(x) c(x, celltype[1]))
p_myelin <- SCP::FeatureStatPlot(
  srt = merfish,
  group.by = "plot_celltype",
  stat.by = c(score_feature_name),
  add_box = FALSE,
  comparisons = my_comparisons,
  bg.by = "plot_celltype",
  palcolor = st_mg_colors_npg,
  bg_palcolor = st_mg_colors_npg
)

ggsave(file.path(dir_out_fig, "spatial_MG_subclusters_UMAP_verified.pdf"), p_umap,
       width = 7, height = 6, bg = "white")
ggsave(file.path(dir_out_fig, "spatial_MG_subclusters_UMAP_verified.png"), p_umap,
       width = 7, height = 6, dpi = 300, bg = "white")
ggsave(file.path(dir_out_fig, "myelin_module_score_FeatureStatPlot_verified.pdf"), p_myelin,
       width = 7, height = 6, bg = "white")
ggsave(file.path(dir_out_fig, "myelin_module_score_FeatureStatPlot_verified.png"), p_myelin,
       width = 7, height = 6, dpi = 300, bg = "white")

message("Single-plot verification outputs: ", dir_out_fig)

# Source: source/project/code_upload/Spatial.r
# Spatial ----------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(spatstat)
library(qs)
library(SCP)

file_master_rds <- 'data/raw_data/merfish_integrated_v2.rds'
file_in_merfish_immune <- 'data/raw_data/merfish_immune.combined_with3.rds'

col_master_celltype <- 'cell.type_manual2'
target_colors <- c("ST_MG0" = "red", "OLG" = "blue", "OPC" = "yellow")
color_other <- "gray80"
merfish <- readRDS(file_in_merfish_immune)
mouse_mg_markers <- qread('../MultiSpecies/data/saved_data/mouse_mg_markers.qs')

celltype_dict <- c(
  'Myelin_MG' = 'ST_MG0',
  'MG_0'      = 'ST_MG1',
  'MG_1'      = 'ST_MG2',
  'MG_2'      = 'ST_MG3',
  'MG_3'      = 'ST_MG4',
  'MG_5'      = 'ST_MG5'
)

merfish@meta.data[['plot_celltype']] <- as.character(merfish@meta.data[['cell.type_manual']])
for (old_name in names(celltype_dict)) {
  new_name <- celltype_dict[[old_name]]
  merfish@meta.data[['plot_celltype']][merfish@meta.data[['plot_celltype']] %in% old_name] <- new_name
}

filtered_data <- mouse_mg_markers %>% filter(.data[['cluster']] == 'cluster2')
int_gene <- intersect(filtered_data[['gene']], rownames(merfish))

merfish <- AddModuleScore(object = merfish,
                          features = list(int_gene),
                          ctrl = 10,
                          name = "myelin")

saveRDS(merfish, file = 'data/processed_data/merfish_immune.combined_with4.rds')

st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35", 
  "ST_MG1" = "#4DBBD5", 
  "ST_MG2" = "#00A087", 
  "ST_MG3" = "#3C5488", 
  "ST_MG4" = "#F39B7F", 
  "ST_MG5" = "#7E6148"
)

# Fig.02 ---------------
p1 <- DimPlot(merfish, group.by = 'plot_celltype', cols = st_mg_colors_npg)

target_order <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")
merfish$plot_celltype <- factor(merfish$plot_celltype, levels = target_order)

celltype <- levels(merfish$plot_celltype)
my_comparisons <- lapply(celltype[-1], function(x) c(x, celltype[1]))

p2 <- FeatureStatPlot(
  srt = merfish, group.by = "plot_celltype",
  stat.by = c("myelin1"), add_box = FALSE,
  comparisons = my_comparisons,
  bg.by = 'plot_celltype',
  palcolor = st_mg_colors_npg,
  bg_palcolor = st_mg_colors_npg
)

combined_plot <- p1 | p2
ggsave(filename = 'results/figure/module_score_myelin.pdf', plot = combined_plot, width = 10, height = 5)

# F4D  myelin1 module-score violins across integrated MERFISH ST_MG0-ST_MG5
# Source: source/project/Reference/combine/scripts/Spatial/scripts/spatial_module_verify_single.r
# Verification copy of spatial_module.r: single-plot reproduction only.
# The original merged object retained in this migration already contains
# plot_celltype and myelin1. No label, score, or source object is changed.

scp_lib <- file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library")
.libPaths(c(scp_lib, .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(SCP)
})

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
file_in_merfish <- file.path(project_dir, "data/processed_data/merfish_immune.combined_with4.rds")
dir_out_fig <- file.path(project_dir, "results/figure/code_verification_module")
dir.create(dir_out_fig, recursive = TRUE, showWarnings = FALSE)

# Original script settings retained verbatim.
col_new_celltype <- "plot_celltype"
score_feature_name <- "myelin1"
st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35",
  "ST_MG1" = "#4DBBD5",
  "ST_MG2" = "#00A087",
  "ST_MG3" = "#3C5488",
  "ST_MG4" = "#F39B7F",
  "ST_MG5" = "#7E6148"
)
target_order <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")

merfish <- readRDS(file_in_merfish)
stopifnot(all(c(col_new_celltype, score_feature_name) %in% colnames(merfish@meta.data)))
merfish$plot_celltype <- factor(merfish$plot_celltype, levels = target_order)

# Blue-box panel D: original DimPlot call, saved as a standalone figure.
p_umap <- DimPlot(merfish, group.by = col_new_celltype, cols = st_mg_colors_npg) +
  ggtitle("Spatial MG subclusters")

# Blue-box panel E: original FeatureStatPlot call, saved as a standalone figure.
celltype <- levels(merfish$plot_celltype)
my_comparisons <- lapply(celltype[-1], function(x) c(x, celltype[1]))
p_myelin <- SCP::FeatureStatPlot(
  srt = merfish,
  group.by = "plot_celltype",
  stat.by = c(score_feature_name),
  add_box = FALSE,
  comparisons = my_comparisons,
  bg.by = "plot_celltype",
  palcolor = st_mg_colors_npg,
  bg_palcolor = st_mg_colors_npg
)

ggsave(file.path(dir_out_fig, "spatial_MG_subclusters_UMAP_verified.pdf"), p_umap,
       width = 7, height = 6, bg = "white")
ggsave(file.path(dir_out_fig, "spatial_MG_subclusters_UMAP_verified.png"), p_umap,
       width = 7, height = 6, dpi = 300, bg = "white")
ggsave(file.path(dir_out_fig, "myelin_module_score_FeatureStatPlot_verified.pdf"), p_myelin,
       width = 7, height = 6, bg = "white")
ggsave(file.path(dir_out_fig, "myelin_module_score_FeatureStatPlot_verified.png"), p_myelin,
       width = 7, height = 6, dpi = 300, bg = "white")

message("Single-plot verification outputs: ", dir_out_fig)

# Source: source/project/code_upload/Spatial.r
# Spatial ----------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(spatstat)
library(qs)
library(SCP)

file_master_rds <- 'data/raw_data/merfish_integrated_v2.rds'
file_in_merfish_immune <- 'data/raw_data/merfish_immune.combined_with3.rds'

col_master_celltype <- 'cell.type_manual2'
target_colors <- c("ST_MG0" = "red", "OLG" = "blue", "OPC" = "yellow")
color_other <- "gray80"
merfish <- readRDS(file_in_merfish_immune)
mouse_mg_markers <- qread('../MultiSpecies/data/saved_data/mouse_mg_markers.qs')

celltype_dict <- c(
  'Myelin_MG' = 'ST_MG0',
  'MG_0'      = 'ST_MG1',
  'MG_1'      = 'ST_MG2',
  'MG_2'      = 'ST_MG3',
  'MG_3'      = 'ST_MG4',
  'MG_5'      = 'ST_MG5'
)

merfish@meta.data[['plot_celltype']] <- as.character(merfish@meta.data[['cell.type_manual']])
for (old_name in names(celltype_dict)) {
  new_name <- celltype_dict[[old_name]]
  merfish@meta.data[['plot_celltype']][merfish@meta.data[['plot_celltype']] %in% old_name] <- new_name
}

filtered_data <- mouse_mg_markers %>% filter(.data[['cluster']] == 'cluster2')
int_gene <- intersect(filtered_data[['gene']], rownames(merfish))

merfish <- AddModuleScore(object = merfish,
                          features = list(int_gene),
                          ctrl = 10,
                          name = "myelin")

saveRDS(merfish, file = 'data/processed_data/merfish_immune.combined_with4.rds')

st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35", 
  "ST_MG1" = "#4DBBD5", 
  "ST_MG2" = "#00A087", 
  "ST_MG3" = "#3C5488", 
  "ST_MG4" = "#F39B7F", 
  "ST_MG5" = "#7E6148"
)

# Fig.02 ---------------
p1 <- DimPlot(merfish, group.by = 'plot_celltype', cols = st_mg_colors_npg)

target_order <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")
merfish$plot_celltype <- factor(merfish$plot_celltype, levels = target_order)

celltype <- levels(merfish$plot_celltype)
my_comparisons <- lapply(celltype[-1], function(x) c(x, celltype[1]))

p2 <- FeatureStatPlot(
  srt = merfish, group.by = "plot_celltype",
  stat.by = c("myelin1"), add_box = FALSE,
  comparisons = my_comparisons,
  bg.by = 'plot_celltype',
  palcolor = st_mg_colors_npg,
  bg_palcolor = st_mg_colors_npg
)

combined_plot <- p1 | p2
ggsave(filename = 'results/figure/module_score_myelin.pdf', plot = combined_plot, width = 10, height = 5)

# F4E  s2r2 MERFISH ROI and magnified inset; OLG blue, OPC yellow, ST OLM red, others gray
# Source: source/project/Reference/combine/scripts/Spatial/scripts/plot.r
# MODULE : Spatial Mapping, Subsetting, and Visualization -------------------

## 1 Parameters --------------------------------------------------------------
# 加载依赖包
library(Seurat)
library(dplyr)
library(ggplot2)

setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/"))

# 主数据路径与列名配置 (更新为 v2 版本及 cell.type_manual2)
file_master_rds <- Sys.getenv("MG_EXTERNAL_022")
col_master_celltype <- 'cell.type_manual2'

# 图片输出目录
dir_out_fig <- 'results/figure/'

# 细胞类型颜色配置
target_colors <- c(
  "ST_MG0" = "red",
  "OLG"    = "blue",
  "OPC"    = "yellow"
)
color_other <- "gray80"

# 绘图基础参数
plot_pt_size <- 0.5
plot_width   <- 6
plot_height  <- 6

# 样本级配置字典 (输出文件统一更新为 _v2 命名)
sample_configs <- list(
  "s2r1" = list(
    in_rds  = Sys.getenv("MG_EXTERNAL_023"),
    out_rds = "data/processed_data/merfish_s2r1_v2.rds",
    coord   = c(x_min = -3000, x_max = -1500, y_min = -8000, y_max = -5300)
  ),
  "s2r2" = list(
    in_rds  = Sys.getenv("MG_EXTERNAL_024"),
    out_rds = "data/processed_data/merfish_s2r2_v2.rds",
    coord   = c(x_min = 800, x_max = 4000, y_min = -10000, y_max = -7500)
  ),
  "s2r3" = list(
    in_rds  = Sys.getenv("MG_EXTERNAL_025"),
    out_rds = "data/processed_data/merfish_s2r3_v2.rds",
    coord   = c(x_min = 100, x_max = 2500, y_min = -8000, y_max = -4500)
  )
)


## 2 Process Data ------------------------------------------------------------
# 读取全局数据
master_merfish <- readRDS(file_master_rds)

# 严格采用新逻辑：提取为主数据框并重命名列
master_celltype_df <- as.data.frame(master_merfish[[col_master_celltype]])
colnames(master_celltype_df) <- 'CellType'

# 释放主数据内存
rm(master_merfish)
gc()

# 初始化存储列表
processed_objects <- list()
subset_objects <- list()

# 批量处理数据 (映射 CellType + 截取坐标子集)
for (smpl in names(sample_configs)) {
  cfg <- sample_configs[[smpl]]
  
  # 1. 读取单样本数据
  obj <- readRDS(cfg$in_rds)
  
  # 2. 映射 CellType (严格采用新逻辑：保留 dataframe 格式赋值)
  smpl_celltype <- master_celltype_df[rownames(master_celltype_df) %in% colnames(obj), , drop = FALSE]
  obj$CellType <- smpl_celltype
  processed_objects[[smpl]] <- obj # 存入列表，留待第4步保存
  
  # 3. 筛选坐标子集
  df_meta <- obj@meta.data %>% select(rotated_x, rotated_y, CellType)
  cells_keep <- df_meta %>% 
    filter(rotated_x >= cfg$coord["x_min"] & rotated_x <= cfg$coord["x_max"] & 
             rotated_y >= cfg$coord["y_min"] & rotated_y <= cfg$coord["y_max"]) %>% 
    rownames()
  
  subset_objects[[smpl]] <- subset(obj, cells = cells_keep)
}


## 3 Generate Plots ----------------------------------------------------------
# 检查并创建图片输出目录
if(!dir.exists(dir_out_fig)) dir.create(dir_out_fig, recursive = TRUE)

plot_list <- list()

# 批量生成图片
for (smpl in names(subset_objects)) {
  obj_sub <- subset_objects[[smpl]]
  
  # 动态生成颜色向量
  current_celltypes <- unique(obj_sub$CellType)
  colors_mapped <- sapply(current_celltypes, function(x) {
    if (x %in% names(target_colors)) return(target_colors[[x]])
    else return(color_other)
  })
  
  # 绘制 DimPlot
  p <- DimPlot(obj_sub, reduction = 'spatial', group.by = 'CellType', pt.size = plot_pt_size) + 
    scale_color_manual(values = colors_mapped) +
    ggtitle(paste("Spatial Plot -", smpl))
  
  plot_list[[smpl]] <- p
}


## 4 Save Results ------------------------------------------------------------
# 批量保存数据和图片
for (smpl in names(sample_configs)) {
  cfg <- sample_configs[[smpl]]
  
  # 保存映射好 CellType 的完整 RDS 对象 (此时已经是 _v2.rds)
  saveRDS(processed_objects[[smpl]], file = cfg$out_rds)
  
  # 保存图片 (PDF & PNG)
  file_base <- paste0(dir_out_fig, smpl, "_plot_v2")
  ggsave(paste0(file_base, ".pdf"), plot = plot_list[[smpl]], width = plot_width, height = plot_height)
  ggsave(paste0(file_base, ".png"), plot = plot_list[[smpl]], width = plot_width, height = plot_height)
}

# Source: source/project/Reference/combine/scripts/Spatial/scripts/plot_zoom.r
# MODULE : Spatial Mapping, Subsetting, and Visualization -------------------

## 1 Parameters --------------------------------------------------------------
# 加载依赖包
library(Seurat)
library(dplyr)
library(ggplot2)

setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/"))

# 主数据路径与列名配置 (更新为 v2 版本及 cell.type_manual2)
file_master_rds <- Sys.getenv("MG_EXTERNAL_022")
col_master_celltype <- 'cell.type_manual2'

# 图片输出目录
dir_out_fig <- 'results/figure/'

# 细胞类型颜色配置
target_colors <- c(
  "ST_MG0" = "red",
  "OLG"    = "blue",
  "OPC"    = "yellow"
)
color_other <- "gray80"

# 绘图基础参数
plot_pt_size <- 0.5
plot_width   <- 4
plot_height  <- 4

# 样本级配置字典 (输出文件统一更新为 _v2 命名)
sample_configs <- list(
  "s2r1" = list(
    in_rds  = Sys.getenv("MG_EXTERNAL_023"),
    out_rds = "data/processed_data/merfish_s2r1_v2.rds",
    coord   = c(x_min = -2000, x_max = -1500, y_min = -8000, y_max = -7500)
  ),
  "s2r2" = list(
    in_rds  = Sys.getenv("MG_EXTERNAL_024"),
    out_rds = "data/processed_data/merfish_s2r2_v2.rds",
    coord   = c(x_min = 1000, x_max = 1500, y_min = -9000, y_max = -8500)
  ),
  "s2r3" = list(
    in_rds  = Sys.getenv("MG_EXTERNAL_025"),
    out_rds = "data/processed_data/merfish_s2r3_v2.rds",
    coord   = c(x_min = 1000, x_max = 1500, y_min = -5500, y_max = -4500)
  )
)


## 2 Process Data ------------------------------------------------------------
# 读取全局数据
master_merfish <- readRDS(file_master_rds)

# 严格采用新逻辑：提取为主数据框并重命名列
master_celltype_df <- as.data.frame(master_merfish[[col_master_celltype]])
colnames(master_celltype_df) <- 'CellType'

# 释放主数据内存
rm(master_merfish)
gc()

# 初始化存储列表
processed_objects <- list()
subset_objects <- list()

# 批量处理数据 (映射 CellType + 截取坐标子集)
for (smpl in names(sample_configs)) {
  cfg <- sample_configs[[smpl]]
  
  # 1. 读取单样本数据
  obj <- readRDS(cfg$in_rds)
  
  # 2. 映射 CellType (严格采用新逻辑：保留 dataframe 格式赋值)
  smpl_celltype <- master_celltype_df[rownames(master_celltype_df) %in% colnames(obj), , drop = FALSE]
  obj$CellType <- smpl_celltype
  processed_objects[[smpl]] <- obj # 存入列表，留待第4步保存
  
  # 3. 筛选坐标子集
  df_meta <- obj@meta.data %>% select(rotated_x, rotated_y, CellType)
  cells_keep <- df_meta %>% 
    filter(rotated_x >= cfg$coord["x_min"] & rotated_x <= cfg$coord["x_max"] & 
             rotated_y >= cfg$coord["y_min"] & rotated_y <= cfg$coord["y_max"]) %>% 
    rownames()
  
  subset_objects[[smpl]] <- subset(obj, cells = cells_keep)
}


## 3 Generate Plots ----------------------------------------------------------
# 检查并创建图片输出目录
if(!dir.exists(dir_out_fig)) dir.create(dir_out_fig, recursive = TRUE)

plot_list <- list()

# 批量生成图片
for (smpl in names(subset_objects)) {
  obj_sub <- subset_objects[[smpl]]
  
  # 动态生成颜色向量
  current_celltypes <- unique(obj_sub$CellType)
  colors_mapped <- sapply(current_celltypes, function(x) {
    if (x %in% names(target_colors)) return(target_colors[[x]])
    else return(color_other)
  })
  
  # 绘制 DimPlot
  p <- DimPlot(obj_sub, reduction = 'spatial', group.by = 'CellType', pt.size = plot_pt_size) + 
    scale_color_manual(values = colors_mapped) +
    ggtitle(paste("Spatial Plot -", smpl))
  
  plot_list[[smpl]] <- p
}


## 4 Save Results ------------------------------------------------------------
# 批量保存数据和图片
for (smpl in names(sample_configs)) {
  cfg <- sample_configs[[smpl]]
  
  # 保存映射好 CellType 的完整 RDS 对象 (此时已经是 _v2.rds)
  saveRDS(processed_objects[[smpl]], file = cfg$out_rds)
  
  # 保存图片 (PDF & PNG)
  file_base <- paste0(dir_out_fig, smpl, "_plot_zoom")
  ggsave(paste0(file_base, ".pdf"), plot = plot_list[[smpl]], width = plot_width, height = plot_height)
  ggsave(paste0(file_base, ".png"), plot = plot_list[[smpl]], width = plot_width, height = plot_height)
}

# F4F  Cross-type centered L(r), OLG versus ST OLM and ST_MG1-ST_MG4
# Source: source/project/code_upload/Spatial.r
# Spatial ----------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(spatstat)
library(qs)
library(SCP)

file_master_rds <- 'data/raw_data/merfish_integrated_v2.rds'
file_in_merfish_immune <- 'data/raw_data/merfish_immune.combined_with3.rds'

col_master_celltype <- 'cell.type_manual2'
target_colors <- c("ST_MG0" = "red", "OLG" = "blue", "OPC" = "yellow")
color_other <- "gray80"
st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35", 
  "ST_MG1" = "#4DBBD5", 
  "ST_MG2" = "#00A087", 
  "ST_MG3" = "#3C5488", 
  "ST_MG4" = "#F39B7F", 
  "ST_MG5" = "#7E6148"
)
obj_l <- readRDS("data/processed_data/merfish_s2r1_v2.rds")
coord_filter <- c(x_min = -3000, x_max = -1500, y_min = -8000, y_max = -5300)

df_l <- obj_l@meta.data %>% 
  select(rotated_x, rotated_y, CellType) %>%
  filter(rotated_x >= coord_filter["x_min"] & rotated_x <= coord_filter["x_max"] & 
           rotated_y >= coord_filter["y_min"] & rotated_y <= coord_filter["y_max"])

df_l$rotated_x <- as.numeric(df_l$rotated_x)
df_l$rotated_y <- as.numeric(df_l$rotated_y)

pp <- ppp(df_l$rotated_x, df_l$rotated_y, 
          window = owin(range(df_l$rotated_x), range(df_l$rotated_y)), 
          marks = factor(df_l$CellType))

centered_L <- function(pp, i, j, r_max, step, nsim) {
  r_seq <- seq(0, r_max, by = step)
  L_sim <- envelope(pp, Lcross, i = i, j = j, nsim = nsim, r = r_seq, 
                    correction = "border", simulate = expression(rlabel(pp)),
                    verbose = FALSE) 
  
  if (is.null(L_sim$lo) || is.null(L_sim$hi) || is.null(L_sim$mmean)) return(NULL)
  
  data.frame(
    r  = r_seq, 
    L  = L_sim$obs - L_sim$mmean, 
    lo = L_sim$lo - L_sim$mmean, 
    hi = L_sim$hi - L_sim$mmean
  )
}

center_cell <- "OLG"
target_cells <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")
results <- list()

for (target in target_cells) {
  if (sum(df_l$CellType == center_cell) >= 10 && sum(df_l$CellType == target) >= 10) {
    res <- centered_L(pp, center_cell, target, r_max = 250, step = 1, nsim = 200)
    if (!is.null(res)) {
      res$cell_type1 <- center_cell
      res$cell_type2 <- target
      results[[paste(center_cell, target, sep = "_vs_")]] <- res
    }
  }
}

df_combined <- bind_rows(results)
df_filtered <- df_combined %>% filter(r >= 5)

# Fig.03 ---------------
p_l_func <- ggplot(df_filtered, aes(x = r)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = cell_type2), alpha = 0.3) +
  geom_line(aes(y = L, color = cell_type2), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ cell_type2, scales = "fixed", ncol = length(target_cells)) + 
  labs(x = "Distance (μm)", y = "Centered L(r)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  scale_color_manual(values = st_mg_colors_npg) +
  scale_fill_manual(values = st_mg_colors_npg)

ggsave(filename = "results/figure/s2r1_OLG_vs_other_MG_subclusters_L_function_plot_single_row.pdf", 
       plot = p_l_func, width = 12, height = 3.5, units = "in", device = cairo_pdf, dpi = 300)


# Supplementary Figure 6 ------------------------------------------------------------

# S6A  Whole s2r2 slice: OLG, OPC, ST-OLM and other cells
# Source: source/project/Reference/combine/scripts/Spatial/scripts/full_slice_four_colour_center_coordinates_1to1.R
project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
data_dir <- file.path(project_dir, "data/processed_data")
out_dir <- file.path(project_dir, "results/figure/code_verification_full_slice_four_colour_center_coordinates_1to1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_colours <- c(
  "Other cells" = "grey80",
  "OLG" = "blue",
  "OPC" = "yellow",
  "ST_MG0" = "red"
)
plot_levels <- names(plot_colours)

for (sample_id in c("s2r1", "s2r2", "s2r3")) {
  obj <- readRDS(file.path(data_dir, paste0("merfish_", sample_id, "_v2.rds")))
  stopifnot(all(c("CellType", "center_x", "center_y") %in% colnames(obj@meta.data)))

  # The original full-slice display labels horizontal axis as y and vertical axis
  # as x. No rotation, mirroring, or use of the Seurat spatial reduction occurs.
  df <- obj@meta.data %>%
    transmute(
      center_x = as.numeric(center_x),
      center_y = as.numeric(center_y),
      four_colour_celltype = if_else(
        CellType %in% c("OLG", "OPC", "ST_MG0"), as.character(CellType), "Other cells"
      )
    ) %>%
    mutate(four_colour_celltype = factor(four_colour_celltype, levels = plot_levels))

  p <- ggplot(df, aes(x = center_y, y = center_x, colour = four_colour_celltype)) +
    geom_point(size = 0.25) +
    scale_colour_manual(values = plot_colours, drop = FALSE) +
    coord_fixed(ratio = 1, expand = FALSE) +
    labs(
      title = paste0("Whole-slice raw-coordinate plot: ", sample_id),
      subtitle = "Original center_y (horizontal) and center_x (vertical); 1:1 coordinate scale",
      x = "y (center_y)", y = "x (center_x)", colour = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )

  stem <- file.path(out_dir, paste0(sample_id, "_whole_slice_OLG_OPC_ST_MG0_other_center_coordinates_1to1"))
  ggsave(paste0(stem, ".pdf"), p, width = 8, height = 8, bg = "white")
  ggsave(paste0(stem, ".png"), p, width = 8, height = 8, dpi = 300, bg = "white")
  rm(obj, df, p)
  gc()
}

# S6B  s2r1 microglial UMAP
# Source: source/project/Reference/combine/scripts/Spatial/scripts/run_per_slice_microglia_umap_from_integrated_logic.R
spatial_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
input_dir <- file.path(spatial_root, "data/processed_data/per_slice_microglia_ST_MG0_MG5")
output_dir <- file.path(spatial_root, "data/processed_data/per_slice_microglia_ST_MG0_MG5_RNA_PCA_UMAP")
figure_dir <- file.path(spatial_root, "results/per_slice_microglia_ST_MG0_MG5_RNA_PCA_UMAP")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Recovered from merfish_immune.combined_with4.rds command log.
seed_use <- 42L
n_pcs <- 50L
umap_dims <- 1:10
cluster_resolutions <- seq(0.1, 2, 0.1)
target_order <- paste0("ST_MG", 0:5)
st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35",
  "ST_MG1" = "#4DBBD5",
  "ST_MG2" = "#00A087",
  "ST_MG3" = "#3C5488",
  "ST_MG4" = "#F39B7F",
  "ST_MG5" = "#7E6148"
)

slices <- c("s2r1", "s2r2", "s2r3")
audit_rows <- list()

for (slice_id in slices) {
  input_file <- file.path(input_dir, paste0("merfish_", slice_id, "_microglia_ST_MG0_MG5_from_v2.rds"))
  output_file <- file.path(output_dir, paste0("merfish_", slice_id, "_microglia_RNA_PCA50_UMAP_PC1_10.rds"))
  pdf_file <- file.path(figure_dir, paste0(slice_id, "_microglia_ST_MG0_MG5_UMAP.pdf"))
  png_file <- file.path(figure_dir, paste0(slice_id, "_microglia_ST_MG0_MG5_UMAP.png"))

  if (!file.exists(input_file)) stop("Input object not found: ", input_file)
  obj <- readRDS(input_file)
  if (!"CellType" %in% colnames(obj@meta.data)) stop("CellType is absent: ", input_file)
  if (!all(as.character(obj$CellType) %in% target_order)) {
    stop("Input contains a non-ST_MG0-ST_MG5 label: ", input_file)
  }

  # There is no integrated assay in a single brain slice; RNA is the valid
  # within-slice analysis assay. The original per-slice object is not modified.
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = min(2000L, nrow(obj)), verbose = FALSE)
  obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)

  # Exact downstream dimensionality and random seeds recovered from the
  # original integrated object's Seurat command log.
  pca_features <- VariableFeatures(obj)
  n_pcs_actual <- min(n_pcs, length(pca_features), ncol(obj) - 1L)
  obj <- RunPCA(obj, features = pca_features, npcs = n_pcs_actual,
                verbose = FALSE, seed.use = seed_use)
  dims_actual <- seq_len(min(max(umap_dims), ncol(Embeddings(obj, "pca"))))
  obj <- FindNeighbors(obj, reduction = "pca", dims = dims_actual, verbose = FALSE)
  obj <- FindClusters(obj, resolution = cluster_resolutions,
                      random.seed = 0L, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = dims_actual,
                 seed.use = seed_use, verbose = FALSE)

  # The original UMAP drawing call is retained; plot_celltype is a copied,
  # ordered alias of the original CellType labels for compatibility with it.
  obj$plot_celltype <- factor(as.character(obj$CellType), levels = target_order)
  p_umap <- DimPlot(obj, group.by = "plot_celltype", cols = st_mg_colors_npg)
  ggsave(pdf_file, p_umap, width = 7, height = 6, bg = "white")
  ggsave(png_file, p_umap, width = 7, height = 6, dpi = 300, bg = "white")
  validObject(obj)
  saveRDS(obj, output_file)

  label_counts <- table(factor(as.character(obj$CellType), levels = target_order))
  audit_rows[[slice_id]] <- tibble(
    slice_id = slice_id,
    input_object = input_file,
    output_object = output_file,
    figure_pdf = pdf_file,
    cells = ncol(obj),
    genes = nrow(obj),
    variable_features = length(pca_features),
    pca_dimensions = ncol(Embeddings(obj, "pca")),
    umap_dimensions = ncol(Embeddings(obj, "umap")),
    selected_label = names(label_counts),
    selected_label_cells = as.integer(label_counts),
    normalization = "LogNormalize; scale.factor=10000",
    pca = paste0("RNA; all variable features; npcs=", n_pcs_actual, "; seed=", seed_use),
    neighbors = paste0("pca; dims=", paste(dims_actual, collapse = ",")),
    clusters = "resolution=0.1,0.2,...,2.0; random.seed=0",
    umap = paste0("pca; dims=", paste(dims_actual, collapse = ","), "; seed=", seed_use,
                  "; Seurat defaults: uwot, n.neighbors=30, min.dist=0.3, metric=cosine")
  )

  rm(obj)
  gc(verbose = FALSE)
}

# S6C  s2r1 myelin1 module scores
# Source: source/project/Reference/combine/scripts/Spatial/three_slice_ST_MG0_MG5_dual_module_score_20260916/scripts/01_mouse_cluster2_module_score_per_slice.R
#!/usr/bin/env Rscript

# Spatial microglia: original reference branch.
# This reproduces the module-score part of spatial_module.r for each supplied
# ST_MG0--ST_MG5 slice, using mouse_mg_markers cluster2 and AddModuleScore
# with ctrl = 10. No input object is overwritten.

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/three_slice_ST_MG0_MG5_dual_module_score_20260916")
result_dir <- file.path(project_dir, "results", "01_mouse_cluster2_module_score")
figure_dir <- file.path(result_dir, "figures")
table_dir <- file.path(result_dir, "tables")
object_dir <- file.path(result_dir, "objects")
log_dir <- file.path(project_dir, "logs")
tmp_dir <- file.path(project_dir, "tmp")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir, XDG_CACHE_HOME = file.path(project_dir, "cache"))
dir.create(Sys.getenv("XDG_CACHE_HOME"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library"), file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(SCP)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260916)

slice_paths <- c(
  s2r1 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r1_microglia_ST_MG0_MG5_from_v2.rds"),
  s2r2 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r2_microglia_ST_MG0_MG5_from_v2.rds"),
  s2r3 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r3_microglia_ST_MG0_MG5_from_v2.rds")
)
mouse_marker_qs <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/saved_data/mouse_mg_markers.qs")
target_order <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")
st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35", "ST_MG1" = "#4DBBD5", "ST_MG2" = "#00A087",
  "ST_MG3" = "#3C5488", "ST_MG4" = "#F39B7F", "ST_MG5" = "#7E6148"
)

mouse_mg_markers <- qread(mouse_marker_qs)
mouse_cluster2_markers <- mouse_mg_markers[as.character(mouse_mg_markers$cluster) == "cluster2", , drop = FALSE]
mouse_cluster2_genes <- unique(as.character(mouse_cluster2_markers$gene))
if (length(mouse_cluster2_genes) == 0L) stop("No cluster2 genes were found in mouse_mg_markers.")

all_slice_intersection <- data.frame()
all_slice_cluster_summary <- data.frame()
all_slice_comparisons <- data.frame()

for (slice_id in names(slice_paths)) {
  slice_object <- readRDS(slice_paths[[slice_id]])
  DefaultAssay(slice_object) <- "RNA"
  if (!"CellType" %in% colnames(slice_object[[]])) stop(slice_id, " lacks CellType.")
  slice_object$plot_celltype <- factor(as.character(slice_object$CellType), levels = target_order)
  if (any(is.na(slice_object$plot_celltype))) stop(slice_id, " contains a CellType outside ST_MG0--ST_MG5.")
  if (!all(target_order %in% levels(droplevels(slice_object$plot_celltype)))) stop(slice_id, " lacks one or more ST_MG groups.")

  slice_gene_intersection <- intersect(mouse_cluster2_genes, rownames(slice_object[["RNA"]]))
  if (length(slice_gene_intersection) == 0L) stop(slice_id, " has no overlap with mouse cluster2 markers.")
  slice_intersection_table <- data.frame(
    slice = slice_id,
    mouse_cluster2_marker = mouse_cluster2_genes,
    present_in_slice_RNA = mouse_cluster2_genes %in% rownames(slice_object[["RNA"]]),
    stringsAsFactors = FALSE
  )
  all_slice_intersection <- rbind(all_slice_intersection, slice_intersection_table)

  slice_object <- AddModuleScore(
    object = slice_object,
    assay = "RNA",
    features = list(slice_gene_intersection),
    ctrl = 10,
    name = "myelin",
    seed = 20260916
  )
  score_column <- "myelin1"
  if (!score_column %in% colnames(slice_object[[]])) stop(slice_id, " did not produce myelin1.")

  comparison_groups <- levels(slice_object$plot_celltype)
  my_comparisons <- lapply(comparison_groups[-1], function(x) c(x, comparison_groups[1]))
  slice_violin <- FeatureStatPlot(
    srt = slice_object,
    assay = "RNA",
    slot = "data",
    group.by = "plot_celltype",
    stat.by = score_column,
    add_box = FALSE,
    comparisons = my_comparisons,
    bg.by = "plot_celltype",
    palcolor = st_mg_colors_npg,
    bg_palcolor = st_mg_colors_npg,
    legend.position = "right",
    xlab = "",
    theme_args = list(base_size = 10)
  )
  slice_violin <- slice_violin & theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.margin = ggplot2::margin(t = 16, r = 4, b = 4, l = 4, unit = "pt")
  )
  slice_violin <- slice_violin + patchwork::plot_annotation(title = paste0(slice_id, ": mouse cluster2 module score"))
  ggsave(file.path(figure_dir, paste0(slice_id, "_mouse_cluster2_myelin1_violin.png")), slice_violin, width = 8.5, height = 6.2, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, paste0(slice_id, "_mouse_cluster2_myelin1_violin.pdf")), slice_violin, width = 8.5, height = 6.2, device = grDevices::cairo_pdf, bg = "white")

  score_values <- as.numeric(slice_object[[score_column]][, 1])
  celltype_values <- as.character(slice_object$plot_celltype)
  slice_cluster_summary <- do.call(rbind, lapply(comparison_groups, function(group_name) {
    group_scores <- score_values[celltype_values == group_name]
    data.frame(
      slice = slice_id,
      cluster = group_name,
      cells = length(group_scores),
      mean_myelin1 = mean(group_scores),
      median_myelin1 = median(group_scores),
      sd_myelin1 = sd(group_scores),
      stringsAsFactors = FALSE
    )
  }))
  all_slice_cluster_summary <- rbind(all_slice_cluster_summary, slice_cluster_summary)

  slice_comparison_statistics <- do.call(rbind, lapply(comparison_groups[-1], function(comparison_cluster) {
    mg0_scores <- score_values[celltype_values == "ST_MG0"]
    other_scores <- score_values[celltype_values == comparison_cluster]
    test_result <- wilcox.test(other_scores, mg0_scores, alternative = "two.sided", exact = FALSE)
    data.frame(
      slice = slice_id,
      reference_cluster = "ST_MG0",
      comparison_cluster = comparison_cluster,
      ST_MG0_cells = length(mg0_scores),
      comparison_cells = length(other_scores),
      ST_MG0_mean_myelin1 = mean(mg0_scores),
      comparison_mean_myelin1 = mean(other_scores),
      wilcoxon_W = unname(test_result$statistic),
      p_value_unadjusted = test_result$p.value,
      stringsAsFactors = FALSE
    )
  }))
  slice_comparison_statistics$p_value_BH_within_slice <- p.adjust(slice_comparison_statistics$p_value_unadjusted, method = "BH")
  slice_comparison_statistics$plot_significance <- ifelse(slice_comparison_statistics$p_value_unadjusted <= 1e-4, "****", ifelse(slice_comparison_statistics$p_value_unadjusted <= 1e-3, "***", ifelse(slice_comparison_statistics$p_value_unadjusted <= 1e-2, "**", ifelse(slice_comparison_statistics$p_value_unadjusted <= 0.05, "*", "ns"))))
  all_slice_comparisons <- rbind(all_slice_comparisons, slice_comparison_statistics)

  slice_cell_audit <- data.frame(
    cell_barcode = colnames(slice_object),
    slice = slice_id,
    CellType = as.character(slice_object$CellType),
    plot_celltype = as.character(slice_object$plot_celltype),
    myelin1 = score_values,
    stringsAsFactors = FALSE
  )
  write.table(slice_cell_audit, file.path(table_dir, paste0(slice_id, "_mouse_cluster2_myelin1_per_cell.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(slice_object, file.path(object_dir, paste0(slice_id, "_ST_MG0_MG5_mouse_cluster2_myelin1_scored.rds")), compress = "gzip")
}

write.table(all_slice_intersection, file.path(table_dir, "mouse_cluster2_marker_intersection_by_slice.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_slice_cluster_summary, file.path(table_dir, "mouse_cluster2_myelin1_summary_by_slice_and_cluster.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_slice_comparisons, file.path(table_dir, "mouse_cluster2_myelin1_Wilcoxon_ST_MG0_vs_other_clusters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parameter_record <- data.frame(
  item = c("reference_script", "marker_reference", "reference_cluster", "score_method", "score_assay", "control_genes", "score_column", "grouping", "reference_group", "violin_method", "test"),
  value = c(
    file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/scripts/spatial_module.r"),
    mouse_marker_qs, "cluster2", "Seurat::AddModuleScore", "RNA", "10", "myelin1", "CellType (ST_MG0--ST_MG5)", "ST_MG0", "SCP::FeatureStatPlot", "two-sided Wilcoxon; raw p displayed as stars, BH p saved per slice"
  ),
  stringsAsFactors = FALSE
)
write.table(parameter_record, file.path(result_dir, "01_mouse_cluster2_parameter_record.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "01_mouse_cluster2_module_score_sessionInfo.txt"))
writeLines(as.character(Sys.time()), file.path(log_dir, "01_mouse_cluster2_module_score_completed_at.txt"))
message("Mouse cluster2 module-score branch completed.")

# S6D  Whole s2r1 slice: OLG, OPC, ST-OLM and other cells
# Source: source/project/Reference/combine/scripts/Spatial/scripts/full_slice_four_colour_center_coordinates_1to1.R
project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
data_dir <- file.path(project_dir, "data/processed_data")
out_dir <- file.path(project_dir, "results/figure/code_verification_full_slice_four_colour_center_coordinates_1to1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_colours <- c(
  "Other cells" = "grey80",
  "OLG" = "blue",
  "OPC" = "yellow",
  "ST_MG0" = "red"
)
plot_levels <- names(plot_colours)

for (sample_id in c("s2r1", "s2r2", "s2r3")) {
  obj <- readRDS(file.path(data_dir, paste0("merfish_", sample_id, "_v2.rds")))
  stopifnot(all(c("CellType", "center_x", "center_y") %in% colnames(obj@meta.data)))

  # The original full-slice display labels horizontal axis as y and vertical axis
  # as x. No rotation, mirroring, or use of the Seurat spatial reduction occurs.
  df <- obj@meta.data %>%
    transmute(
      center_x = as.numeric(center_x),
      center_y = as.numeric(center_y),
      four_colour_celltype = if_else(
        CellType %in% c("OLG", "OPC", "ST_MG0"), as.character(CellType), "Other cells"
      )
    ) %>%
    mutate(four_colour_celltype = factor(four_colour_celltype, levels = plot_levels))

  p <- ggplot(df, aes(x = center_y, y = center_x, colour = four_colour_celltype)) +
    geom_point(size = 0.25) +
    scale_colour_manual(values = plot_colours, drop = FALSE) +
    coord_fixed(ratio = 1, expand = FALSE) +
    labs(
      title = paste0("Whole-slice raw-coordinate plot: ", sample_id),
      subtitle = "Original center_y (horizontal) and center_x (vertical); 1:1 coordinate scale",
      x = "y (center_y)", y = "x (center_x)", colour = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )

  stem <- file.path(out_dir, paste0(sample_id, "_whole_slice_OLG_OPC_ST_MG0_other_center_coordinates_1to1"))
  ggsave(paste0(stem, ".pdf"), p, width = 8, height = 8, bg = "white")
  ggsave(paste0(stem, ".png"), p, width = 8, height = 8, dpi = 300, bg = "white")
  rm(obj, df, p)
  gc()
}

# S6E  s2r3 microglial UMAP
# Source: source/project/Reference/combine/scripts/Spatial/scripts/run_per_slice_microglia_umap_from_integrated_logic.R
spatial_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
input_dir <- file.path(spatial_root, "data/processed_data/per_slice_microglia_ST_MG0_MG5")
output_dir <- file.path(spatial_root, "data/processed_data/per_slice_microglia_ST_MG0_MG5_RNA_PCA_UMAP")
figure_dir <- file.path(spatial_root, "results/per_slice_microglia_ST_MG0_MG5_RNA_PCA_UMAP")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Recovered from merfish_immune.combined_with4.rds command log.
seed_use <- 42L
n_pcs <- 50L
umap_dims <- 1:10
cluster_resolutions <- seq(0.1, 2, 0.1)
target_order <- paste0("ST_MG", 0:5)
st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35",
  "ST_MG1" = "#4DBBD5",
  "ST_MG2" = "#00A087",
  "ST_MG3" = "#3C5488",
  "ST_MG4" = "#F39B7F",
  "ST_MG5" = "#7E6148"
)

slices <- c("s2r1", "s2r2", "s2r3")
audit_rows <- list()

for (slice_id in slices) {
  input_file <- file.path(input_dir, paste0("merfish_", slice_id, "_microglia_ST_MG0_MG5_from_v2.rds"))
  output_file <- file.path(output_dir, paste0("merfish_", slice_id, "_microglia_RNA_PCA50_UMAP_PC1_10.rds"))
  pdf_file <- file.path(figure_dir, paste0(slice_id, "_microglia_ST_MG0_MG5_UMAP.pdf"))
  png_file <- file.path(figure_dir, paste0(slice_id, "_microglia_ST_MG0_MG5_UMAP.png"))

  if (!file.exists(input_file)) stop("Input object not found: ", input_file)
  obj <- readRDS(input_file)
  if (!"CellType" %in% colnames(obj@meta.data)) stop("CellType is absent: ", input_file)
  if (!all(as.character(obj$CellType) %in% target_order)) {
    stop("Input contains a non-ST_MG0-ST_MG5 label: ", input_file)
  }

  # There is no integrated assay in a single brain slice; RNA is the valid
  # within-slice analysis assay. The original per-slice object is not modified.
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = min(2000L, nrow(obj)), verbose = FALSE)
  obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)

  # Exact downstream dimensionality and random seeds recovered from the
  # original integrated object's Seurat command log.
  pca_features <- VariableFeatures(obj)
  n_pcs_actual <- min(n_pcs, length(pca_features), ncol(obj) - 1L)
  obj <- RunPCA(obj, features = pca_features, npcs = n_pcs_actual,
                verbose = FALSE, seed.use = seed_use)
  dims_actual <- seq_len(min(max(umap_dims), ncol(Embeddings(obj, "pca"))))
  obj <- FindNeighbors(obj, reduction = "pca", dims = dims_actual, verbose = FALSE)
  obj <- FindClusters(obj, resolution = cluster_resolutions,
                      random.seed = 0L, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = dims_actual,
                 seed.use = seed_use, verbose = FALSE)

  # The original UMAP drawing call is retained; plot_celltype is a copied,
  # ordered alias of the original CellType labels for compatibility with it.
  obj$plot_celltype <- factor(as.character(obj$CellType), levels = target_order)
  p_umap <- DimPlot(obj, group.by = "plot_celltype", cols = st_mg_colors_npg)
  ggsave(pdf_file, p_umap, width = 7, height = 6, bg = "white")
  ggsave(png_file, p_umap, width = 7, height = 6, dpi = 300, bg = "white")
  validObject(obj)
  saveRDS(obj, output_file)

  label_counts <- table(factor(as.character(obj$CellType), levels = target_order))
  audit_rows[[slice_id]] <- tibble(
    slice_id = slice_id,
    input_object = input_file,
    output_object = output_file,
    figure_pdf = pdf_file,
    cells = ncol(obj),
    genes = nrow(obj),
    variable_features = length(pca_features),
    pca_dimensions = ncol(Embeddings(obj, "pca")),
    umap_dimensions = ncol(Embeddings(obj, "umap")),
    selected_label = names(label_counts),
    selected_label_cells = as.integer(label_counts),
    normalization = "LogNormalize; scale.factor=10000",
    pca = paste0("RNA; all variable features; npcs=", n_pcs_actual, "; seed=", seed_use),
    neighbors = paste0("pca; dims=", paste(dims_actual, collapse = ",")),
    clusters = "resolution=0.1,0.2,...,2.0; random.seed=0",
    umap = paste0("pca; dims=", paste(dims_actual, collapse = ","), "; seed=", seed_use,
                  "; Seurat defaults: uwot, n.neighbors=30, min.dist=0.3, metric=cosine")
  )

  rm(obj)
  gc(verbose = FALSE)
}

# S6F  s2r3 myelin1 module scores
# Source: source/project/Reference/combine/scripts/Spatial/three_slice_ST_MG0_MG5_dual_module_score_20260916/scripts/01_mouse_cluster2_module_score_per_slice.R
#!/usr/bin/env Rscript

# Spatial microglia: original reference branch.
# This reproduces the module-score part of spatial_module.r for each supplied
# ST_MG0--ST_MG5 slice, using mouse_mg_markers cluster2 and AddModuleScore
# with ctrl = 10. No input object is overwritten.

project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/three_slice_ST_MG0_MG5_dual_module_score_20260916")
result_dir <- file.path(project_dir, "results", "01_mouse_cluster2_module_score")
figure_dir <- file.path(result_dir, "figures")
table_dir <- file.path(result_dir, "tables")
object_dir <- file.path(result_dir, "objects")
log_dir <- file.path(project_dir, "logs")
tmp_dir <- file.path(project_dir, "tmp")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir, XDG_CACHE_HOME = file.path(project_dir, "cache"))
dir.create(Sys.getenv("XDG_CACHE_HOME"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_library"), file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(SCP)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260916)

slice_paths <- c(
  s2r1 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r1_microglia_ST_MG0_MG5_from_v2.rds"),
  s2r2 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r2_microglia_ST_MG0_MG5_from_v2.rds"),
  s2r3 = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/data/processed_data/per_slice_microglia_ST_MG0_MG5/merfish_s2r3_microglia_ST_MG0_MG5_from_v2.rds")
)
mouse_marker_qs <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/saved_data/mouse_mg_markers.qs")
target_order <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")
st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35", "ST_MG1" = "#4DBBD5", "ST_MG2" = "#00A087",
  "ST_MG3" = "#3C5488", "ST_MG4" = "#F39B7F", "ST_MG5" = "#7E6148"
)

mouse_mg_markers <- qread(mouse_marker_qs)
mouse_cluster2_markers <- mouse_mg_markers[as.character(mouse_mg_markers$cluster) == "cluster2", , drop = FALSE]
mouse_cluster2_genes <- unique(as.character(mouse_cluster2_markers$gene))
if (length(mouse_cluster2_genes) == 0L) stop("No cluster2 genes were found in mouse_mg_markers.")

all_slice_intersection <- data.frame()
all_slice_cluster_summary <- data.frame()
all_slice_comparisons <- data.frame()

for (slice_id in names(slice_paths)) {
  slice_object <- readRDS(slice_paths[[slice_id]])
  DefaultAssay(slice_object) <- "RNA"
  if (!"CellType" %in% colnames(slice_object[[]])) stop(slice_id, " lacks CellType.")
  slice_object$plot_celltype <- factor(as.character(slice_object$CellType), levels = target_order)
  if (any(is.na(slice_object$plot_celltype))) stop(slice_id, " contains a CellType outside ST_MG0--ST_MG5.")
  if (!all(target_order %in% levels(droplevels(slice_object$plot_celltype)))) stop(slice_id, " lacks one or more ST_MG groups.")

  slice_gene_intersection <- intersect(mouse_cluster2_genes, rownames(slice_object[["RNA"]]))
  if (length(slice_gene_intersection) == 0L) stop(slice_id, " has no overlap with mouse cluster2 markers.")
  slice_intersection_table <- data.frame(
    slice = slice_id,
    mouse_cluster2_marker = mouse_cluster2_genes,
    present_in_slice_RNA = mouse_cluster2_genes %in% rownames(slice_object[["RNA"]]),
    stringsAsFactors = FALSE
  )
  all_slice_intersection <- rbind(all_slice_intersection, slice_intersection_table)

  slice_object <- AddModuleScore(
    object = slice_object,
    assay = "RNA",
    features = list(slice_gene_intersection),
    ctrl = 10,
    name = "myelin",
    seed = 20260916
  )
  score_column <- "myelin1"
  if (!score_column %in% colnames(slice_object[[]])) stop(slice_id, " did not produce myelin1.")

  comparison_groups <- levels(slice_object$plot_celltype)
  my_comparisons <- lapply(comparison_groups[-1], function(x) c(x, comparison_groups[1]))
  slice_violin <- FeatureStatPlot(
    srt = slice_object,
    assay = "RNA",
    slot = "data",
    group.by = "plot_celltype",
    stat.by = score_column,
    add_box = FALSE,
    comparisons = my_comparisons,
    bg.by = "plot_celltype",
    palcolor = st_mg_colors_npg,
    bg_palcolor = st_mg_colors_npg,
    legend.position = "right",
    xlab = "",
    theme_args = list(base_size = 10)
  )
  slice_violin <- slice_violin & theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.margin = ggplot2::margin(t = 16, r = 4, b = 4, l = 4, unit = "pt")
  )
  slice_violin <- slice_violin + patchwork::plot_annotation(title = paste0(slice_id, ": mouse cluster2 module score"))
  ggsave(file.path(figure_dir, paste0(slice_id, "_mouse_cluster2_myelin1_violin.png")), slice_violin, width = 8.5, height = 6.2, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, paste0(slice_id, "_mouse_cluster2_myelin1_violin.pdf")), slice_violin, width = 8.5, height = 6.2, device = grDevices::cairo_pdf, bg = "white")

  score_values <- as.numeric(slice_object[[score_column]][, 1])
  celltype_values <- as.character(slice_object$plot_celltype)
  slice_cluster_summary <- do.call(rbind, lapply(comparison_groups, function(group_name) {
    group_scores <- score_values[celltype_values == group_name]
    data.frame(
      slice = slice_id,
      cluster = group_name,
      cells = length(group_scores),
      mean_myelin1 = mean(group_scores),
      median_myelin1 = median(group_scores),
      sd_myelin1 = sd(group_scores),
      stringsAsFactors = FALSE
    )
  }))
  all_slice_cluster_summary <- rbind(all_slice_cluster_summary, slice_cluster_summary)

  slice_comparison_statistics <- do.call(rbind, lapply(comparison_groups[-1], function(comparison_cluster) {
    mg0_scores <- score_values[celltype_values == "ST_MG0"]
    other_scores <- score_values[celltype_values == comparison_cluster]
    test_result <- wilcox.test(other_scores, mg0_scores, alternative = "two.sided", exact = FALSE)
    data.frame(
      slice = slice_id,
      reference_cluster = "ST_MG0",
      comparison_cluster = comparison_cluster,
      ST_MG0_cells = length(mg0_scores),
      comparison_cells = length(other_scores),
      ST_MG0_mean_myelin1 = mean(mg0_scores),
      comparison_mean_myelin1 = mean(other_scores),
      wilcoxon_W = unname(test_result$statistic),
      p_value_unadjusted = test_result$p.value,
      stringsAsFactors = FALSE
    )
  }))
  slice_comparison_statistics$p_value_BH_within_slice <- p.adjust(slice_comparison_statistics$p_value_unadjusted, method = "BH")
  slice_comparison_statistics$plot_significance <- ifelse(slice_comparison_statistics$p_value_unadjusted <= 1e-4, "****", ifelse(slice_comparison_statistics$p_value_unadjusted <= 1e-3, "***", ifelse(slice_comparison_statistics$p_value_unadjusted <= 1e-2, "**", ifelse(slice_comparison_statistics$p_value_unadjusted <= 0.05, "*", "ns"))))
  all_slice_comparisons <- rbind(all_slice_comparisons, slice_comparison_statistics)

  slice_cell_audit <- data.frame(
    cell_barcode = colnames(slice_object),
    slice = slice_id,
    CellType = as.character(slice_object$CellType),
    plot_celltype = as.character(slice_object$plot_celltype),
    myelin1 = score_values,
    stringsAsFactors = FALSE
  )
  write.table(slice_cell_audit, file.path(table_dir, paste0(slice_id, "_mouse_cluster2_myelin1_per_cell.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(slice_object, file.path(object_dir, paste0(slice_id, "_ST_MG0_MG5_mouse_cluster2_myelin1_scored.rds")), compress = "gzip")
}

write.table(all_slice_intersection, file.path(table_dir, "mouse_cluster2_marker_intersection_by_slice.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_slice_cluster_summary, file.path(table_dir, "mouse_cluster2_myelin1_summary_by_slice_and_cluster.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_slice_comparisons, file.path(table_dir, "mouse_cluster2_myelin1_Wilcoxon_ST_MG0_vs_other_clusters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parameter_record <- data.frame(
  item = c("reference_script", "marker_reference", "reference_cluster", "score_method", "score_assay", "control_genes", "score_column", "grouping", "reference_group", "violin_method", "test"),
  value = c(
    file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial/scripts/spatial_module.r"),
    mouse_marker_qs, "cluster2", "Seurat::AddModuleScore", "RNA", "10", "myelin1", "CellType (ST_MG0--ST_MG5)", "ST_MG0", "SCP::FeatureStatPlot", "two-sided Wilcoxon; raw p displayed as stars, BH p saved per slice"
  ),
  stringsAsFactors = FALSE
)
write.table(parameter_record, file.path(result_dir, "01_mouse_cluster2_parameter_record.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "01_mouse_cluster2_module_score_sessionInfo.txt"))
writeLines(as.character(Sys.time()), file.path(log_dir, "01_mouse_cluster2_module_score_completed_at.txt"))
message("Mouse cluster2 module-score branch completed.")

# S6G  Whole s2r3 slice: OLG, OPC, ST-OLM and other cells
# Source: source/project/Reference/combine/scripts/Spatial/scripts/full_slice_four_colour_center_coordinates_1to1.R
project_dir <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Spatial")
data_dir <- file.path(project_dir, "data/processed_data")
out_dir <- file.path(project_dir, "results/figure/code_verification_full_slice_four_colour_center_coordinates_1to1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_colours <- c(
  "Other cells" = "grey80",
  "OLG" = "blue",
  "OPC" = "yellow",
  "ST_MG0" = "red"
)
plot_levels <- names(plot_colours)

for (sample_id in c("s2r1", "s2r2", "s2r3")) {
  obj <- readRDS(file.path(data_dir, paste0("merfish_", sample_id, "_v2.rds")))
  stopifnot(all(c("CellType", "center_x", "center_y") %in% colnames(obj@meta.data)))

  # The original full-slice display labels horizontal axis as y and vertical axis
  # as x. No rotation, mirroring, or use of the Seurat spatial reduction occurs.
  df <- obj@meta.data %>%
    transmute(
      center_x = as.numeric(center_x),
      center_y = as.numeric(center_y),
      four_colour_celltype = if_else(
        CellType %in% c("OLG", "OPC", "ST_MG0"), as.character(CellType), "Other cells"
      )
    ) %>%
    mutate(four_colour_celltype = factor(four_colour_celltype, levels = plot_levels))

  p <- ggplot(df, aes(x = center_y, y = center_x, colour = four_colour_celltype)) +
    geom_point(size = 0.25) +
    scale_colour_manual(values = plot_colours, drop = FALSE) +
    coord_fixed(ratio = 1, expand = FALSE) +
    labs(
      title = paste0("Whole-slice raw-coordinate plot: ", sample_id),
      subtitle = "Original center_y (horizontal) and center_x (vertical); 1:1 coordinate scale",
      x = "y (center_y)", y = "x (center_x)", colour = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )

  stem <- file.path(out_dir, paste0(sample_id, "_whole_slice_OLG_OPC_ST_MG0_other_center_coordinates_1to1"))
  ggsave(paste0(stem, ".pdf"), p, width = 8, height = 8, bg = "white")
  ggsave(paste0(stem, ".png"), p, width = 8, height = 8, dpi = 300, bg = "white")
  rm(obj, df, p)
  gc()
}

# S6H  OLM proportion in three observations
# Source: source/recovered/experimental/S4D_OLM_proportion.R
suppressPackageStartupMessages({
  library(ggplot2)
})
output_dir <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "S4D")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
olm_data <- data.frame(
  group = "OLM",
  value = c(0.08244, 0.10910, 0.08240),
  point_x = c(0.94, 1.00, 1.06)
)
plot_olm <- ggplot(olm_data, aes(x = 1, y = value)) +
  geom_boxplot(
    aes(group = 1),
    width = 0.38,
    outlier.shape = NA,
    fill = "#1786C0",
    alpha = 0.35,
    color = "#1786C0",
    linewidth = 0.9,
    staplewidth = 0.55
  ) +
  geom_point(
    aes(x = point_x),
    shape = 21,
    size = 2.6,
    stroke = 0.9,
    fill = "#1786C0",
    color = "black"
  ) +
  scale_x_continuous(
    breaks = 1,
    labels = "OLM",
    limits = c(0.5, 1.5),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    name = expression(Trim2^"+" * Iba1^"+" * OLM/Iba1^"+" ~ MG),
    breaks = seq(0, 0.12, 0.02),
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL) +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(linewidth = 0.9),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    plot.margin = margin(8, 8, 7, 7)
  ) +
  coord_cartesian(ylim = c(0, 0.12))
print(plot_olm)
ggsave(
  file.path(output_dir, "OLM_三值箱线图.png"),
  plot_olm,
  width = 3.0,
  height = 4.3,
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(output_dir, "OLM_三值箱线图.pdf"),
  plot_olm,
  width = 3.0,
  height = 4.3,
  device = cairo_pdf
)
cat("输入值（%）：", paste(olm_data$value * 100, collapse = ", "), "\n")
cat("输出目录：", output_dir, "\n", sep = "")


# Figure 5 ------------------------------------------------------------

# F5A  Cell-cycle phase on the human microglia UMAP
# Source: source/project/code_upload/Fig.01.r
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

# Source: source/project/Reference/combine/scripts/04.module.r
# MODULE：cell cycle------------------
tmp_obj <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

DefaultAssay(tmp_obj) <- "RNA"
tmp_obj <- CellCycleScoring(tmp_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# 活性递进配色
phase_colors <- c("G1" = "#1F77B4", "S" = "#FF7F0E", "G2M" = "#2CA02C")
p_cycle <- DimPlot(tmp_obj, group.by = 'Phase', cols = phase_colors)
ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/results/figure/cellcycle_dimplot.pdf"),p_cycle,width = 6,height = 5)

# F5B  Monocle3 pseudotime trajectory
# Source: source/project/code_upload/Pseudotime.r
sco <- qread('data/processed_data/mg_cca_meta_v2.qs')
data <- GetAssayData(sco, assay = 'RNA', layer  = 'counts')
cell_metadata <- sco@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))

cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, method = "PCA")
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = 'PCA')

int.embed <- Embeddings(sco, reduction = "umap")
int.embed <- int.embed[rownames(cds@int_colData$reducedDims$UMAP), ]
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds, k = 100, cluster_method = 'louvain')
cds <- learn_graph(cds, learn_graph_control = list(euclidean_distance_ratio = 0.8))
cds <- order_cells(cds) 

qsave(cds, 'scripts/Pseudotime/data/processed_data/save_cds.qs')

# Fig.01 ---------------
p_traj <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                     label_leaves = FALSE, label_branch_points = FALSE) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  NoAxes() +
  ggtitle("Trajectory Pseudotime") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("scripts/Pseudotime/results/figure/Trajectory_Pseudotime.pdf", plot = p_traj, width = 5, height = 4)
ggsave("scripts/Pseudotime/results/figure/Trajectory_Pseudotime.png", plot = p_traj, width = 5, height = 4)

# Source: source/project/Reference/combine/scripts/Pseudotime/scripts/module_pseudotime.r
#读取数据----------------------
sco <- qread('data/processed_data/mg_cca_meta_v2.qs')

#获取表达矩阵----------------
data <- GetAssayData(sco, assay = 'RNA', layer  = 'counts')

#细胞注释信息------------
cell_metadata <- sco@meta.data

#基因名----------------
gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))

#创建一个新的 cell_data_set 对象-----------------
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

#预处理cds-------------------
cds <- preprocess_cds(cds, method = "PCA")
#降维------------------
cds <- reduce_dimension(cds, reduction_method = "UMAP",preprocess_method = 'PCA')
#画图
# p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="sub_clusters2", show_trajectory_graph = FALSE) + ggtitle('cds.umap')

#将seurat对象的UMAP导入
int.embed <- Embeddings(sco, reduction = "umap")
#排序
int.embed <- int.embed[rownames(cds@int_colData$reducedDims$UMAP),]
#导入-------------------
cds@int_colData$reducedDims$UMAP <- int.embed

#聚类分区，不同分区的细胞会进行单独的轨迹分析
cds <- cluster_cells(cds,k = 100,cluster_method = 'louvain')

#构建细胞轨迹
cds <- learn_graph(cds, learn_graph_control = list(euclidean_distance_ratio = 0.8))

#发育轨迹（拟时序）排列细胞。
cds <- order_cells(cds)

qsave(cds,'scripts/Pseudotime/data/processed_data/save_cds.qs')
cds <- qread('scripts/Pseudotime/data/processed_data/save_cds.qs')

#画拟时序umap图
p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
               label_leaves = FALSE,  label_branch_points = FALSE) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  NoAxes() +
  ggtitle("Trajectory_Pseudotime") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) # 居中标题
p
ggsave("scripts/Pseudotime/results/figure/Trajectory_Pseudotime.pdf", plot = p, width = 5, height = 4)
ggsave("scripts/Pseudotime/results/figure/Trajectory_Pseudotime.png", plot = p, width = 5, height = 4)

# F5C  Pseudotime heatmap, cell-type bar, and selected GO labels
# Source: source/project/code_upload/Pseudotime.r
Track_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 50)
write.csv(Track_genes, "scripts/Pseudotime/results/table/Track_genes_all.csv", row.names = FALSE)
qsave(Track_genes, 'scripts/Pseudotime/data/processed_data/Track_genes.qs')

genes <- row.names(subset(Track_genes, q_value == 0 & morans_I > 0.15))
mat <- pre_pseudotime_matrix(cds_obj = cds, gene_list = genes)
qsave(mat, 'scripts/Pseudotime/data/processed_data/mat.qs')

ck <- clusterData(exp = mat, cluster.method = "kmeans", cluster.num = 5)
qsave(ck, 'scripts/Pseudotime/data/processed_data/ck.qs')

my_base_colors <- rev(c("#810000", "#CE1212", "#F05454", "#ffd06f",
                        "#ffe6b7", "#aadce0", "#72bcd5", "#528fad",
                        "#376795", "#1e466e"))

my_anno_colors <- c(
  "MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A",
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628"
)

# Fig.02 ---------------
pdf('scripts/Pseudotime/results/figure/monocle3_2.pdf', height = 10, width = 8, onefile = FALSE)
visCluster(object = ck,
           add.sampleanno = FALSE,
           plot.type = "both",
           cluster_rows = FALSE,
           cluster.order = c(1, 2, 5, 4, 3),
           markGenes = sample(rownames(mat), 30, replace = FALSE),
           ht.col.list = list(
             col_range = seq(-2, 2, length.out = 10), 
             col_color = my_base_colors
           ),
           line.side = "left")
dev.off()

enrich <- enrichCluster(object = ck,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 100,
                        seed = 24)

groups <- unique(enrich$group)
for (g in groups) {
  tmp_data <- enrich |> dplyr::filter(group == g)
  file_name <- paste0("scripts/Pseudotime/results/table/GO_Result_", g, ".csv")
  write.csv(tmp_data, file = file_name, row.names = FALSE)
}

palette <- c("Grays", "Light Grays", "Blues 2", "Blues 3", "Purples 2")

gglist <- lapply(seq_along(unique(enrich$group)), function(x){
  tmp <- enrich |> dplyr::filter(group == unique(enrich$group)[x]) |>
    dplyr::arrange(desc(pvalue))
  
  tmp$Description <- factor(tmp$Description, levels = tmp$Description)
  
  p <- ggplot(tmp) +
    geom_col(aes(x = -log10(pvalue), y = Description, fill = -log10(pvalue)), width = 0.75) +
    geom_line(aes(x = log10(ratio), y = as.numeric(Description)), color = "grey50") +
    geom_point(aes(x = log10(ratio), y = Description), size = 3, color = "orange") +
    theme_bw() +
    scale_y_discrete(position = "right", labels = function(x) stringr::str_wrap(x, width = 40)) +
    scale_x_continuous(sec.axis = sec_axis(~., name = "log10(ratio)")) +
    colorspace::scale_fill_binned_sequential(palette = palette[x]) +
    ylab("")
  
  return(p)
})

names(gglist) <- paste("C", 1:5, sep = "")

# Fig.03 ---------------
pdf('scripts/Pseudotime/results/figure/Trajectory_Heatmap_with_GO.pdf', height = 10, width = 14, onefile = FALSE)
visCluster(object = ck,
           add.sampleanno = FALSE,
           plot.type = "both",
           cluster_rows = FALSE,
           cluster.order = c(1, 2, 5, 4, 3),
           markGenes = sample(rownames(mat), 30, replace = FALSE),
           ht.col.list = list(
             col_range = seq(-2, 2, length.out = 10), 
             col_color = my_base_colors
           ),
           line.side = "left",
           ggplot.panel.arg = c(4.5, 0.5, 16, "grey90", NA),
           gglist = gglist)
dev.off()

# Fig.04 ---------------
pt_values <- pseudotime(cds)
pt_values <- pt_values[is.finite(pt_values)]
ordered_barcodes <- names(sort(pt_values))
sorted_cell_types <- cds@colData[ordered_barcodes, "CellType1"]

df_bar <- data.frame(
  pos = 1:length(sorted_cell_types), 
  y = 1,
  CellType = sorted_cell_types
) 

p_bar <- ggplot(df_bar, aes(x = pos, y = y, fill = CellType)) +
  geom_raster() + 
  scale_fill_manual(values = my_anno_colors) + 
  theme_void() + 
  theme(
    legend.position = "none", 
    plot.margin = margin(0, 0, 0, 0) 
  )

ggsave("scripts/Pseudotime/results/figure/cell_type_bar_only.pdf", plot = p_bar, width = 8, height = 0.4, limitsize = FALSE)

# Source: source/project/Reference/combine/scripts/Pseudotime/scripts/module_pseudotime.r
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=50)
# Console中运行trace('calculateLW', edit = T, where = asNamespace("monocle3"))
# 在出现的文档中Matrix::rBind 改为rbind，或者将rbind 改为Matrix::rBind

#导出
write.csv(Track_genes, "scripts/Pseudotime/results/table/Track_genes_all.csv", row.names = F)
#保存
qsave(Track_genes,'scripts/Pseudotime/data/processed_data/Track_genes.qs')
Track_genes <- qread('scripts/Pseudotime/data/processed_data/Track_genes.qs')

#拟时基因热图
genes <- row.names(subset(Track_genes, q_value == 0 & morans_I > 0.15))

mat <- pre_pseudotime_matrix(cds_obj = cds,
                             gene_list = genes)

#保存
qsave(mat,'scripts/Pseudotime/data/processed_data/mat.qs')
mat <- qread('scripts/Pseudotime/data/processed_data/mat.qs')

head(mat[1:5,1:5])


# # -----------
cds <- qread('scripts/Pseudotime/data/processed_data/save_cds.qs')
mat <- qread('scripts/Pseudotime/data/processed_data/mat.qs')
# 
# #选取合适的clusters
# 
# library(ClusterGVis)
# 
# #kmeans
# ck <- clusterData(exp = mat,
#                   cluster.method = "kmeans",
#                   cluster.num = 5)
# 
# qsave(ck,'scripts/Pseudotime/data/processed_data/ck.qs')
ck <- qread('scripts/Pseudotime/data/processed_data/ck.qs')



my_base_colors <- rev(c("#810000", "#CE1212", "#F05454", "#ffd06f",
                        "#ffe6b7", "#aadce0", "#72bcd5", "#528fad",
                        "#376795", "#1e466e"))


heatmap_cells <- colnames(mat)
my_cell_types <- cds@colData[heatmap_cells, "CellType1"]
my_cell_types <- as.character(my_cell_types)
print(length(my_cell_types) == ncol(mat))
unique_types <- unique(my_cell_types)

my_anno_colors <- c(
  "MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A",
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628"
)

names(my_anno_colors) <- c("MG2", "MG5", "MG0", "MG3", "MG4", "MG1")
my_expanded_cell_colors <- my_anno_colors[my_cell_types]
print(length(my_expanded_cell_colors) == 21244)

# add line annotation
pdf('scripts/Pseudotime/results/figure/monocle3_2.pdf',height = 10,width = 8,onefile = F)
visCluster(object = ck,
           add.sampleanno = F,
           plot.type = "both",
           cluster_rows = FALSE,
           cluster.order = c(1, 2, 5, 4, 3),
           markGenes = sample(rownames(mat),30,replace = F),
           ht.col.list = list(
             col_range = seq(-2, 2, length.out = 10), 
             col_color = my_base_colors
           ),
           line.side = "left",
           )
dev.off()

# 带GO的热图--------------
enrich <- enrichCluster(object = ck,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 100,
                        seed = 24)

# 获取所有不重复的 cluster (group) 名称
groups <- unique(enrich$group)

# 循环遍历每一个 cluster
for (g in groups) {
  # 提取当前 cluster 的数据
  tmp_data <- enrich |> dplyr::filter(group == g)
  
  # 设置要保存的文件名，例如 "GO_Result_C1.csv"
  file_name <- paste0(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/Pseudotime/results/table/GO_Result_"), g, ".csv")
  
  # 导出为 CSV 文件
  write.csv(tmp_data, file = file_name, row.names = FALSE)
  
  cat("已保存:", file_name, "\n")
}
pdf('scripts/Pseudotime/results/figure/Trajectory_Heatmap_with_GO.pdf', height = 10, width = 14, onefile = F)

visCluster(object = ck,
           add.sampleanno = F,
           plot.type = "both",
           cluster_rows = FALSE,
           cluster.order = c(1, 2, 5, 4, 3),
           markGenes = sample(rownames(mat), 30, replace = F),
           ht.col.list = list(
             col_range = seq(-2, 2, length.out = 10), 
             col_color = my_base_colors
           ),
           line.side = "left",
           ggplot.panel.arg = c(4.5,0.5,16,"grey90",NA),
           gglist = gglist
)
dev.off()


# 手动画注释的bar----------------------
library(ggplot2)

# 1. 提取所有细胞的拟时间值
pt_values <- pseudotime(cds)

# 2. 过滤掉没有算在轨迹里的细胞（拟时间为 Infinite 的细胞）
# 通常 pre_pseudotime_matrix 内部也会自动丢弃这些细胞
pt_values <- pt_values[is.finite(pt_values)]

# 3. 【核心步骤】按照拟时间从小到大，对细胞条形码（Barcode）进行排序
ordered_barcodes <- names(sort(pt_values))

# 建议运行这行检查一下，应该是 TRUE
print(paste("排序后的细胞数与 mat 列数是否一致:", length(ordered_barcodes) == ncol(mat)))

# 4. 按照这个排序好的细胞顺序，去提取对应的细胞类型
# 这样提取出来的细胞类型，就和 mat 里的列 100% 一一对应了！
sorted_cell_types <- cds@colData[ordered_barcodes, "CellType1"]

# 5. 构建绘图数据框
df_bar <- data.frame(
  pos = 1:length(sorted_cell_types), # x 轴 1 到两万多
  y = 1,
  # 加上 "MG" 前缀，并转化为 factor，匹配你的颜色字典
  CellType =sorted_cell_types) 

# 6. 画图
p_bar <- ggplot(df_bar, aes(x = pos, y = y, fill = CellType)) +
  geom_raster() + # geom_raster 处理两万个点极快
  scale_fill_manual(values = my_anno_colors) + 
  theme_void() + 
  theme(
    legend.position = "none", 
    plot.margin = margin(0, 0, 0, 0) 
  )

# 7. 保存为又长又扁的 PDF
ggsave("scripts/Pseudotime/results/figure/cell_type_bar_only.pdf", 
       plot = p_bar, 
       width = 8,     # 宽度建议和你的热图保持一致
       height = 0.4, 
       limitsize = FALSE)

# F5D  TRIM2, RNF220, MBP, MOBP, ST18 and PLP1 along pseudotime
# Source: source/project/code_upload/Pseudotime.r
# Fig.05 ---------------
genes_to_plot <- c('TRIM2','RNF220','ST18','MBP','BCAS1','MOBP','PLP1') 

p_jitter <- monocle3::plot_genes_in_pseudotime(
  cds[genes_to_plot, ], 
  color_cells_by = "CellType1", 
  min_expr = 0.5, 
  ncol = 7,
  panel_order = genes_to_plot 
) + 
  scale_colour_manual(values = my_anno_colors)

ggsave("scripts/Pseudotime/results/figure/Genes_Jitterplot.pdf", plot = p_jitter, width = 16, height = 2)

# Source: source/project/Reference/combine/scripts/Pseudotime/scripts/module_pseudotime.r
# 指定基因画jitter plot-----------
# genes <- c('TRIM2','RNF220','ST18','MBP','BCAS1','MOBP','PLP1')
genes <- c('ACSL1','CTSL','MOG','CTNNA2','FGF14')

my_anno_colors <- c(
  "MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A",
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628"
)

p <- monocle3::plot_genes_in_pseudotime(
  cds[genes,], 
  color_cells_by = "CellType1", 
  min_expr = 0.5, 
  ncol = 7,
  panel_order = genes  # <--- 核心参数：直接指定分面顺序
) + 
  scale_colour_manual(values = my_anno_colors)

ggsave("scripts/Pseudotime/results/figure/Genes_Jitterplot_2.pdf", plot = p, width = 16, height = 2)

# F5E  Mouse-myelin debris stimulation of HMC3 experimental schematic

# F5F  Iba1/Trim2/MOG immunofluorescence in control versus myelin-exposed HMC3

# F5G  TRIM2/MBP/MOG qPCR treatment-versus-control barplots

# F5H  OLM gene heatmap in snRNA-seq and six-sample bulk RNA-seq
# Source: source/project/Reference/combine/scripts/Heatmap.r
# MODULE 9: Microglia Other Function Gene Heatmap ------------------------------
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/"))
rm(list=ls())

## 9.1 Parameters --------------------------------------------------------------
input_file      <- 'data/processed_data/mg_cca_meta_v2.qs'
output_fig_dir  <- 'results/figure/'
output_filename <- paste0(output_fig_dir, 'bulk_gene_heatmap.pdf')

dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

target_genes <- c("PLLP", "FAM13C", "ITGA2", "NCKAP5", "TTLL11", 
                  "ZNF608", "AK5", "CORO2B", "TRIM9", "DIXDC1", 
                  "GAB1", "ARNT2", "KIAA1958", "TNRC6C", "RNF144A", 
                  "SH3TC2", "TRPC1", "PLEKHB1", "TTLL7", "TRIM2", 
                  "ENOSF1", "PRKD1", "LRRC1", "SEMA4D", "FAXDC2", 
                  "SLC13A3")


## 9.2 Process Data ------------------------------------------------------------
library(Seurat)
library(pheatmap)
library(qs)

if(!file.exists(input_file)) stop("Input file not found.")
tmp_obj <- qread(input_file)

# Calculate Average Expression
avg_exp_list <- AverageExpression(
  tmp_obj,
  group.by = "CellType1",
  assays = "RNA"
)
avg_mat <- as.data.frame(avg_exp_list$RNA)

# Filter and Sort Genes
# Ensure we only define genes that exist in the matrix and maintain input order
valid_genes <- intersect(target_genes, rownames(avg_mat))
plot_data   <- avg_mat[valid_genes, , drop = FALSE]
plot_data   <- plot_data[match(valid_genes, rownames(plot_data)), ]


## 9.3 Generate & Save Plot ----------------------------------------------------
pdf(output_filename, width = 3, height = 3.5)

pheatmap(
  plot_data,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  show_rownames = TRUE,
  legend = TRUE
)

dev.off()
# MODULE 9: Microglia Other Function Gene Heatmap ------------------------------
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/"))
rm(list=ls())

## 9.1 Parameters --------------------------------------------------------------
output_fig_dir  <- 'results/figure/'
output_filename <- paste0(output_fig_dir, 'bulk_heatmap.pdf')

dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

target_genes <- c("PLLP", "FAM13C", "ITGA2", "NCKAP5", "TTLL11", 
                  "ZNF608", "AK5", "CORO2B", "TRIM9", "DIXDC1", 
                  "GAB1", "ARNT2", "KIAA1958", "TNRC6C", "RNF144A", 
                  "SH3TC2", "TRPC1", "PLEKHB1", "TTLL7", "TRIM2", 
                  "ENOSF1", "PRKD1", "LRRC1", "SEMA4D", "FAXDC2", 
                  "SLC13A3")


## 9.2 Process Data ------------------------------------------------------------
library(Seurat)
library(pheatmap)
library(qs)

vsd <- readRDS(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/processed_data/vsd.rds"))
mat <- assay(vsd)[target_genes, ]

# Filter and Sort Genes
# Ensure we only define genes that exist in the matrix and maintain input order
valid_genes <- intersect(target_genes, rownames(mat))
plot_data   <- mat[valid_genes, , drop = FALSE]
plot_data   <- plot_data[match(valid_genes, rownames(plot_data)), ]


## 9.3 Generate & Save Plot ----------------------------------------------------
pdf(output_filename, width = 3, height = 3.5)

pheatmap(
  plot_data,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  show_rownames = TRUE,
  legend = TRUE
)

dev.off()

message("Module 9 Done.")

# Source: source/project/Reference/combine/scripts/BULK/scripts/DESeq2.r
# 1. 数据加载---------------------
raw_counts <- read.csv(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/processed_data/All_gene_counts.list"),sep = '\t',header=TRUE,row.names=1)
raw_counts <- raw_counts[,c(-4,-8)]

# 4. 数据预处理与质控----------------
col_data <- read.csv(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/saved_data/sample_info.txt"),sep = '\t',header=TRUE,row.names = 1)
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = col_data,
  design = ~ group)

# 过滤低表达基因
keep <- rowSums(counts(dds) >= 10) >= 4 # 我们的过滤标准（至少4个样本中计数≥10）在统计效能和基因保留之间取得了良好平衡
dds <- dds[keep,]
cat("过滤后保留", nrow(dds), "个基因\n")

# 方差稳定变换
vsd <- vst(dds, blind = FALSE) # vst（方差稳定变换）：计算速度更快，特别适用于大样本量。rlog（正则化对数变换）：计算速度较慢，但对于小样本量可能表现稍好。
saveRDS(vsd,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/processed_data/vsd.rds"))

# F5I  GO barplot of five terms from preselected workbook
# Source: source/project/Reference/combine/scripts/20260519_go_barplot.r
library(readxl)
library(ggplot2)
library(dplyr)

go_data <- read_excel(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/saved_data/20260519_go.xlsx"))

p <- ggplot(go_data, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_col(width = 0.65) +
  scale_fill_gradient(low = "#2E86AB", high = "#A2BCE0") +
  labs(title = "GO Enrichment Bar Plot", x = "Gene Count", y = "GO Term", fill = "Adjusted p-value") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10), panel.grid = element_blank())

ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/results/figure/20260519_go_barplot.pdf"), plot = p, width = 8, height = 4, dpi = 300)
ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/results/figure/20260519_go_barplot.png"), plot = p, width = 8, height = 4, dpi = 300)


# Supplementary Figure 7 ------------------------------------------------------------

# S7A  Three control HMC3 fields: Iba1/Trim2/MOG/DAPI

# S7B  Three myelin-exposed HMC3 fields: Iba1/Trim2/MOG/DAPI


# Figure 6 ------------------------------------------------------------

# F6A  GFP-HMC3 myelin conditioning and corpus-callosum transplantation schematic

# F6B  GFP-HMC3 fluorescence and viral transduction violin plot
# Source: source/project/Reference/combine/scripts/Boxplot/scripts/sheet5_vlnplot.r
library(openxlsx)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rstatix)

# MODULE : Group Comparisons Violin Plot -------------------

## 1 Parameters  ---------------------------------------------------
input_filepath  <- "data/processed_data/sheet5.xlsx"
plot_name       <- "sheet5_violin"
output_filepath <- paste0("results/figure/", plot_name, ".pdf")

var_group_name  <- "Group"
var_value_name  <- "Value"

mycolor      <- c("#0077c1", "#00a99e", "#6bc72b")
label_x      <- ""
label_y      <- "value"
plot_margins <- c(0.2, 0.2, 0.2, 0.2) 

stat_method    <- "t.test"
my_symnum      <- list(
  cutpoints = c(0, 0.001, 0.01, 0.05, 1),
  symbols   = c("***", "**", "*", "ns")
)

# 【修改1：调整画板宽高】竖着的图需要更高一点的空间放显著性括号
out_width  <- 2   
out_height <- 5   

## 2 Process Data ------------------------------------------------------------
df <- read.xlsx(input_filepath)

dt <- melt(df,
           id.vars       = NULL,
           measure.vars  = colnames(df),
           variable.name = var_group_name,
           value.name    = var_value_name)

# ### 2.5 Calculate Statistics  ---------------------------------------
# stat.test <- dt %>%
#   t_test(Value ~ Group) %>%
#   add_xy_position(x = "Group")
# 
# stat.test$p.signif <- symnum(stat.test$p, 
#                              corr = FALSE, na = FALSE, 
#                              cutpoints = my_symnum$cutpoints, 
#                              symbols = my_symnum$symbols)


## 3 Generate Plots ----------------------------------------------------------
mytheme2 <- theme_classic() +
  theme(
    text              = element_text(family = "sans", colour = "black", size = 12),
    axis.line         = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks.length = unit(1.5, units = "mm"),
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
    legend.position   = "none",
    plot.margin       = unit(x = plot_margins, units = "inches")
  )

p_base <- ggviolin(dt, 
                   x     = var_group_name, 
                   y     = var_value_name,
                   color = var_group_name,
                   width = 0.6,               # 小提琴图稍微放宽一点更好看
                   fill  = var_group_name,
                   add   = "jitter",          # 依然保留散点
                   add.params = list(
                     shape = 21,            # 设置点的形状为 21 号（可描边的圆）
                     size  = 1.5,           # 点的大小
                     color = "black",       # 点的描边颜色（黑色）
                     fill  = var_group_name # 点的填充颜色（映射到分组，与小提琴图一致）
                   )) + 
  scale_colour_manual(values = alpha(mycolor, 1)) +
  scale_fill_manual(values = alpha(mycolor, 0.3)) +
  labs(x = label_x, y = label_y) +
  mytheme2


p_final <- p_base

# # 【修改3：清理掉用于翻转的参数】
# p_final <- p_base + stat_pvalue_manual(
#   stat.test,
#   label        = "p.signif",
#   tip.length   = 0.02,
#   bracket.size = 0.5
#   # 删除了 angle = -90 和 coord.flip = TRUE
# )

print(p_final)

## 4 Save Results ------------------------------------------------------------
out_dir <- dirname(output_filepath)
if(!dir.exists(out_dir) && out_dir != ".") {
  dir.create(out_dir, recursive = TRUE)
}

ggsave(filename = output_filepath, 
       plot     = p_final, 
       width    = out_width, 
       height   = out_height)

# F6C  Olig2/GFP/DAPI corpus-callosum images after HMC3-control or myelin-exposed HMC3 graft

# F6D  Olig2-positive OLG fluorescent area fraction after graft
# Source: source/recovered/experimental/F6D_OLG_area_fraction.R
library(openxlsx)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rstatix)
input_filepath <- file.path(Sys.getenv("MG_PLOT_INPUT_ROOT", "data/plot_inputs"), "sp8_m0.xlsx")
plot_name       <- "sp8_m0"
output_filepath <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "F6D", paste0(plot_name, ".pdf"))
var_group_name  <- "Group"
var_value_name  <- "Value"
mycolor      <- c("#0077c1", "#00a99e", "#6bc72b")
label_x      <- ""
label_y      <- "value"
plot_margins <- c(0.2, 0.2, 0.2, 0.2) 
stat_method    <- "t.test"
my_symnum      <- list(
  cutpoints = c(0, 0.001, 0.01, 0.05, 1),
  symbols   = c("***", "**", "*", "ns")
)
out_width  <- 2   
out_height <- 3   
df <- read.xlsx(input_filepath)
dt <- melt(df,
           id.vars       = NULL,
           measure.vars  = colnames(df),
           variable.name = var_group_name,
           value.name    = var_value_name)
stat.test <- dt %>%
  t_test(Value ~ Group) %>%
  add_xy_position(x = "Group")
stat.test$p.signif <- symnum(stat.test$p,
                             corr = FALSE, na = FALSE,
                             cutpoints = my_symnum$cutpoints,
                             symbols = my_symnum$symbols)
mytheme2 <- theme_classic() +
  theme(
    text              = element_text(family = "sans", colour = "black", size = 12),
    axis.line         = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks.length = unit(1.5, units = "mm"),
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
    legend.position   = "none",
    plot.margin       = unit(x = plot_margins, units = "inches")
  )
p_base <- ggboxplot(dt, 
                    x     = var_group_name, 
                    y     = var_value_name,
                    color = var_group_name,
                    width = 0.4,
                    fill  = var_group_name,
                    add   = "jitter",
                    add.params = list(size = 0.8)) + 
  scale_colour_manual(values = alpha(mycolor, 1)) +
  scale_fill_manual(values = alpha(mycolor, 0.3)) +
  labs(x = label_x, y = label_y) +
  # coord_flip() + # <--- 使用 coord_flip() 翻转整个画布
  mytheme2
p_final <- p_base + stat_pvalue_manual(
  stat.test,
  label        = "p.signif",
  tip.length   = 0.02,
  bracket.size = 0.5
  # 删除了 angle = -90 和 coord.flip = TRUE
)
print(p_final)
out_dir <- dirname(output_filepath)
if(!dir.exists(out_dir) && out_dir != ".") {
  dir.create(out_dir, recursive = TRUE)
}
ggsave(filename = output_filepath, 
       plot     = p_final, 
       width    = out_width, 
       height   = out_height)

# F6E  Young C57/young SAMP8/aging SAMP8 schematic and TEM images

# F6F  G-ratios comparing young C57, young SAMP8 and aging SAMP8
# Source: source/recovered/experimental/F6F_G_ratio_age_groups.R
library(openxlsx)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rstatix)
input_filepath <- file.path(Sys.getenv("MG_PLOT_INPUT_ROOT", "data/plot_inputs"), "sheet1.xlsx")
plot_name       <- "sheet1"
output_filepath <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "F6F", paste0(plot_name, ".pdf"))
var_group_name  <- "Group"
var_value_name  <- "Value"
mycolor      <- c("#0077c1", "#00a99e", "#6bc72b")
label_x      <- ""
label_y      <- "value"
plot_margins <- c(0.2, 0.2, 0.2, 0.2) 
stat_method    <- "t.test"
my_symnum      <- list(
  cutpoints = c(0, 0.001, 0.01, 0.05, 1),
  symbols   = c("***", "**", "*", "ns")
)
out_width  <- 6
out_height <- 2
df <- read.xlsx(input_filepath)
dt <- melt(df,
           id.vars       = NULL,
           measure.vars  = colnames(df),
           variable.name = var_group_name,
           value.name    = var_value_name)
stat.test <- dt %>%
  t_test(Value ~ Group, 
         # comparisons = my_comparisons
  ) %>%
  add_xy_position(x = "Group")
stat.test$p.signif <- symnum(stat.test$p, 
                             corr = FALSE, na = FALSE, 
                             cutpoints = my_symnum$cutpoints, 
                             symbols = my_symnum$symbols)
mytheme2 <- theme_classic() +
  theme(
    text              = element_text(family = "sans", colour = "black", size = 12),
    axis.line         = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks.length = unit(1.5, units = "mm"),
    legend.position   = "none",
    plot.margin       = unit(x = plot_margins, units = "inches")
  )
p_base <- ggboxplot(dt, 
                    x     = var_group_name, 
                    y     = var_value_name,
                    color = var_group_name,
                    width = 0.4,
                    fill  = var_group_name,
                    add   = "jitter",
                    add.params = list(size = 0.8)) + 
  scale_colour_manual(values = alpha(mycolor, 1)) +
  scale_fill_manual(values = alpha(mycolor, 0.3)) +
  labs(x = label_x, y = label_y) +
  coord_flip() + # <--- 使用 coord_flip() 翻转整个画布
  mytheme2
p_final <- p_base + stat_pvalue_manual(
  stat.test,
  label        = "p.signif",
  tip.length   = 0.02,
  bracket.size = 0.5,
  angle        = -90,  # 让星号顺时针旋转90度，与括号平行
  coord.flip   = TRUE  # <--- 关键参数：告诉它画布被翻转了，星号和括号会自动对齐
)
print(p_final)
out_dir <- dirname(output_filepath)
if(!dir.exists(out_dir) && out_dir != ".") {
  dir.create(out_dir, recursive = TRUE)
}
ggsave(filename = output_filepath, 
       plot     = p_final, 
       width    = out_width, 
       height   = out_height)

# F6G  OLM transplantation/remyelination schematic and CC-saline/CC-HMC3/CC-OLM TEM

# F6H  G-ratios for control, HMC3 and OLM transplantation groups
# Source: source/recovered/experimental/F6H_G_ratio_transplant.R
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggbeeswarm)
  library(patchwork)
})
input_file <- file.path(Sys.getenv("MG_PLOT_INPUT_ROOT", "data/plot_inputs"), "2026-7-29-统计图1.xlsx")
output_dir <- file.path(Sys.getenv("MG_OUTPUT_ROOT", "outputs"), "F6H")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
raw_g <- read_excel(input_file, sheet = "Sheet1", range = "A1:J146", col_names = FALSE)
group_cols <- c(3, 7, 10)
group_names <- as.character(unlist(raw_g[6, group_cols]))
g_ratio <- lapply(seq_along(group_cols), function(i) {
  tibble(
    group = group_names[i],
    value = suppressWarnings(as.numeric(unlist(raw_g[8:nrow(raw_g), group_cols[i]])))
  )
}) |>
  bind_rows() |>
  filter(!is.na(value)) |>
  mutate(group = factor(group, levels = group_names))
pair_index <- combn(group_names, 2, simplify = FALSE)
pairwise_tests <- lapply(pair_index, function(z) {
  test <- t.test(value ~ group, data = filter(g_ratio, group %in% z))
  tibble(group1 = z[1], group2 = z[2], p_value = test$p.value)
}) |>
  bind_rows() |>
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adjusted <= 0.001 ~ "***",
      p_adjusted <= 0.01 ~ "**",
      p_adjusted <= 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    significance_raw = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
print(pairwise_tests |> select(group1, group2, p_value, p_adjusted))
summary_g <- g_ratio |>
  group_by(group) |>
  summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    q1 = quantile(value, 0.25),
    q3 = quantile(value, 0.75),
    .groups = "drop"
  )
palette_g <- setNames(c("#0077c1", "#00a99e", "#6bc72b"), group_names)
set.seed(20260729)
p_g <- ggplot(g_ratio, aes(x = value, y = group, color = group)) +
  geom_boxplot(
    aes(fill = group),
    width = 0.40, outlier.shape = NA,
    linewidth = 0.60, staplewidth = 0
  ) +
  geom_point(
    size = 0.80,
    position = position_jitter(
      width = 0, height = 0.15, seed = 20260729
    )
  ) +
  scale_color_manual(values = palette_g, guide = "none") +
  scale_fill_manual(
    values = scales::alpha(palette_g, 0.30),
    guide = "none"
  ) +
  scale_x_continuous(
    name = "G-ratio",
    limits = c(0.48, 1.025),
    breaks = seq(0.5, 0.9, 0.1),
    expand = expansion(mult = c(0.01, 0))
  ) +
  labs(y = NULL) +
  theme_classic() +
  theme(
    text = element_text(
      family = "sans", colour = "black", size = 12
    ),
    axis.line = element_line(linewidth = 0.60, colour = "black"),
    axis.ticks = element_line(linewidth = 0.60, colour = "black"),
    axis.ticks.length = grid::unit(1.5, units = "mm"),
    legend.position = "none",
    plot.margin = grid::unit(
      c(0.2, 0.2, 0.2, 0.2), units = "inches"
    )
  ) +
  coord_cartesian(clip = "off")
bracket_x <- c(0.972, 0.998, 0.950)
add_comparison_brackets <- function(plot_object, labels) {
  for (i in seq_len(nrow(pairwise_tests))) {
    y1 <- match(pairwise_tests$group1[i], group_names)
    y2 <- match(pairwise_tests$group2[i], group_names)
    x <- bracket_x[i]
    plot_object <- plot_object +
      annotate("segment", x = x, xend = x, y = y1, yend = y2, linewidth = 0.5) +
      annotate("segment", x = x - 0.007, xend = x, y = y1, yend = y1, linewidth = 0.5) +
      annotate("segment", x = x - 0.007, xend = x, y = y2, yend = y2, linewidth = 0.5) +
      annotate(
        "text", x = x + 0.009, y = (y1 + y2) / 2,
        label = labels[i], angle = -90,
        size = 3.5
      )
  }
  plot_object
}
p_g_raw <- add_comparison_brackets(p_g, pairwise_tests$significance_raw)
p_g <- add_comparison_brackets(p_g, pairwise_tests$significance)
ggsave(file.path(output_dir, "图1_G-ratio_散点箱线图.png"), p_g, width = 6, height = 2, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "图1_G-ratio_散点箱线图.pdf"), p_g, width = 6, height = 2, device = cairo_pdf)
ggsave(file.path(output_dir, "图1B_G-ratio_未校正P值_散点箱线图.png"), p_g_raw, width = 6, height = 2, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "图1B_G-ratio_未校正P值_散点箱线图.pdf"), p_g_raw, width = 6, height = 2, device = cairo_pdf)
write.csv(summary_g, file.path(output_dir, "G-ratio_描述统计.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(pairwise_tests, file.path(output_dir, "G-ratio_Welch_t检验_BH校正.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(summary_g)
print(pairwise_tests)


# Supplementary Figure 8 ------------------------------------------------------------

# S8A  Three representative TEM fields each for CC-saline, CC-HMC3 and CC-HMC3-mye

