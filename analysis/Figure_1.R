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


# Figure 1

# F1A  Study-design schematic and embedded human microglia UMAP
# Kind: mixed_manual_computational
mg_register_panel("F1A", list(
  list(source="source/project/code_upload/Fig.01.r", ranges="77-107", requested_ranges="77-107", role="", code=quote({
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
  }))
))

# F1B  Nuclei and genes per nucleus by Condition
# Kind: computational
mg_register_panel("F1B", list(
  list(source="source/figure_settings/F1B_Condition.R", ranges="1-46", requested_ranges="1-46", role="Original plotting expressions with the figure-specific grouping", code=quote({
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

  }))
))

# F1C  Six-study dataset summary table
# Kind: manual_summary_table
mg_register_panel("F1C", NULL)

# F1D  Six functional AddModuleScore UMAPs
# Kind: computational
mg_register_panel("F1D", list(
  list(source="source/project/code_upload/Fig.01.r", ranges="545-643", requested_ranges="545-643", role="", code=quote({
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
  }))
))

# F1E_dot  36 selected genes (six per MG state), dot plot component
# Kind: computational
mg_register_panel("F1E_dot", list(
  list(source="source/project/Reference/combine/scripts/20260519_dotplot.r", ranges="24-65", requested_ranges="24-65", role="", code=quote({
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
  }))
))

# F1E_heatmap  Top50 marker heatmap component
# Kind: computational
mg_register_panel("F1E_heatmap", list(
  list(source="source/project/code_upload/Fig.01.r", ranges="193-234", requested_ranges="193-234", role="", code=quote({
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
  }))
))

# F1E_GO  MG-state GO text blocks
# Kind: computational
mg_register_panel("F1E_GO", list(
  list(source="source/project/code_upload/Fig.01.r", ranges="236-428", requested_ranges="236-428", role="", code=quote({
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
  }))
))

# F1E_SCENIC  Regulon activity heatmap component
# Kind: computational
mg_register_panel("F1E_SCENIC", list(
  list(source="source/project/code_upload/SCENIC.r", ranges="1-134", requested_ranges="1-134", role="", code=quote({
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
  }))
))
