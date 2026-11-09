setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/"))
rm(list=ls())

# 0. Global Setup --------------------------------------------------------------
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


## 重新分配cluster名称----------
tmp_obj <- qread('data/processed_data/mg_cca_meta.qs') 

tmp_obj@meta.data <- tmp_obj@meta.data %>%
  mutate(CellType1 = case_when(
    integrated_snn_res.0.3 %in% c("4", "6") ~ "MG0", # 髓鞘相关 (Myelin-associated)
    integrated_snn_res.0.3 %in% c("8", "9") ~ "MG1", # 突触相关 (Synaptic-associated)
    integrated_snn_res.0.3 %in% c("0", "2") ~ "MG2", # 稳态 (Homeostatic)
    integrated_snn_res.0.3 %in% c("5")      ~ "MG3", # 疾病相关 (Disease-associated)
    integrated_snn_res.0.3 %in% c("1","7")  ~ "MG4", # 炎症相关 (Inflammatory-associated)
    integrated_snn_res.0.3 %in% c("3")      ~ "MG5", # 增殖相关 (Proliferation-associated)

    TRUE ~ integrated_snn_res.0.3 
  ))
tmp_obj$CellType1 <- as.factor(tmp_obj$CellType1)
unique(tmp_obj$CellType1)
qsave(tmp_obj,'data/processed_data/mg_cca_meta_v2.qs')

# MODULE: DimPlot---------------------
cluster_col <- 'CellType1'
manual_colors <- c(
  "MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A", 
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628"
)
(p1 <- DimPlot(tmp_obj,group.by = cluster_col,cols = manual_colors,label = T))
(p2 <- DimPlot(tmp_obj,group.by = cluster_col,cols = manual_colors))
ggsave(paste0('results/figure/',cluster_col,'_dimplot.pdf'),plot = p1,width = 6,height = 5)
ggsave(paste0('results/figure/',cluster_col,'_dimplot.png'),plot = p1,width = 6,height = 5)

# MODULE 1: FindAllMarkers (Calculation) ---------------------------------------

## 1.1 Parameters --------------------------------------------------------------
cluster_col <- 'CellType1'      
min_pct     <- 0.25                          
logfc_thr   <- 0.25                          

# Auto-filenames
res_suffix  <- gsub("integrated_snn_", "", cluster_col)
qs_path     <- paste0('data/saved_data/', res_suffix, '_markers.qs')
xlsx_path   <- paste0('results/table/', res_suffix, '_markers.xlsx')

## 1.2 Run Analysis ------------------------------------------------------------
tmp_obj <- qread('data/processed_data/mg_cca_meta_v2.qs') 
Idents(tmp_obj) <- cluster_col

message(paste("Calculating Markers for:", cluster_col))
all_markers <- FindAllMarkers(tmp_obj, only.pos = TRUE, min.pct = min_pct, logfc.threshold = logfc_thr)

## 1.3 Save --------------------------------------------------------------------
qsave(all_markers, qs_path)
write.xlsx(all_markers, file = xlsx_path, overwrite = TRUE)
message("Module 1 Done.")



############# ********************* #####################
# MODULE 2: DotPlot (Visualization) --------------------------------------------

## 2.1 Parameters --------------------------------------------------------------
cluster_col <- 'CellType1'      
top_n       <- 5                            

# Auto-filenames
res_suffix  <- gsub("integrated_snn_", "", cluster_col)
marker_path <- paste0('data/saved_data/', res_suffix, '_markers.qs') 
png_path    <- paste0("results/figure/", res_suffix, "_top", top_n, "_dotplot.png")
pdf_path    <- paste0("results/figure/", res_suffix, "_top", top_n, "_dotplot.pdf")

## 2.2 Process Data ------------------------------------------------------------
tmp_obj <- qread('data/processed_data/mg_cca_meta_v2.qs')
if(!file.exists(marker_path)) stop("Marker file not found. Run Module 1 first.")
markers <- qread(marker_path)

top_markers <- markers %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup()

genes_list <- unique(top_markers$gene)

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
  labs(title = paste0("Top ", top_n, " Genes in ", res_suffix), x = NULL, y = NULL) +
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



############# ********************* #####################
# MODULE 3: Heatmap (Visualization) --------------------------------------------

## 3.1 Parameters --------------------------------------------------------------
cluster_col <- 'CellType1'      
top_n       <- 50                            

# Auto-filenames
res_suffix  <- gsub("integrated_snn_", "", cluster_col)
marker_path <- paste0('data/saved_data/', res_suffix, '_markers.qs')
png_path    <- paste0("results/figure/", res_suffix, "_top", top_n, "_heatmap.png")
pdf_path    <- paste0("results/figure/", res_suffix, "_top", top_n, "_heatmap.pdf")

## 3.2 Process Data ------------------------------------------------------------
tmp_obj <- qread('data/processed_data/mg_cca_meta_v2.qs')
if(!file.exists(marker_path)) stop("Marker file not found. Run Module 1 first.")
markers <- qread(marker_path)

top_markers <- markers %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup()

DefaultAssay(tmp_obj) <- 'RNA'
avg_exp_obj <- AverageExpression(tmp_obj, group.by = cluster_col, assays = "RNA")
heatmap_data <- avg_exp_obj$RNA[top_markers$gene, ]

## 3.3 Plotting ----------------------------------------------------------------
p_heat_obj <- pheatmap(heatmap_data, 
                       cluster_cols = FALSE, cluster_rows = FALSE,
                       scale = "row", show_rownames = FALSE, legend = T, 
                       silent = TRUE, angle_col = "0",border_color = NA)
p_heat_gg <- as.ggplot(p_heat_obj)

## 3.4 Save --------------------------------------------------------------------
ggsave(png_path, p_heat_gg, width = 5, height = 18, limitsize = FALSE)
ggsave(pdf_path, p_heat_gg, width = 5, height = 18, limitsize = FALSE)
message("Module 3 Done.")



############# ********************* #####################
# MODULE 4: GO Enrichment (Functional Analysis) --------------------------------

## 4.1 Parameters --------------------------------------------------------------
cluster_col <- 'CellType1'       
top_n_go    <- 10                           

# Manual Colors (Modify based on your real cluster names)
manual_colors <- c(
  "MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A", 
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628"
)

# Auto-filenames
res_suffix   <- gsub("integrated_snn_", "", cluster_col)
marker_path  <- paste0('data/saved_data/', res_suffix, '_markers.qs')
go_save_path <- paste0("data/saved_data/", res_suffix, "_res_go_list.qs")
excel_path   <- paste0('results/table/', res_suffix, '_GO_enrichment.xlsx')
merged_pdf   <- paste0("results/figure/", res_suffix, "_GO_Combined_Text.pdf")
bar_plot_dir <- paste0("results/figure/", res_suffix, "_individual_GO_barplots/")
dir.create(bar_plot_dir, recursive = TRUE, showWarnings = FALSE)


## 4.2 Process Data (Run GO & Assign Colors) -----------------------------------
library(RColorBrewer)
library(scales)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(openxlsx)
library(qs)
library(dplyr)
library(ggplot2)
library(patchwork)

if(!file.exists(marker_path)) stop("Marker file not found. Run Module 1 first.")
markers <- qread(marker_path)

# 1. Define Cluster List
all_clusters <- as.character(sort(unique(markers$cluster)))

# 2. Assign Colors
final_cols <- manual_colors[all_clusters]
if (any(is.na(final_cols))) {
  n_miss <- length(all_clusters)
  auto_cols <- if(n_miss <= 12) brewer.pal(max(n_miss, 3), "Set3")[1:n_miss] else colorRampPalette(brewer.pal(12, "Set3"))(n_miss)
  final_cols <- setNames(auto_cols, all_clusters)
  for(cls in names(manual_colors)) if(cls %in% all_clusters) final_cols[cls] <- manual_colors[cls]
}
cluster_colors <- final_cols[!is.na(final_cols)]

# 3. GO Calculation
if (file.exists(go_save_path)) {
  res_go <- qread(go_save_path)
} else {
  res_go <- list()
  for (i in all_clusters) {
    genes <- markers$gene[markers$cluster == i & markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.5]
    if(length(genes) > 0){
      gene_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
      if (nrow(gene_entrez) > 0) {
        ego <- enrichGO(gene = gene_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                        pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
        if (!is.null(ego) && nrow(ego) > 0) res_go[[i]] <- ego
      }
    }
  }
  qsave(res_go, go_save_path)
}

if (length(res_go) > 0) write.xlsx(as.data.frame(merge_result(res_go)), file = excel_path, overwrite = TRUE)


## 4.3 Individual Bar Plots ----------------------------------------------------
go_df <- as.data.frame(merge_result(res_go))

if(nrow(go_df) > 0) {
  go_df$Cluster <- as.character(go_df$Cluster)
  valid_clusters <- unique(go_df$Cluster)
  
  for (i in valid_clusters) {
    plot_data <- go_df %>% 
      filter(Cluster == i) %>% 
      arrange(p.adjust) %>% 
      head(top_n_go) %>% 
      mutate(logP = -log10(p.adjust)) %>%
      mutate(Description = str_wrap(Description, width = 60)) # 防止文字过长
    
    # Lock Factor Levels
    plot_data$Description <- factor(plot_data$Description, levels = rev(plot_data$Description))
    
    # Plot (Using your logic: Text on bars, no axis text)
    p <- ggplot(plot_data, aes(x = logP, y = Description)) +
      geom_bar(stat = "identity", width = 0.8, alpha = 0.8, fill = cluster_colors[[i]]) +
      geom_text(aes(x = 0.05, label = Description), size = 5, hjust = 0, color = "black") + 
      theme_minimal() +
      labs(title = paste0("Cluster ", i, " Top ", top_n_go, " GO Terms"), x = "-log10(p.adjust)", y = NULL) +
      theme(
        axis.text.y = element_blank(), # Hide Y axis text
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      scale_x_continuous(expand = c(0, 0, 0.4, 0)) # Extend X axis for text space
    
    safe_name <- gsub("[^a-zA-Z0-9]", "_", i) 
    ggsave(paste0(bar_plot_dir, "Cluster_", safe_name, "_barplot.pdf"), plot = p, width = 8, height = 5)
  }
}


## 4.4 Combined Text Block Plot ------------------------------------------------
plot_list <- list()

# Loop through ALL clusters to handle "No Enrichment" cases
for (i in all_clusters) {
  
  plot_data <- go_df %>% 
    filter(Cluster == i) %>% 
    arrange(p.adjust) %>% 
    head(top_n_go)
  
  if (nrow(plot_data) > 0) {
    plot_data <- plot_data %>% mutate(y_pos = rev(seq_len(n())), x_pos = 0.05)
    
    p_block <- ggplot(plot_data) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = cluster_colors[[i]], alpha = 0.7) +
      geom_text(aes(x = x_pos, y = y_pos, label = Description), size = 5, hjust = 0, color = "black", fontface = "bold") +
      theme_void() +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0.5, top_n_go + 0.5), expand = c(0, 0)) +
      theme(plot.margin = margin(t = 2, r = 0, b = 2, l = 0, unit = "pt"))
    
  } else {
    # Handle clusters with no results
    p_block <- ggplot() +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = top_n_go + 1, fill = "grey95") +
      annotate("rect", xmin = 0, xmax = 0.02, ymin = 0, ymax = top_n_go + 1, fill = cluster_colors[[i]]) +
      annotate("text", x = 0.5, y = (top_n_go + 1)/2, label = paste(i, ": No Enrichment"), size = 5, color = "grey50") +
      theme_void() +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0.5, top_n_go + 0.5), expand = c(0, 0))
  }
  
  plot_list[[i]] <- p_block
}

if(length(plot_list) > 0){
  combined_plot <- wrap_plots(plot_list, ncol = 1) & theme(plot.margin = margin(0,0,0,0))
  ggsave(merged_pdf, combined_plot, width = 6, height = 18, limitsize = FALSE)
}
message("Module 4 Done.")

############# ********************* #####################

# MODULE: Single-cell Sample & Composition Analysis----------------

## 1. Parameters & Environment -------------------------------------------------
group_col <- 'Study'
type_col  <- 'CellType1'
output_dir <- "./results/figure/"

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

manual_colors <- c(
  "MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A", 
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628"
)

## 2. DimPlot Split by Meta ----------------------------------------------------
p1 <- DimPlot(tmp_obj, reduction = "umap", group.by = group_col, 
              raster = FALSE, label = FALSE, pt.size = 0.1) +
  ggtitle(group_col) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2.5))) +
  tidydr::theme_dr(xlength = 0.2, ylength = 0.2,
                   arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")) &
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none",
        axis.title.x = element_text(hjust = 0), axis.title.y = element_text(hjust = 0))

p2 <- DimPlot(tmp_obj, group.by = group_col, split.by = group_col) + 
  NoAxes() + theme(plot.title = element_blank())

pp12 <- (ggplotify::as.ggplot(p1) | ggplotify::as.ggplot(p2)) + plot_layout(widths = c(1, 6))

ggsave(paste0(output_dir, "UMAP_Split_", group_col, ".pdf"), pp12, width = 18, height = 3)

## 3. Nuclei Number Statistics -------------------------------------------------
plot_df_num <- tmp_obj@meta.data %>%
  group_by(Group = .data[[group_col]]) %>%
  summarise(Count = n(), Mean_Genes = round(mean(nFeature_RNA, na.rm = TRUE))) %>%
  mutate(Label = paste0(Count, " (", Mean_Genes, ")")) %>%
  # 先按数量排序（数量少的在 Level 前面，数量多的在后面）
  mutate(Group = reorder(Group, Count)) %>%
  # 强制将 "Unknown" 移动到第一个 Level（在 coord_flip 中对应最下方）
  mutate(Group = forcats::fct_relevel(Group, intersect("Unknown", levels(Group))))

p_stat <- ggplot(plot_df_num, aes(x = Group, y = Count, fill = Group)) +
  geom_col(color = "black", width = 0.8) +
  coord_flip() +
  geom_text(aes(label = Label), hjust = 0, nudge_y = 200, size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = NULL, y = NULL, title = "Number of nuclei (mean genes per nucleus)") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "plain"),
        axis.text = element_text(size = 12, color = "black"))

ggsave(paste0(output_dir, "Bar_Numbers_", group_col, ".pdf"), p_stat, width = 8, height = 6)

## 4. Cell Type Percentage -----------------------------------------------------
plot_df_pct <- tmp_obj@meta.data %>%
  # 统一使用 Group 变量名
  group_by(Group = .data[[group_col]], CellType = .data[[type_col]]) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Group) %>%
  mutate(
    Proportion = Count / sum(Count),
    Label = ifelse(Proportion > 0.03, paste0(round(Proportion * 100, 0), "%"), ""),
    # 继承 tmp_obj 中已经定义好的 factor 顺序，不再按数量 reorder
    Group = factor(Group, levels = levels(tmp_obj@meta.data[[group_col]])),
    # 确保 Unknown 始终在最底部 (Level 1)
    Group = forcats::fct_relevel(Group, intersect("Unknown", levels(Group)))
  )

p_stack <- ggplot(plot_df_pct, aes(x = Group, y = Count, fill = CellType)) +
  geom_col(position = "fill", color = "black", width = 0.8) +
  coord_flip() +
  geom_text(aes(label = Label), position = position_fill(vjust = 0.5), size = 3.5) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = manual_colors) +
  labs(x = NULL, y = "Proportion", fill = "Cell Type",
       title = paste0("Cell Type Composition per ", group_col)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "right")

ggsave(paste0(output_dir, "Bar_Composition_", group_col, ".pdf"), p_stack, width = 8, height = 3)
message("Module Done.")

############# ********************* #####################

# MODULE: PYSCENIC----------------
############# ********************* #####################

# MODULE 5: Microglia Subpopulation Scoring (AddModuleScore) -------------------

## 5.1 Parameters --------------------------------------------------------------
input_file      <- 'data/processed_data/mg_cca_meta_v2.qs'
output_fig_dir  <- 'results/figure/'
combined_save   <- paste0(output_fig_dir, 'combined_MG_ModuleScore_FeaturePlot.pdf')

dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

# 更新后的基因集：加入了三个新的MG3相关通路
mg_marker_list <- list(
  Myelination = c('PLP1', 'NFASC', 'BCAS1', 'UGT8', 'ABCA2', 'MBP', 'CNTN2', 'ANK2', 'FA2H', 'MYRF', 'DLG1', 'JAM3', 'MAG', 'SH3TC2', 'MAL', 'PLLP', 'ASPA', 'NDRG1', 'PTN', 'PMP22', 'TPPP', 'NKX6-2', 'QKI', 'SERINC5', 'MTMR2', 'CD9', 'NTRK2', 'SIRT2', 'DICER1'),
  Axonogenesis = c('NRXN1', 'CNTN1', 'PAK3', 'DCLK1', 'ROBO2', 'KALRN', 'MAP2', 'CHL1', 'SPTBN4', 'CTNNA2', 'ATP8A2', 'DCC', 'EPHA6', 'NRCAM', 'NTRK2', 'PARD3', 'CDH2', 'LRRC4C', 'CHN1', 'MAP1B', 'NRXN3', 'NFIB', 'ROBO1', 'ANK3', 'DAB1', 'NCAM1', 'PTPRO', 'APP'),
  Mononuclear_Cell_Proliferation = c('MEF2C', 'DOCK8', 'ST6GAL1', 'IL6ST', 'CCND3', 'SYK', 'CD86', 'TGFBR2', 'CARD11'),
  Signal_Transduction_By_P53_Class_Mediator = c('RPL11', 'RPS20', 'RPS15', 'RPL23', 'RPL37', 'RPL26', 'RPS7', 'RPL5', 'CD74'),
  Receptor_Mediated_Endocytosis = c('RAB31', 'STAB1', 'GRB2', 'SDCBP', 'IGF2R', 'AP2A2', 'FCHO2', 'ITSN1'),
  Macrophage_Activation = c('MFHAS1', 'GPR137B', 'VSIG4', 'HSPD1', 'CTSC', 'HAVCR2')
)

plot_pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(10)


## 5.2 Process Data ------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(qs)
library(cowplot)
library(patchwork)

if(!file.exists(input_file)) stop("Input file not found.")
tmp_obj <- qread(input_file)

score_features <- c()
for (id in names(mg_marker_list)) {
  score_name <- paste0(id, "_Score")
  tmp_obj <- AddModuleScore(
    object = tmp_obj,
    features = list(mg_marker_list[[id]]),
    name = score_name
  )
  score_features <- c(score_features, paste0(score_name, "1"))
}


## 5.3 Generate Plots ----------------------------------------------------------
plots <- lapply(names(mg_marker_list), function(id) {
  feature_name <- paste0(id, "_Score1")
  
  FeaturePlot(tmp_obj, features = feature_name, order = TRUE, cols = plot_pal, raster = FALSE) +
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


## 5.4 Save Results ------------------------------------------------------------
if(length(plots) > 0) {
  # Combined Plot
  temp_p <- FeaturePlot(tmp_obj, features = score_features[1], cols = plot_pal) + 
    theme(legend.text = element_text(color = "transparent"))
  common_legend <- cowplot::get_legend(temp_p + theme(legend.position = "right"))
  
  # 因为现在有9个图，改为3列的网格布局(3x3)更为美观
  combined_p <- wrap_plots(plots, ncol = 3)
  final_plot <- cowplot::plot_grid(combined_p, common_legend, rel_widths = c(1, 0.08))
  
  # 相应地增加了图的高度和宽度
  ggsave(combined_save, final_plot, height = 10, width = 14, bg = '#ffffff')
}

# Individual Plots
for (id in names(plots)) {
  p_single <- plots[[id]] + theme(legend.position = "right", legend.text = element_text(size = 8))
  ggsave(paste0(output_fig_dir, 'ModuleScore_', id, '.pdf'), p_single, height = 5, width = 6, bg = '#ffffff')
}

message("Module 5 Done.")


# MODULE 6: ROE (Relative Observed to Expected) Analysis -----------------------

## 6.1 Parameters --------------------------------------------------------------
input_file     <- 'data/processed_data/mg_cca_meta_v2.qs'
output_fig_dir <- 'results/figure/'

var1 <- "CellType1"
var2 <- "Region"        

output_file <- paste0(output_fig_dir, "ROE_", var1, "_", var2, ".pdf")
dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)


## 6.2 Process Data ------------------------------------------------------------
library(qs)
library(pheatmap)
library(RColorBrewer)

if(!file.exists(input_file)) stop("Input file not found.")
tmp_obj <- qread(input_file)

tab <- table(tmp_obj[[var1, drop=TRUE]], tmp_obj[[var2, drop=TRUE]])

if("Unknown" %in% colnames(tab)){
  tab <- tab[, colnames(tab) != "Unknown", drop = FALSE]
}

row_sum <- rowSums(tab)
col_sum <- colSums(tab)
total   <- sum(tab)

expected <- outer(row_sum, col_sum) / total
roe <- as.matrix(tab / expected)

rows_to_plot <- rownames(roe)
cols_to_plot <- colnames(roe)


## 6.3 Generate & Save Plots ---------------------------------------------------
plot_pal <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(100))

pdf(output_file, width = 6, height = 6)
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
  # scale = 'row',
  main = paste0("ROE: ", var1, " vs ", var2),
  angle_col = "45"
)
dev.off()

message("Module 6 Done.")

# MODULE 7: Literature Signature Scoring & DotPlot -----------------------------

## 7.1 Parameters --------------------------------------------------------------
input_file      <- 'data/processed_data/mg_cca_meta_v2.qs'
geneset_file    <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/saved_data/Genesets_literature.xlsx")
output_fig_dir  <- 'results/figure/'
output_filename <- paste0(output_fig_dir, 'Literature_Signature_DotPlot.pdf')

dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

# Mapping: Name_in_Seurat = Name_in_Excel
target_signatures_map <- c(
  'DAM.ALL'           = 'up_DAM_Keren_shaul_Thrupp',
  'Human_DAM_Silvin'  = 'Human_Adult_DAMs_Silvin',
  'Human_YAM_Silvin'  = 'Human_Adult_YAMs_Silvin',
  'Human_DIM_Silvin'  = 'Human_Adult_DIMs_Silvin'
)

plot_colors <- c("darkblue", "white", "darkred")


## 7.2 Process Data ------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(qs)
library(openxlsx)

if(!file.exists(input_file)) stop("Input file not found.")
tmp_obj <- qread(input_file)

# Prepare Genesets
genesets_df <- read.xlsx(geneset_file)
genesets_df$genesets <- paste0(genesets_df$Population, "_", genesets_df$Study)

# Extract Markers
mg_marker_list <- lapply(target_signatures_map, function(sig_id) {
  return(genesets_df$Gene[genesets_df$genesets == sig_id])
})

# Add Module Scores
score_features <- c()
for (name in names(mg_marker_list)) {
  tmp_obj <- AddModuleScore(
    object = tmp_obj,
    features = list(mg_marker_list[[name]]),
    name = name
  )
  score_features <- c(score_features, paste0(name, "1"))
}


## 7.3 Generate Plots ----------------------------------------------------------
p <- DotPlot(tmp_obj, 
             features = score_features, 
             dot.scale = 10, 
             group.by = 'CellType1') +
  scale_colour_gradient2(low = plot_colors[1], mid = plot_colors[2], high = plot_colors[3]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = NULL, y = NULL)


## 7.4 Save Results ------------------------------------------------------------
ggsave(output_filename, p, height = 5, width = 4.5, bg = '#ffffff')

message("Module 7 Done.")

# MODULE 8: Microglia Myelination Gene Heatmap ---------------------------------

## 8.1 Parameters --------------------------------------------------------------
input_file      <- 'data/processed_data/mg_cca_meta_v2.qs'
output_fig_dir  <- 'results/figure/'
output_filename <- paste0(output_fig_dir, 'MG_myelination_gene_heatmap.pdf')

dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

myelination_genes <- c(
  "ST18", "SLC44A1", "PHLPP1", "PTPRD", "NCAM2", "PLP1", "NCAM1", "KCNMB4",
  "GNAO1", "UGT8", "MOBP", "MYRF", "BCAS1", "MBP", "MOG", "RNF220",
  "TRIM2", "CLDN11", "FA2H", "DLG1", "SH3TC2", "PTN"
)


## 8.2 Process Data ------------------------------------------------------------
library(Seurat)
library(pheatmap)
library(qs)
library(dplyr)

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
# Select genes present in data and reorder to match input list
valid_genes <- intersect(myelination_genes, rownames(avg_mat))
plot_data   <- avg_mat[valid_genes, , drop = FALSE]
plot_data   <- plot_data[match(valid_genes, rownames(plot_data)), ]


## 8.3 Generate & Save Plot ----------------------------------------------------
pdf(output_filename, width = 3, height = 3)
pheatmap(
  plot_data,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  show_rownames = TRUE,
  legend = FALSE
)
dev.off()

message("Module 8 Done.")

# MODULE 9: Microglia Other Function Gene Heatmap ------------------------------

## 9.1 Parameters --------------------------------------------------------------
input_file      <- 'data/processed_data/mg_cca_meta_v2.qs'
output_fig_dir  <- 'results/figure/'
output_filename <- paste0(output_fig_dir, 'MG_other_function_gene_heatmap.pdf')

dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

target_genes <- c(
  "CTNNA3", "CNTNAP4", "SPOCK1", "CADM2", "TJP1", "CNTN2", "CTNND2",
  "RNF220", "C10orf90", "FBXL7", "TRIM2", "PPM1H", "RNF144A", 
  "UNC5C", "TTLL7", "MAGI2", "NTM", "NFASC", "ENPP2", 
  "DNAJC6", "DLC1", "LPGAT1"
)


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
pdf(output_filename, width = 3, height = 3)

pheatmap(
  plot_data,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  show_rownames = TRUE,
  legend = FALSE
)

dev.off()

message("Module 9 Done.")


# MODULE 10: GSEA Analysis (Single Cluster) ------------------------------------

## 10.1 Parameters -------------------------------------------------------------
input_file      <- 'data/processed_data/mg_cca_meta_v2.qs'
dir_data        <- "scripts/GSEA/data/"
dir_res         <- "scripts/GSEA/results/"
custom_set_file <- paste0(dir_data, 'my_gene_set.xlsx')
go_gmt_file     <- paste0(dir_data, "c5.go.v2023.2.Hs.symbols.gmt")

target_group_by <- 'CellType1'
target_cluster  <- 'MG0'  # Modify this to change the cluster

output_qs       <- paste0(dir_data, 'GSEA_kk_', target_cluster, '.qs')
output_txt      <- paste0(dir_res, target_cluster, "_combined_genesets_GSEA_result.txt")
output_pdf_p    <- paste0(dir_res, target_cluster, "_combined_GSEA_select_plot.pdf")
output_pdf_nop  <- paste0(dir_res, target_cluster, "_combined_GSEA_select_plot_no_pvalue.pdf")

dir.create(dir_res, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_data, recursive = TRUE, showWarnings = FALSE)

target_go_terms <- c(
  'GOBP_POSITIVE_REGULATION_OF_MYELINATION', 'GOBP_REGULATION_OF_MYELINATION',
  'GOBP_OLIGODENDROCYTE_DIFFERENTIATION', 'GOBP_CELL_JUNCTION_ORGANIZATION',
  'GOBP_ADAPTIVE_IMMUNE_RESPONSE', 'GOBP_FATTY_ACID_TRANSPORT',
  'GOBP_MYELIN_ASSEMBLY', 'GOBP_OLIGODENDROCYTE_DEVELOPMENT',
  'GOBP_LIPID_CATABOLIC_PROCESS'
)

selected_plot_ids <- c(
  'oligodendrocyte development',
  'oligodendrocyte differentiation',
  'myelin assembly',
  'cell junction organization'
  # 'positive regulation of myelination'
)
# plot_colors <- c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D23")
plot_colors <- c("#D51F26", "#272E6A", "#208A42", "#89288F")


## 10.2 Process Data -----------------------------------------------------------
library(Seurat)
library(clusterProfiler)
library(openxlsx)
library(dplyr)
library(qs)
library(GseaVis)

if(!file.exists(input_file)) stop("Input file not found.")
tmp_obj <- qread(input_file)

# Prepare Gene Sets
gene_sets_raw <- read.xlsx(custom_set_file, colNames = FALSE, na.strings = "")
temp_gmt_path <- paste0(dir_data, 'temp_my_genesets.gmt')
write.table(gene_sets_raw, temp_gmt_path, sep = "\t", row.names = F, col.names = F, quote = F)

custom_gmt <- read.gmt(temp_gmt_path)
custom_gmt$gene <- toupper(custom_gmt$gene)

gmt_go <- read.gmt(go_gmt_file)
filtered_go <- gmt_go %>%
  filter(term %in% target_go_terms) %>%
  mutate(term = gsub("GOBP_", "", term)) %>%
  mutate(term = tolower(gsub("_", " ", term)))

combined_gmt <- rbind(custom_gmt, filtered_go)

# Prepare Ranked Gene List
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

# Run GSEA
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


## 10.3 Generate Plots ---------------------------------------------------------
p_gsea_pval <- gseaNb(
  object    = kk, 
  geneSetID = selected_plot_ids,
  subPlot   = 3, 
  curveCol  = plot_colors,
  addPval   = TRUE,
  pvalSize  = 5,
  base_size = 15
)

p_gsea_clean <- gseaNb(
  object    = kk, 
  geneSetID = selected_plot_ids,
  subPlot   = 3, 
  curveCol  = plot_colors,
  addPval   = FALSE,
  pvalSize  = 5,
  base_size = 15
)


## 10.4 Save Results -----------------------------------------------------------
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

message(paste0("Module 10 Done for cluster: ", target_cluster))

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
