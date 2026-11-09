setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/"))
rm(list=ls())

# 0. Libraries -----------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Seurat)
library(Matrix)
library(cowplot)
library(reshape2)
library(SingleCellExperiment)
library(umap)
library(ggthemes)
library(dplyr)
library(patchwork)
library(ggsci)
library(scales)
library(tidyverse)
library(magick)
library(magrittr)
library(devtools)
library(ComplexHeatmap)
library(ggplotify)
library(circlize)
library(pheatmap)
library(harmony)
library(qs)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)

# load data-------------------------------
tmp_obj <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta.qs"))
DimPlot(tmp_obj,group.by = 'integrated_snn_res.0.1')

# FindAllMarkers-----------------
# 1. Parameters ----------------------------------------------------------------
cluster_col <- 'integrated_snn_res.0.3'      
min_pct     <- 0.25                          
logfc_thr   <- 0.25                          

# Auto-generate filenames
res_suffix  <- gsub("integrated_snn_", "", cluster_col)
qs_path     <- paste0('data/saved_data/', res_suffix, '_markers.qs')
xlsx_path   <- paste0('results/table/', res_suffix, '_markers.xlsx')

# 2. Run Analysis --------------------------------------------------------------
Idents(tmp_obj) <- cluster_col

sub_clusters_markers <- FindAllMarkers(
  tmp_obj,
  only.pos = TRUE,
  min.pct = min_pct,
  logfc.threshold = logfc_thr
)

# 3. Save Results --------------------------------------------------------------
qsave(sub_clusters_markers, qs_path)
write.xlsx(sub_clusters_markers, file = xlsx_path, overwrite = TRUE)

message(paste("Markers saved to:", qs_path))


# dotplot---------------
# — 1. Parameters ----------------------------------------------------------------
top_n_num   <- 10                            # 取 Top 几
cluster_col <- 'integrated_snn_res.0.1'      # 分组列名

# Auto-generate filenames
res_suffix <- gsub("integrated_snn_", "", cluster_col)
excel_path <- paste0('results/table/', res_suffix, '_top', top_n_num, '_markers.xlsx')
png_path   <- paste0("results/figure/", res_suffix, "_top", top_n_num, "_dotplot.png")
pdf_path   <- paste0("results/figure/", res_suffix, "_top", top_n_num, "_dotplot.pdf")

# — 2. Process Data --------------------------------------------------------------
sub_clusters_markers <- qread('data/saved_data/sub_clusters_markers.qs')

top_markers <- sub_clusters_markers %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = top_n_num) %>%
  ungroup()

write.xlsx(top_markers, file = excel_path, overwrite = TRUE)
genes_list <- unique(top_markers$gene)

# — 3. Visualization -------------------------------------------------------------
DefaultAssay(tmp_obj) <- 'RNA'
p_base <- DotPlot(tmp_obj, features = genes_list, group.by = cluster_col)

p_final <- p_base$data %>%
  dplyr::rename(pct_exp = pct.exp, avg_exp = avg.exp.scaled, Gene = features.plot, Cluster = id) %>%
  ggplot(aes(x = Gene, y = Cluster, size = pct_exp, color = avg_exp)) +
  geom_point() +
  scale_color_gradient2(low = "#000000", mid = "#007965", high = "#ffcc29", midpoint = 0) +
  theme_classic() +
  labs(title = paste0("Top ", top_n_num, " Genes in ", res_suffix), x = NULL, y = NULL) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

# — 4. Save ----------------------------------------------------------------------
ggsave(png_path, plot = p_final, width = 12, height = 3, limitsize = FALSE)
ggsave(pdf_path, plot = p_final, width = 12, height = 3, limitsize = FALSE)


# heatmap--------------------
# 1. Parameters ----------------------------------------------------------------
top_n       <- 20                            
cluster_col <- 'integrated_snn_res.0.1'      

# Auto-generate filenames
res_suffix  <- gsub("integrated_snn_", "", cluster_col)
png_path    <- paste0("results/figure/", res_suffix, "_top", top_n, "_heatmap.png")
pdf_path    <- paste0("results/figure/", res_suffix, "_top", top_n, "_heatmap.pdf")

# 2. Process Data --------------------------------------------------------------
tmp_obj <- qread('data/processed_data/mg_cca_meta.qs')
sub_clusters_markers <- qread('data/saved_data/sub_clusters_markers.qs')

# Filter Top N
top_markers <- sub_clusters_markers %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup()

# Calculate Average Expression
DefaultAssay(tmp_obj) <- 'RNA'
avg_exp_obj <- AverageExpression(tmp_obj, group.by = cluster_col, assays = "RNA")
avg_mat <- avg_exp_obj$RNA

# Match & Sort Matrix by Marker Order
heatmap_data <- avg_mat[top_markers$gene, ]

# 3. Visualization -------------------------------------------------------------
p_heatmap <- pheatmap(heatmap_data, 
                      cluster_cols = FALSE, 
                      cluster_rows = FALSE,
                      scale = "row",        
                      show_rownames = FALSE, 
                      legend = FALSE,         
                      silent = TRUE,
                      angle_col = 0)

p_final <- as.ggplot(p_heatmap)

# 4. Save ----------------------------------------------------------------------
ggsave(png_path, plot = p_final, width = 5, height = 18, limitsize = FALSE)
ggsave(pdf_path, plot = p_final, width = 5, height = 18, limitsize = FALSE)




# GO--------------------------
# 1. Parameters ----------------------------------------------------------------
cluster_col <- 'integrated_snn_res.0.1'      
top_n_go    <- 5                             # Number of GO terms to plot
p_val_cut   <- 0.05
q_val_cut   <- 0.05

# Define Colors based on biological identity (0-DAM, 1-HOM, 2-MAM, 3-RIB, 4-SYM)
cluster_colors <- c(
  "0" = "#E41A1C", # Red (Disease/Alert)
  "1" = "#377EB8", # Blue (Homeostatic/Calm)
  "2" = "#4DAF4A", # Green (Myelin/Phagocytosis)
  "3" = "#984EA3", # Purple (Ribosome/High Activity)
  "4" = "#FF7F00"  # Orange (Synapse/Connection)
)

# Auto-generate filenames
res_suffix  <- gsub("integrated_snn_", "", cluster_col)
go_save_path <- paste0("data/saved_data/", res_suffix, "_res_go_list.qs")
excel_path   <- paste0('results/table/', res_suffix, '_GO_enrichment.xlsx')
merged_pdf   <- paste0("results/figure/", res_suffix, "_GO_Combined_TextBlocks.pdf")
merged_png   <- paste0("results/figure/", res_suffix, "_GO_Combined_TextBlocks.png")
dir.create("results/figure/cluster_go/", recursive = TRUE, showWarnings = FALSE)

# 2. GO Calculation ------------------------------------------------------------
# Load Markers
sub_clusters_markers <- qread('data/saved_data/sub_clusters_markers.qs')

# Filter Significant Markers
filt_markers <- sub_clusters_markers %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>% # Only Upregulated
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  ungroup()

# Run Loop (Check if result exists to avoid re-running slow steps)
if (file.exists(go_save_path)) {
  message("Loading existing GO results...")
  res_go <- qread(go_save_path)
} else {
  message("Running GO Enrichment...")
  res_go <- list()
  clusters <- levels(as.factor(filt_markers$cluster))
  
  for (i in clusters) {
    message(paste("Processing Cluster:", i))
    genes <- filt_markers$gene[filt_markers$cluster == i]
    
    # ID Conversion (Symbol -> Entrez)
    gene_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    
    if (nrow(gene_entrez) > 0) {
      ego <- enrichGO(
        gene          = gene_entrez$ENTREZID,
        OrgDb         = org.Hs.eg.db,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = p_val_cut,
        qvalueCutoff  = q_val_cut,
        readable      = TRUE
      )
      if (!is.null(ego)) res_go[[i]] <- ego
    }
  }
  qsave(res_go, go_save_path)
}

# Save Merged Table
if (length(res_go) > 0) {
  merge_res_go <- merge_result(res_go)
  write.xlsx(as.data.frame(merge_res_go), file = excel_path, overwrite = TRUE)
}

# 3. Individual Bar Plots ------------------------------------------------------
# Convert to Data Frame for plotting
go_df <- as.data.frame(merge_result(res_go))

for (i in names(cluster_colors)) {
  # Filter Data
  plot_data <- go_df %>%
    filter(Cluster == i) %>%
    arrange(p.adjust) %>%
    head(top_n_go) %>%
    mutate(logP = -log10(p.adjust))
  
  if (nrow(plot_data) == 0) next
  
  # Lock Factor Levels for Order
  plot_data$Description <- factor(plot_data$Description, levels = rev(plot_data$Description))
  
  # Plot
  p <- ggplot(plot_data, aes(x = logP, y = Description)) +
    geom_bar(stat = "identity", width = 0.8, alpha = 0.8, fill = cluster_colors[i]) +
    geom_text(aes(x = 0.05, label = Description), size = 7, hjust = 0, color = "black") +
    theme_minimal() +
    labs(title = paste0("Cluster ", i, " Top ", top_n_go, " GO Terms"), x = "-log10(p.adjust)", y = NULL) +
    theme(
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    scale_x_continuous(expand = c(0, 0, 0.2, 0))
  
  # Save Individual Plot
  ggsave(paste0('results/figure/cluster_go/', res_suffix, '_Cluster_', i, '.pdf'), 
         plot = p, width = 8, height = 5)
}

# 4. Combined Text Block Plot (Patchwork) --------------------------------------
plot_list <- list()
valid_clusters <- intersect(names(cluster_colors), unique(go_df$Cluster))

for (i in valid_clusters) {
  
  plot_data <- go_df %>%
    filter(Cluster == i) %>%
    arrange(p.adjust) %>%
    head(top_n_go)
  
  # Coordinate setup for text
  plot_data <- plot_data %>%
    mutate(y_pos = rev(seq_len(nrow(plot_data))), x_pos = 0.05)
  
  p_block <- ggplot(plot_data) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
             fill = cluster_colors[i], alpha = 0.7) +
    geom_text(aes(x = x_pos, y = y_pos, label = Description), 
              size = 5, hjust = 0, color = "black", fontface = "bold") +
    theme_void() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, top_n_go + 0.5), expand = c(0, 0)) +
    theme(plot.margin = margin(t = 2, r = 0, b = 2, l = 0, unit = "pt"))
  
  plot_list[[i]] <- p_block
}

# Combine
combined_plot <- wrap_plots(plot_list, ncol = 1) & theme(plot.margin = margin(0,0,0,0))

# Dynamic Height: 1.5 inch per cluster
ggsave(merged_pdf, combined_plot, width = 4, height = 18, limitsize = FALSE)
ggsave(merged_png, combined_plot, width = 4, height = 18, dpi = 300, limitsize = FALSE)

message("GO Analysis and Plotting Completed.")


manual_colors <- c(
  "HOM" = "#377EB8", "DAM" = "#E41A1C", "MAM" = "#4DAF4A", 
  "RAM" = "#984EA3", "IAM" = "#E6AB02", "SAM" = "#A65628"
)

DimPlot(tmp_obj,group.by = 'CellType1',cols = manual_colors,label = T)




# 小胶质细胞亚群定义 (Microglia Sub-cluster Definition) ----------
# 将聚类结果转换为字符型，方便后续重命名
tmp_obj$sub_clusters <- as.character(tmp_obj$integrated_snn_res.0.2)

tmp_obj@meta.data %<>% 
  mutate(sub_clusters = case_when(
    integrated_snn_res.0.2 == "0"  ~ "HOM-1", # 稳态小胶质细胞 1 (Homeostatic 1)
    integrated_snn_res.0.2 == "3"  ~ "HOM-2", # 稳态小胶质细胞 2 (Homeostatic 2)
    integrated_snn_res.0.2 == "6"  ~ "DAM-1", # 疾病相关小胶质细胞 1 (Disease-Associated 1)
    integrated_snn_res.0.2 == "9"  ~ "DAM-2", # 疾病相关小胶质细胞 2 (Disease-Associated 2)
    integrated_snn_res.0.2 == "1"  ~ "ARM-1", # 激活反应性小胶质细胞 1 (Activated Response 1)
    integrated_snn_res.0.2 == "7"  ~ "ARM-2", # 激活反应性小胶质细胞 2 (Activated Response 2)
    integrated_snn_res.0.2 == "4"  ~ "PIM",   # 促炎亚群 (Pro-inflammatory)
    integrated_snn_res.0.2 == "2"  ~ "MAM",   # 髓鞘相关亚群 (Myelin-Associated)
    integrated_snn_res.0.2 == "5"  ~ "SYM",   # 突触调节亚群 (Synapse-Regulating)
    integrated_snn_res.0.2 == "8"  ~ "SRM",   # 应激反应亚群 (Stress-Responsive)
    integrated_snn_res.0.2 == "10" ~ "BAM",   # 边界相关巨噬细胞 (Border-Associated Macrophages)
    integrated_snn_res.0.2 == "11" ~ "CRM",   # 趋化因子富集 (Chemokine-Rich)
    integrated_snn_res.0.2 == "12" ~ "NAM"    # 神经元相关亚群 (Neuron-Associated)
  ))

# 统计每个亚群的细胞数量，检查是否有遗漏 (NA)
table(tmp_obj$sub_clusters)

qsave(tmp_obj,'data/processed_data/mg_sub_clusters.qs')

# --- 颜色设置 (Color Palette) ---
# 为功能相近的亚群设置同色系但不同深浅的颜色
mg_colors <- c(
  # 稳态蓝系：展现基础生理状态
  "HOM-1" = "#4E79A7", # 中蓝色
  "HOM-2" = "#A0CBE8", # 浅青蓝色
  
  # 疾病红/橙系：突出病理相关性
  "DAM-1" = "#FB8072", # 珊瑚红
  "DAM-2" = "#FCBBA1", # 浅肉粉色
  
  # 激活紫系：表现免疫反应过渡
  "ARM-1" = "#807DBA", # 深紫色
  "ARM-2" = "#BEBADA", # 浅紫色
  
  # 其他功能特化亚群
  "PIM"   = "#FDB462", # 橙色 (促炎)
  "MAM"   = "#B3DE69", # 草绿色 (髓鞘)
  "SYM"   = "#FCCDE5", # 浅粉色 (突触)
  "SRM"   = "#D9D9D9", # 灰色 (应激/热休克)
  "BAM"   = "#BC80BD", # 丁香紫 (巨噬细胞)
  "CRM"   = "#E41A1C", # 正红色 (趋化因子富集)
  "NAM"   = "#FFED6F"  # 柠檬黄 (神经元交互)
)


# --- 绘图 (UMAP Visualization) ---
library(ggplot2)

p <- DimPlot(tmp_obj, 
             group.by = 'sub_clusters', 
             cols = mg_colors,      # 应用自定义色系
             reduction = "umap",    # 使用 UMAP 降维结果
             label = TRUE,          # 在图上显示亚群名称
             label.size = 6,        # 标签字体大小
             repel = TRUE) +        # 自动排布标签防止重叠
  labs(title = "Microglia Functional Sub-clusters", # 图片标题
       x = "UMAP_1", y = "UMAP_2") +                 # 坐标轴名称
  theme_bw() +                                       # 使用白底黑框主题
  theme(
    panel.grid = element_blank(),                    # 移除背景网格线，使图面整洁
    axis.text = element_text(color = "black"),       # 坐标轴文字设为黑色
    panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"), # 加深外框
    plot.title = element_text(hjust = 0.5, face = "bold"), # 标题居中并加粗
    legend.text = element_text(size = 10)            # 图例文本大小
  )

# 在屏幕上打印图片
print(p)

# 如果需要保存图片用于论文，取消下面一行的注释
ggsave("results/figure/Microglia_Detailed_Subclusters.pdf", p, width = 8, height = 6)



library(ggplot2)

p <- DimPlot(tmp_obj, 
             group.by = 'sub_clusters', 
             cols = mg_colors,      # 应用自定义色系
             reduction = "umap",    # 使用 UMAP 降维结果
             label = TRUE,          # 在图上显示亚群名称
             label.size = 0,        # 标签字体大小
             repel = TRUE) +        # 自动排布标签防止重叠
  labs(title = "Microglia Functional Sub-clusters", # 图片标题
       x = "UMAP_1", y = "UMAP_2") +                 # 坐标轴名称
  theme_bw() +                                       # 使用白底黑框主题
  theme(
    panel.grid = element_blank(),                    # 移除背景网格线，使图面整洁
    axis.text = element_text(color = "black"),       # 坐标轴文字设为黑色
    panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"), # 加深外框
    plot.title = element_text(hjust = 0.5, face = "bold"), # 标题居中并加粗
    legend.text = element_text(size = 10)            # 图例文本大小
  )
p
ggsave("results/figure/Microglia_Detailed_Subclusters_no_label.pdf", p, width = 8, height = 6)


# ----------------
# 小胶质myelination基因集热图heatmap--------------
gene_list <- c("ST18", "SLC44A1", "PHLPP1", "PTPRD", "NCAM2", "PLP1", "NCAM1", "KCNMB4",
               "GNAO1", "UGT8", "MOBP", "MYRF", "BCAS1", "MBP", "MOG", "OPALIN",
               "QKI", "CLDN11", "FA2H", "DLG1", "SH3TC2", "PTN")
Exp_av <- AverageExpression(tmp_obj,
                            group.by = "sub_clusters",
                            assays = "RNA")

Exp_av <- Exp_av[[1]]
Exp_av <- as.data.frame(Exp_av)
Exp_av$gene <- rownames(Exp_av)
result <- Exp_av %>%
  filter(gene %in% gene_list)

result <- result[match(gene_list, result$gene),]
#现在的result中有一列是gene列，所以不能直接放在pheatmap中进行绘图，#因为里面不全是数字了，所以要删除那一列
result <- subset(result, select = -gene)

pdf("results/figure/MG_myelination_gene_heatmap.pdf", width = 4, height = 4)
pheatmap(result, cluster_cols = FALSE, cluster_rows = FALSE,scale = "row",show_rownames = T,legend = F)
dev.off()

# cluster_cols = FALSE: 这个参数决定是否对列进行聚类。如果设置为FALSE，那么列就不会被聚类。
# cluster_rows = FALSE: 这个参数决定是否对行进行聚类。如果设置为FALSE，那么行就不会被聚类。
# scale = "row": 这个参数决定是否对数据进行缩放，以及如何进行缩放。如果设置为"row"，那么每一行的数据都会被缩放，使得每一行的平均值为0，标准差为1。这样可以帮助你更好地比较不同基因的表达模式。

# 小胶质其他功能基因集热图heatmap--------------
gene_list <- c("CTNNA3", "CNTNAP4", "SPOCK1", "CADM2", "TJP1", "CNTN2", "CTNND2",
               "RNF220", "C10orf90", "FBXL7", "TRIM2", "PPM1H", "RNF144A", 
               "UNC5C","TTLL7", "MAGI2", "NTM", "NFASC", "ENPP2", 
               "DNAJC6", "DLC1", "LPGAT1")

Exp_av <- AverageExpression(tmp_obj,
                            group.by = "sub_clusters",
                            assays = "RNA")

Exp_av <- Exp_av[[1]]
Exp_av <- as.data.frame(Exp_av)
Exp_av$gene <- rownames(Exp_av)
result <- Exp_av %>%
  filter(gene %in% gene_list)

result <- result[match(gene_list, result$gene),]
#现在的result中有一列是gene列，所以不能直接放在pheatmap中进行绘图，#因为里面不全是数字了，所以要删除那一列
result <- subset(result, select = -gene)

pdf("results/figure/MG_other_function_gene_heatmap.pdf", width = 4, height = 4)
pheatmap(result, cluster_cols = FALSE, cluster_rows = FALSE,scale = "row",show_rownames = T,legend = F)
dev.off()

# cluster_cols = FALSE: 这个参数决定是否对列进行聚类。如果设置为FALSE，那么列就不会被聚类。
# cluster_rows = FALSE: 这个参数决定是否对行进行聚类。如果设置为FALSE，那么行就不会被聚类。
# scale = "row": 这个参数决定是否对数据进行缩放，以及如何进行缩放。如果设置为"row"，那么每一行的数据都会被缩放，使得每一行的平均值为0，标准差为1。这样可以帮助你更好地比较不同基因的表达模式。