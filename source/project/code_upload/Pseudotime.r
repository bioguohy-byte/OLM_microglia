# Pseudotime -------------------
library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)
library(stringr)
library(qs)
library(ClusterGVis)
library(org.Hs.eg.db)
library(colorspace)

R.utils::setOption("clusterProfiler.download.method", 'auto')

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