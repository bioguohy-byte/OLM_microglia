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

# Fig.02 ---------------
pdf(paste0(out_dir, "stroke_cellchat_circle_olig.pdf"), width = 5, height = 5)
netVisual_circle(net_weight_no_self,
                 arrow.size = 0.5, arrow.width = 0.8, vertex.weight = 13,
                 weight.scale = TRUE, label.edge = FALSE,
                 color.use = cell_colors, title.name = "Interaction weights/strength",
                 sources.use = cluster_olig)
invisible(dev.off())

# Fig.03 ---------------
pdf(paste0(out_dir, "stroke_cellchat_circle_subset.pdf"), width = 5, height = 5)
netVisual_circle(mat_subset,
                 arrow.size = 1, arrow.width = 0.8, vertex.weight = 13,
                 weight.scale = TRUE, label.edge = FALSE,
                 color.use = new_col_subset, title.name = "Interaction weights/strength",
                 edge.width.max = 15)
invisible(dev.off())

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