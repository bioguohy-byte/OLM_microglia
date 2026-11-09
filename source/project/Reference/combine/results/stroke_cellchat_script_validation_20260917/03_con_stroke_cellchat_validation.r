# Validation-only copy of the computational and circle-plot sections of:
#   combine/scripts/CellChat/scripts/02.con_stroke_cellchat.r
# Original objects and figures are not overwritten.

library(CellChat)
library(Seurat)
library(qs)
library(ggplot2)
library(patchwork)
library(igraph)
library(ComplexHeatmap)

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
validation_root <- file.path(project_root, "combine/results/stroke_cellchat_script_validation_20260917")
input_file <- file.path(validation_root, "stroke_cellchat_reproduced.qs")
cellchat_output <- file.path(validation_root, "con_stroke_cellchat_reproduced.qs")
figure_output <- file.path(validation_root, "stroke_cellchat_circle_multi_reproduced.pdf")

# Original Cell-Cell Contact CellChat calculation.
seu <- qread(input_file)
# Compatibility-only translation of the original `slot = "data"` call for
# SeuratObject >= 5; it selects the same normalized RNA expression matrix.
seu.data.input <- GetAssayData(seu, assay = "RNA", layer = "data")
seu.identity <- subset(seu@meta.data, select = "cellchat_type")
cellchat <- createCellChat(object = seu.data.input, meta = seu.identity, group.by = "cellchat_type")
cellchat <- addMeta(cellchat, meta = seu.identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = c("Cell-Cell Contact"))
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
qsave(cellchat, cellchat_output)

# Exact circle-plot settings used for the multi-MG source panel in the original.
clusters_source_multi <- c("MG0", "MG1", "MG2", "MG3", "MG4", "MG5")
cell_colors <- c(
  "MG0" = "#4DAF4A", "MG1" = "#A65628", "MG2" = "#377EB8",
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG5" = "#8B0000",
  "Neuron" = "#E31A1C", "Astrocyte" = "#FF7F00", "OLIG" = "#1F78B4",
  "OPC" = "#A6CEE3", "Ependymal" = "#6A3D9A", "Vascular" = "#FB9A99",
  "T cell" = "#FDBF6F"
)
net_weight_no_self <- cellchat@net$weight
diag(net_weight_no_self) <- 0

pdf(figure_output, width = 5, height = 5)
netVisual_circle(
  net_weight_no_self,
  arrow.size = 0.5, arrow.width = 0.8, vertex.weight = 13,
  weight.scale = TRUE, label.edge = FALSE, color.use = cell_colors,
  title.name = "Interaction weights/strength", sources.use = clusters_source_multi
)
invisible(dev.off())

cat("Reproduced CellChat object:\n", cellchat_output, "\n", sep = "")
cat("Reproduced circle plot:\n", figure_output, "\n", sep = "")
cat("Node order:\n", paste(levels(cellchat@idents), collapse = " | "), "\n", sep = "")
