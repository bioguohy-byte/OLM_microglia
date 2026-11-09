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


# Figure 4

# F4A  Cell-cell-contact circle network, MG0-MG5 sources and major brain-cell targets
# Kind: computational
mg_register_panel("F4A", list(
  list(source="source/project/code_upload/CellChat.r", ranges="1-93", requested_ranges="1-93", role="prepared-upload setup, CellChat calculation and complete circle graphics device", code=quote({
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
  }))
))

# F4B  Three-block selected ligand-receptor bubbles: OLM->major cells, MG->OLIG, OLIG->MG
# Kind: computational
mg_register_panel("F4B", list(
  list(source="source/project/code_upload/CellChat.r", ranges="1-8;38-83;113-146", requested_ranges="1-8;38-83;113-146", role="historical bubble code; cellchat object comes from the preceding calculation or its saved checkpoint", code=quote({
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
  }))
))

# F4C  Integrated three-slice MERFISH microglial UMAP ST-OLM/ST_MG1-ST_MG5
# Kind: computational
mg_register_panel("F4C", list(
  list(source="source/project/Reference/combine/scripts/Spatial/scripts/spatial_module_verify_single.r", ranges="1-62", requested_ranges="1-62", role="existing checkpoint verification/redraw implementation", code=quote({
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
  })),
  list(source="source/project/code_upload/Spatial.r", ranges="1-16;78-134", requested_ranges="1-16;78-134", role="Historical with3-to-with4 processing block", code=quote({
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
  }))
))

# F4D  myelin1 module-score violins across integrated MERFISH ST_MG0-ST_MG5
# Kind: computational
mg_register_panel("F4D", list(
  list(source="source/project/Reference/combine/scripts/Spatial/scripts/spatial_module_verify_single.r", ranges="1-62", requested_ranges="1-62", role="existing checkpoint verification/redraw implementation", code=quote({
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
  })),
  list(source="source/project/code_upload/Spatial.r", ranges="1-16;78-134", requested_ranges="1-16;78-134", role="Historical with3-to-with4 processing block", code=quote({
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
  }))
))

# F4E  s2r2 MERFISH ROI and magnified inset; OLG blue, OPC yellow, ST OLM red, others gray
# Kind: computational
mg_register_panel("F4E", list(
  list(source="source/project/Reference/combine/scripts/Spatial/scripts/plot.r", ranges="1-128", requested_ranges="1-128", role="complete overview setup and plot", code=quote({
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
  })),
  list(source="source/project/Reference/combine/scripts/Spatial/scripts/plot_zoom.r", ranges="1-128", requested_ranges="1-128", role="complete zoom implementation; manual montage remains separate", code=quote({
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
  }))
))

# F4F  Cross-type centered L(r), OLG versus ST OLM and ST_MG1-ST_MG4
# Kind: computational
mg_register_panel("F4F", list(
  list(source="source/project/code_upload/Spatial.r", ranges="1-16;106-113;136-205", requested_ranges="1-16;106-113;136-205", role="", code=quote({
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
  }))
))
