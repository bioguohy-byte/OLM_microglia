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


# Figure 2

# F2A  Condition and healthy-age MG0-MG5 100% composition
# Kind: computational
mg_register_panel("F2A", list(
  list(source="source/project/Reference/combine/scripts/age_MG0_MG5_composition_disease_and_Health_age_from_original_module4.R", ranges="1-101", requested_ranges="1-101", role="", code=quote({
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
  }))
))

# F2B  23 myelin/OLM-associated genes across MG0-MG5
# Kind: computational
mg_register_panel("F2B", list(
  list(source="source/project/Reference/combine/scripts/20260520_dotplot.r", ranges="24-65", requested_ranges="24-65", role="", code=quote({
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
  }))
))

# F2C  TRIM2 RNF220 MBP MOG FeaturePlots
# Kind: computational
mg_register_panel("F2C", list(
  list(source="source/project/Reference/combine/scripts/20260520_featureplot.r", ranges="18-30", requested_ranges="18-30", role="", code=quote({
input_file      <- 'data/processed_data/mg_cca_meta_v2.qs'
output_fig_dir  <- 'results/figure/'
save_path   <- paste0(output_fig_dir, '20260520_featurePlot.pdf')

dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)

gene_list <- c('TRIM2','RNF220','MBP','MOG')

plot_pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(10)
tmp_obj <- qread(input_file)

plot <- FeaturePlot(tmp_obj,features = gene_list,cols = plot_pal, raster = FALSE)
ggsave(save_path, plot, width = 7.5, height = 6, limitsize = FALSE)
  }))
))

# F2D  Nine-term OLM GSEA
# Kind: computational
mg_register_panel("F2D", list(
  list(source="source/project/code_upload/GSEA.r", ranges="1-146", requested_ranges="1-146", role="", code=quote({
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
  }))
))
