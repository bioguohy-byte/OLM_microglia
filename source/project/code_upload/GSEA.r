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