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

message("Module 9 Done.")



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
