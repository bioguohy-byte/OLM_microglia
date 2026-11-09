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
############# ********************* #####################
# MODULE 2: DotPlot (Visualization) --------------------------------------------

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