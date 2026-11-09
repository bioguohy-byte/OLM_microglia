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
message("Module 2 Done.")
