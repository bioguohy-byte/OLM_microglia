rm(list=ls())

# 0. library --------------------------
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
library(enrichplot)
library(GseaVis)


# 跨物种小胶质细胞整合---------------
# 提取CPZ小鼠小胶质细胞（白质）
mouse_wm <- qread(Sys.getenv("MG_EXTERNAL_019"))
mouse_wm_mg <- subset(mouse_wm,subset = integrated_snn_res.0.2 %in% c('1','5','17'))
mouse_wm_mg$Species <- as.factor("Mouse")

# 提取衰老猕猴小胶质细胞（海马）
monkey_hip <- qread(Sys.getenv("MG_EXTERNAL_020"))
monkey_hip_mg <- subset(monkey_hip,subset = integrated_snn_res.0.1 %in% c('0','9','11'))
monkey_hip_mg$Species <- as.factor("Monkey")

# 小胶质图谱的小胶
human_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))
human_mg$Species <- as.factor("Human")

