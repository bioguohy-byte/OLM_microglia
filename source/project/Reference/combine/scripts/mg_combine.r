setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/"))
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


# subset mg------------------
Su_v2 <- qread('Su_et_al./data/processed_data/Su_CSC_alldata_v2.qs')
Su_mg <- subset(Su_v2,subset = CellType %in% c('MG'))
qsave(Su_mg,'Su_et_al./data/processed_data/Su_mg.qs')

Chia_v2 <- qread('Chia_et_al./data/processed_data/Chia_v2.qs')
Chia_mg <- subset(Chia_v2,subset = CellType %in% c('Microglia'))
qsave(Chia_mg,'Chia_et_al./data/processed_data/Chia_mg.qs')

Fullard_v2 <- qread('Fullard_et_al./data/processed_data/Fullard_v2.qs')
Fullard_mg <- subset(Fullard_v2,subset = CellType %in% c('Mic'))
qsave(Fullard_mg,'Fullard_et_al./data/processed_data/Fullard_mg.qs')

Jakel_v2 <- qread('Jakel_et_al./data/processed_data/Jakel_v2.qs')
Jakel_mg <- subset(Jakel_v2,subset = CellType %in% c('Immune_cells','ImOlGs','Microglia_Macrophages'))
qsave(Jakel_mg,'Jakel_et_al./data/processed_data/Jakel_mg.qs')

Grubman_v2 <- qread('Grubman_et_al./data/processed_data/Grubman_v2.qs')
Grubman_mg <- subset(Grubman_v2,subset = CellType %in% c('Microglia'))
qsave(Grubman_mg,'Grubman_et_al./data/processed_data/Grubman_mg.qs')



test1 <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/MG_Integrated_Harmony_Batch.qs"))
DimPlot(test1)



test2 <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/MG_Integrated_CCA_Study.qs"))
DimPlot(test2)
FeaturePlot(test2,features = 'PLP1')
DefaultAssay(test2) <- 'RNA'

test3 <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/MG_Integrated_CCA_Study_2.qs"))
DimPlot(test3)
FeaturePlot(test3,features = 'PLP1')
DefaultAssay(test3) <- 'RNA'