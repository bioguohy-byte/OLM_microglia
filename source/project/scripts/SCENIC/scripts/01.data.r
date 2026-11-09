##文件准备
# pyscenic的输入文件是行为基因名，列为细胞ID的矩阵，所以在seurat对象中导出矩阵的时候需要转置一下，可以用标准化矩阵，也可以用counts矩阵，影响不大！
#表达矩阵、meta----R中进行

setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "scripts/SCENIC/"))
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

#加载数据
# 小胶质细胞
tmp_obj <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "data/processed_data/14.mg_combined_meta.qs"))

write.csv(t(as.matrix(tmp_obj@assays$RNA@counts)), file = "data/tmp_obj.csv") #其中前面的t就是转置的意思


