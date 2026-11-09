# meta ---------------------------------------------------------------
setwd(Sys.getenv("MG_PROJECT_ROOT"))
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
library(DoubletFinder)
library(qs)

#------------------------------------------------
tmp_obj <- readRDS('data/processed_data/11.integrate_cca_CellType_v1.rds')


#-----------------------------------------------
# — 1.1 给数据加上Age信息
tmp_obj$Sample <- as.character(tmp_obj$sample)
tmp_obj$Age <- as.character(tmp_obj$Sample)

tmp_obj$Age[tmp_obj$Age%in% c('s_18.5ys')]<-"18.5ys"
tmp_obj$Age[tmp_obj$Age%in% c('s_40ys')]<-"40ys"
tmp_obj$Age[tmp_obj$Age%in% c('s_50ys')]<-"50ys"
tmp_obj$Age[tmp_obj$Age%in% c('s_60ys')]<-"60ys"
tmp_obj$Age[tmp_obj$Age%in% c('s_86ys')]<-"86ys"
tmp_obj$Age[tmp_obj$Age%in% c('s_88ys')]<-"88ys"
tmp_obj$Age[tmp_obj$Age%in% c('s_0.6ys')]<-"0.6ys"
tmp_obj$Age[tmp_obj$Age%in% c('s_89_1ys','s_89_2ys')]<-"89ys"
tmp_obj$Age[tmp_obj$Age%in% c('s_89ys')]<-"89ys"
tmp_obj$Age[tmp_obj$Age%in% c('s_73_1ys','s_73_2ys')]<-"73ys"
tmp_obj$Age[tmp_obj$Age%in% c('s_77_1ys','s_77_2ys')]<-"77ys"
tmp_obj$Age[tmp_obj$Age%in% c('s_79_1ys','s_79_2ys')]<-"79ys"
tmp_obj$Age[tmp_obj$Age%in% c('s_1.2ys')]<-"1.2ys"
tmp_obj$Age[tmp_obj$Age%in% c('s_80M_1ys','s_80M_2ys','s_80F_1ys','s_80F_2ys')]<-"80ys"
tmp_obj$Age[tmp_obj$Age%in% c('s_3.2ys')]<-"3.2ys"
tmp_obj$Age[tmp_obj$Age%in% c('l_4D')]<-"0.01ys"
tmp_obj$Age[tmp_obj$Age%in% c('l_32y')]<-"32ys"
tmp_obj$Age[tmp_obj$Age%in% c('l_48y')]<-"48ys"
tmp_obj$Age[tmp_obj$Age%in% c('l_50y')]<-"50ys"
tmp_obj$Age[tmp_obj$Age%in% c('l_56y')]<-"56ys"
tmp_obj$Age[tmp_obj$Age%in% c('l_60y')]<-"60ys"
tmp_obj$Age[tmp_obj$Age%in% c('l_64f','l_64m')]<-"64ys"
tmp_obj$Age[tmp_obj$Age%in% c('l_68y')]<-"68ys"


# — 1.2 给数据加上Stage信息
tmp_obj$Stage <- as.character(tmp_obj$Age)
tmp_obj$Stage[tmp_obj$Stage%in% c("0.01ys","0.6ys","1.2ys","3.2ys")]<-"Child"
tmp_obj$Stage[tmp_obj$Stage%in% c("18.5ys","32ys","40ys","50ys","56ys")]<-"Adult"
tmp_obj$Stage[tmp_obj$Stage%in% c('48ys')]<-"Stroke"
tmp_obj$Stage[tmp_obj$Stage%in% c("60ys","64ys","68ys","73ys","77ys","79ys","80ys","86ys","88ys","89ys")]<-"Aging"


# 定义 Age 的从小到大的顺序
age_levels <- c(
  "0.01ys", "0.6ys", "1.2ys", "3.2ys", 
  "18.5ys", "32ys", "40ys", "48ys", "50ys", "56ys", 
  "60ys", "64ys", "68ys", "73ys", "77ys", "79ys", "80ys", "86ys", "88ys", "89ys"
)

# 定义 Stage 的逻辑顺序
stage_levels <- c("Child", "Adult", "Aging", "Stroke")


tmp_obj$Age <- factor(tmp_obj$Age, levels = age_levels)
tmp_obj$Stage <- factor(tmp_obj$Stage, levels = stage_levels)

#-------------------------------------------------------------------------------
cols_to_remove <- grep("snn_res", colnames(tmp_obj@meta.data), value = TRUE)
cols_to_remove <- c(cols_to_remove, "percent.mito", "sample")
tmp_obj@meta.data <- tmp_obj@meta.data[, !colnames(tmp_obj@meta.data) %in% cols_to_remove]
colnames(tmp_obj@meta.data)

#-------------------------------------------------------------------------------
qsave(tmp_obj, 'data/processed_data/12.integrate_cca_CellType_meta_v1.qs', nthreads = 8)

saveRDS(tmp_obj,'data/processed_data/12.integrate_cca_CellType_meta_v1.rds')


