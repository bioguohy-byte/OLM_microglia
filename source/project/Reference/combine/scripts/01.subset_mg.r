setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/"))
rm(list=ls())

# 0. Libraries -----------------------------------------------------------------
library(Seurat)
library(harmony)
library(qs)
library(dplyr)
library(patchwork)
library(ggplot2)

library(future) # 新增
plan("multicore", workers = 50) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限

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

Yao_v2 <- qread('Yao_et_al./data/processed_data/Yao_v2.qs')
Yao_mg <- subset(Yao_v2,subset = CellType %in% c('Micro'))
qsave(Yao_mg,'Yao_et_al./data/processed_data/Yao_mg.qs')
