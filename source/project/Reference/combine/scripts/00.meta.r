setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/"))
rm(list=ls())

# # 多线程运行 -------------------
# library(future)
# plan(multisession, workers = 15)
# nbrOfWorkers()
# options(future.globals.maxSize = 1999999 * 1923^2)

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
library(qs)

# process metadata------------
# — Su_et_al.----------------
Su <- qread('Su_et_al./data/processed_data/Su_CSC_alldata.qs')

Su_v2 <- Su
Su_v2@meta.data <- Su@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")]

Su_v2$Batch <- as.factor(Su$Batch)
Su_v2$Sample <- as.factor(Su$Sample)
Su_v2$Condition <- as.factor("Healthy")
Su_v2$Gender <- as.factor(Su$Sex)
Su_v2$Age <- as.factor(Su$Age)
Su_v2$Region <- as.factor("Hippocampus")
Su_v2$CellType <- as.factor(Su$MajorCellTypes)
Su_v2$Study <- as.factor("Su GSE185553, GSE198323 GSE199243")
qsave(Su_v2,'Su_et_al./data/processed_data/Su_CSC_alldata_v2.qs')

# — Chia_et_al.----------------
Chia <- qread('Chia_et_al./data/processed_data/Chia.qs')

Chia_v2 <- Chia
Chia_v2@meta.data <- Chia@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")]

Chia_v2$Batch <- as.factor(Chia$sample)
Chia_v2$Sample <- as.factor(Chia$sample)
Chia_v2$Condition <- as.factor(Chia$group)
Chia_v2$Gender <- as.factor(Chia$sex)
Chia_v2$Age <- as.factor('Aged')
Chia_v2$Region <- as.factor('Prefrontal cortex')
Chia_v2$CellType <- as.factor(Chia$cell.type)
Chia_v2$Study <- as.factor("Chia GSE303823")
qsave(Chia_v2,'Chia_et_al./data/processed_data/Chia_v2.qs')


# — Fullard_et_al.----------------
Fullard <- qread('Fullard_et_al./data/processed_data/Fullard.qs')

Fullard_v2 <- Fullard
Fullard_v2@meta.data <- Fullard@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")]

Fullard_v2$Batch <- as.factor(Fullard$pid)
Fullard_v2$Sample <- as.factor(Fullard$subject_id2)
Fullard_v2$Condition <- as.factor(Fullard$clinical)
Fullard_v2$Gender <- as.factor('Unknown')
Fullard_v2$Age <- as.factor('Unknown')
Fullard_v2$Region <- as.factor(Fullard$tissue)
Fullard_v2$CellType <- as.factor(Fullard$celltype)
Fullard_v2$Study <- as.factor("Fullard GSE164485")
qsave(Fullard_v2,'Fullard_et_al./data/processed_data/Fullard_v2.qs')

# — Jakel_et_al.----------------
Jakel <- qread('Jakel_et_al./data/processed_data/Jakel.qs')

Jakel_v2 <- Jakel
Jakel_v2@meta.data <- Jakel@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")]

Jakel_v2$Batch <- as.factor(Jakel$sample)
Jakel_v2$Sample <- as.factor(Jakel$donor_ID)
Jakel_v2$Condition <- as.factor(Jakel$Condition)
Jakel_v2$Gender <- as.factor(Jakel$donor_gender)
Jakel_v2$Age <- as.factor(Jakel$donor_age)
Jakel_v2$Region <- as.factor(Jakel$region)
Jakel_v2$CellType <- as.factor(Jakel$Celltypes)
Jakel_v2$Study <- as.factor("Jakel GSE118257")
qsave(Jakel_v2,'Jakel_et_al./data/processed_data/Jakel_v2.qs')

# — Grubman_et_al.----------------
Grubman <- qread('Grubman_et_al./data/processed_data/Grubman.qs')

Grubman_v2 <- Grubman
Grubman_v2@meta.data <- Grubman@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")]

Grubman_v2$Batch <- as.factor(Grubman$batch)
Grubman_v2$Sample <- as.factor(Grubman$patient)
Grubman_v2$Condition <- as.factor(Grubman$batchCond)
Grubman_v2$Gender <- as.factor(Grubman$sex)
Grubman_v2$Age <- as.factor("Unknown")
Grubman_v2$Region <- as.factor("Entorhinal Cortex")
Grubman_v2$CellType <- as.factor(Grubman$CellType)
Grubman_v2$Study <- as.factor("Grubman GSE138852")
qsave(Grubman_v2,'Grubman_et_al./data/processed_data/Grubman_v2.qs')


# — Yao_et_al.----------------
Yao <- qread('Yao_et_al./data/processed_data/Seuobj.CellType.v1.qs')
Yao_v2 <- Yao
Yao_v2@meta.data <- Yao@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt")]

Yao_v2$Batch <- as.factor(Yao$Batch)
Yao_v2$Sample <- as.factor(Yao$Sample)
Yao_v2$Condition <- as.factor(Yao$Condition)
Yao_v2$Gender <- as.factor(Yao$Gender)
Yao_v2$Age <- Yao$Age
Yao_v2$Region <- as.factor(Yao$Region)
Yao_v2$CellType <- as.factor(Yao$CellType)
Yao_v2$Study <- as.factor("Yao et al., 2024")
qsave(Yao_v2,'Yao_et_al./data/processed_data/Yao_v2.qs')
