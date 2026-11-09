setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/"))
rm(list=ls())

# 0. Libraries -----------------------------------------------------------------
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

tmp_obj <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_dims5.qs"))

# metadata--------------------
# Study--------
tmp_obj@meta.data <- tmp_obj@meta.data %>%
  mutate(Study = case_when(
    Study == "Su GSE185553, GSE198323 GSE199243" ~ "Su et al., 2022",
    Study == "Chia GSE303823"                      ~ "Chia et al., 2025",
    Study == "Fullard GSE164485"                   ~ "Fullard et al., 2021",
    Study == "Jakel GSE118257"                     ~ "Jakel et al., 2019",
    Study == "Grubman GSE138852"                   ~ "Grubman et al., 2019",
    TRUE ~ Study
  ))

# 验证结果
unique(tmp_obj$Study)

# Region-----------
tmp_obj@meta.data <- tmp_obj@meta.data %>%
  mutate(Region = case_when(
    Region == "Prefrontal cortex"              ~ "Prefrontal Cortex",
    Region == "PFC"                            ~ "Prefrontal Cortex",
    Region == "Hippocampus"                    ~ "Hippocampus",
    Region == "medulla"                        ~ "Medulla",
    Region == "choroidplexus"                  ~ "Choroid Plexus",
    Region == "White matter"                   ~ "White Matter",
    Region == "Entorhinal Cortex"              ~ "Entorhinal Cortex",
    TRUE ~ Region
  ))

# 验证结果
unique(tmp_obj$Region)

# Condition-----------
tmp_obj@meta.data <- tmp_obj@meta.data %>%
  mutate(Condition = case_when(
    Condition %in% c("CTRL", "Ctrl", "ct", "Healthy") ~ "Control",
    Condition %in% c("AD", "ADD") ~ "AD",
    Condition %in% c("PDD") ~ "PDD",
    Condition %in% c("DLB") ~ "DLB",
    Condition == "MS" ~ "MS",
    Condition %in% c("COVID", "COVID_HIV") ~ "COVID",
    TRUE ~ Condition
  ))

# 验证结果
unique(tmp_obj$Condition)

# Gender-----------
tmp_obj@meta.data <- tmp_obj@meta.data %>%
  mutate(Gender = case_when(
    Gender == "Female" ~ "F",
    Gender == "Male" ~ "M",
    TRUE ~ Gender
  ))

# 验证结果
unique(tmp_obj$Gender)
tmp_obj$Gender <- factor(tmp_obj$Gender, levels = c("F", "M", "Unknown"))

# Age-----------
tmp_obj@meta.data <- tmp_obj@meta.data %>%
  mutate(Age = case_when(
    Age %in% c("0", "1", "2", "3","4","6","13","14","15","16","0.01y") ~ "Child",
    Age %in% c("18", "27","39","50","59","35yr","37yr","44yr","49yr","57yr","32y","48y","50y","56y") ~ "Adult",
    Age %in% c("60","66","71","85","86","88","90","92","95","60yr","64yr","82yr","Aged","60y","64y","68y") ~ "Aged",
    TRUE ~ Age
  ))

# 验证结果
unique(tmp_obj$Age)



##################---------------
levels_condition <- c("Control", "AD", "PDD", "DLB", "MS", "Stroke","COVID")
levels_age       <- c("Child", "Adult", "Aged", "Unknown")
levels_gender    <- c("F", "M", "Unknown")

tmp_obj@meta.data <- tmp_obj@meta.data %>%
  mutate(
    Condition = factor(Condition, levels = intersect(levels_condition, unique(Condition))),
    Age       = factor(Age,       levels = intersect(levels_age, unique(Age))),
    Gender    = factor(Gender,    levels = intersect(levels_gender, unique(Gender))),
    
    Study     = factor(Study),
    Region    = factor(Region)
  )

str(tmp_obj@meta.data[, c("Study", "Region", "Condition", "Gender", "Age")])
qsave(tmp_obj, file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta.qs"))
