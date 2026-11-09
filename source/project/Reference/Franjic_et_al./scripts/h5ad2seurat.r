# 用h5ad的adata.X（counts）转换成seurat的counts---------------
rm(list = ls())
library(Seurat)
# library(qs)
library(reticulate)
library(hdf5r)
library(sceasy)
library(BiocParallel)
register(MulticoreParam(workers = 10, progressbar = TRUE)) 

# 工作目录----------
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference"))

# h5ad data ------------
h5ad_data_path <- "/Franjic_et_al./data/processed_data/human_brain_hippocampus_Franjic_2021_10x/processedData/annot.h5ad"
seurat_data_path <- "/Franjic_et_al./data/processed_data/human_brain_hippocampus_Franjic_2021_10x/processedData/Franjic.rds"

# 在R语言中加载python环境
# 我的conda环境的路径
use_condaenv(Sys.getenv("MG_EXTERNAL_001"))

# Seurat to AnnData
sceasy::convertFormat(h5ad_data_path, 
                      from="anndata", 
                      to="seurat",
                      main_layer = "counts", #X -> counts 默认是把adata.X赋予给seurat中的counts的
                      outFile=seurat_data_path)
