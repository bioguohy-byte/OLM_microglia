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


print("----------DoubletFinder----------------")
seu.obj <- readRDS('data/processed_data/05.no_integrate_pre_cluster.rds')
tmp_obj<-seu.obj
tmp_meta<-tmp_obj@meta.data
all_DL<-data.frame()
summary<-data.frame()
tmp_Plist<-list()
for (i in unique(tmp_obj$orig.ident)) {
  # i<-"NTES"
  subobj<-subset(tmp_obj, subset = orig.ident %in% i)
  DoubletRate = ncol(subobj)*8*1e-6-0.002
  # DoubletRate = 0.04
  
  ## Pre-process Seurat object (standard)
  seu_df <- NormalizeData(subobj)
  seu_df <- FindVariableFeatures(seu_df, selection.method = "vst", nfeatures = 2000)
  seu_df <- ScaleData(seu_df)
  seu_df <- RunPCA(seu_df)
  seu_df <- RunUMAP(seu_df, dims = 1:10)
  
  ## pK Identification (no ground-truth)
  sweep.res.list_kidney <- paramSweep_v3(seu_df, PCs = 1:10, sct = FALSE)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  pK_bcmvn <- bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)] %>% as.character() %>% as.numeric()
  
  ## Homotypic Doublet Proportion Estimate
  annotations <- tmp_meta[tmp_meta$orig.ident%in%i,]$Raw_type
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(DoubletRate*nrow(seu_df@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies
  seu_df <- doubletFinder_v3(seu_df, PCs = 1:10, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
  
  p<-DimPlot(seu_df, reduction = "umap", group.by = grep("^DF",colnames(seu_df@meta.data),value = T))+labs(title = i)
  tmp_Plist[[i]]<-p
  
  DL<-seu_df@meta.data[,c(grep("^DF",colnames(seu_df@meta.data),"orig.ident"),value = T),drop=F] %>%  
    {colnames(.)<-c("if_DL","orig.ident");.}
  
  all_DL<-rbind(all_DL,DL)
  summary<-rbind(summary,data.frame(sample=i,nCells=nrow(seu_df@meta.data),
                                    Raw_Doublet_rate=DoubletRate,Raw_Doublet=nExp_poi,
                                    Homotypic_rate=homotypic.prop,Hetertypic_Doublet=nExp_poi.adj))
  
}

print("----------Save Results----------------")
save(all_DL,summary,tmp_Plist,file = paste("data/processed_data/06.no_integrate_DoubletFinder.result.RData",sep = ""))


# --------------------------------
load("data/processed_data/06.no_integrate_DoubletFinder.result.RData")

summary

Reduce("+",tmp_Plist)+plot_layout(ncol = 4,guides = "collect")


tmp_obj@meta.data %<>% {.$DoubletFinder<-"Singlet";.} %<>% 
  {.[rownames(all_DL[all_DL$if_DL%in%"Doublet",]),]$DoubletFinder<-"Doublet";.}

saveRDS(tmp_obj,file = 'data/processed_data/07.no_integrate_after_df.rds')

p1 <- DimPlot(tmp_obj, reduction = "umap", group.by = "DoubletFinder",pt.size = 0.2,raster = F)
p2 <- DimPlot(tmp_obj, reduction = "umap", group.by = "RNA_snn_res.0.4",pt.size = 0.2,raster = F,label = T,repel = T)
p3 <- DimPlot(tmp_obj, reduction = "umap", group.by = "Raw_type",pt.size = 0.2,raster = F,label = T)

p_df <- p1|p2|p3
ggsave("results/figure/no_integrate_after_df.png",plot = p_df,width = 25,height = 8)
ggsave("results/figure/no_integrate_after_df.pdf",plot = p_df,width = 25,height = 8)
