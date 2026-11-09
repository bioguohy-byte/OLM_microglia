setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Yao_et_al./"))
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
library(DoubletFinder)

# 多线程运行
# # 终端运行时选择 multicore
library(future) # 新增
plan("multicore", workers = 50) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限

# # rstudio server运行时选择 multisession
# library(future)
# plan(multisession, workers = 6)
# nbrOfWorkers()
# options(future.globals.maxSize = 1999999 * 1924^2)

# 1. 加载数据 -------------------------
l_4D <- Read10X_h5("data/raw_data/4D_filtered_feature_bc_matrix.h5", use.names = T)
l_32y <- Read10X_h5("data/raw_data/32Y_filtered_feature_bc_matrix.h5", use.names = T)
l_48y <- Read10X_h5("data/raw_data/48Y_filtered_feature_bc_matrix.h5", use.names = T)
l_50y <- Read10X_h5("data/raw_data/50Y_filtered_feature_bc_matrix.h5", use.names = T)
l_56y <- Read10X_h5("data/raw_data/56Y_filtered_feature_bc_matrix.h5", use.names = T)
l_60y <- Read10X_h5("data/raw_data/60Y_filtered_feature_bc_matrix.h5", use.names = T)
l_64f <- Read10X_h5("data/raw_data/64F_filtered_feature_bc_matrix.h5", use.names = T)
l_64m <- Read10X_h5("data/raw_data/64M_filtered_feature_bc_matrix.h5", use.names = T)
l_68y <- Read10X_h5("data/raw_data/68Y_filtered_feature_bc_matrix.h5", use.names = T)

# 2. 创建seurat对象 -----------------
l_4D <- CreateSeuratObject(l_4D, project = "l_4D", min.cells = 3, min.features = 200)
l_48y <- CreateSeuratObject(l_48y, project = "l_48y", min.cells = 3, min.features = 200)
l_32y <- CreateSeuratObject(l_32y, project = "l_32y", min.cells = 3, min.features = 200)
l_50y <- CreateSeuratObject(l_50y, project = "l_50y", min.cells = 3, min.features = 200)
l_56y <- CreateSeuratObject(l_56y, project = "l_56y", min.cells = 3, min.features = 200)
l_60y <- CreateSeuratObject(l_60y, project = "l_60y", min.cells = 3, min.features = 200)
l_64f <- CreateSeuratObject(l_64f, project = "l_64f", min.cells = 3, min.features = 200)
l_64m <- CreateSeuratObject(l_64m, project = "l_64m", min.cells = 3, min.features = 200)
l_68y <- CreateSeuratObject(l_68y, project = "l_68y", min.cells = 3, min.features = 200)

alldata <- merge(x=l_4D,
                 y=c(l_32y,l_48y,l_50y,l_56y,l_60y,l_64f,l_64m,l_68y),
                 add.cell.ids=c("l_4D","l_32y","l_48y","l_50y","l_56y",
                                "l_60y","l_64f","l_64m","l_68y"))

# 3. 添加metadata信息----------------------------------
meta <- read.csv("meta/sample_information.csv",
                 header = TRUE,
                 check.names = FALSE,
                 sep = ",")

tmp_meta<-alldata@meta.data %>% rownames_to_column() %>%
  left_join(.,meta,by=c("orig.ident"="Sample")) %>%
  column_to_rownames(var = "rowname") %>%
  dplyr::select(.,!1:3)

tmp_obj <- alldata %>%
  AddMetaData(.,tmp_meta) %>%
  NormalizeData() %>%
  FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(., features = VariableFeatures(.), npcs = 30, verbose = T) %>%
  FindNeighbors(.,dims = 1:30) %>%
  FindClusters(.,resolution = seq(0.1,1,0.1)) %>%
  RunUMAP(., seed.use = 42, reduction = "pca", dims = 1:30)

tmp_obj$Sample <- tmp_obj$orig.ident

# qsave(tmp_obj,"data/processed_data/Raw_data.for_check_batch_effect.qs")


# 4. 绘制meta umap------------------
tmp_obj <- qread('data/processed_data/Raw_data.for_check_batch_effect.qs')

(p1<-DimPlot(tmp_obj,group.by = c("orig.ident","Batch","Sample","Condition","Age","Gender","Region"),label = T))
(p<-DimPlot(tmp_obj,group.by = "Batch",label = T,repel = T)+NoLegend())
(p2<-DimPlot(tmp_obj,group.by = "Batch",split.by = "Batch",ncol = 6))
(p|p2)+plot_layout(widths = c(1,2)) #右图的宽度是左图的 2 倍

tmp_Idents<-"RNA_snn_res.0.2"
tmp_label=13
(p1<-DimPlot(tmp_obj, reduction = "umap",label.size = 5,
             group.by = tmp_Idents,
             # cols = Major_color[levels(tmp_Sobj$Major),]$Color,
             raster = F,label = T, pt.size = 0.1,repel = T)+
    guides(color = guide_legend(ncol = 1,override.aes = list(size = 2.5)))+
    tidydr::theme_dr(xlength = 0.2,ylength = 0.2,
                     arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
    theme(panel.grid = element_blank(),plot.title = element_blank(),
          legend.position = "none",
          axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
          axis.title.y.right = element_text(size = tmp_label,hjust = 0.5,face = "bold")))

# 5. 绘制marker dotplot--------------------
# 定义各群marker
(tmp_Marker<-data.frame(Gname=c("AQP4","ITGB4","ETNPPL","SLC1A2","RFX4","GLI3","GRAMD2B","GFAP",
                                "VCAN","C1QL1","SUSD5","COL9A1",
                                "MOG","MOBP","SLC5A11","PICK1","HAPLN2","CNDP1","MYRF","CPM","TF","TMC6","CDK18",
                                "PTPRC", # Immune
                                "CSF1R","MRC1","IL18","C3","SLC2A5","BIN2","DOCK2","ENTPD1",
                                "CFAP54","DNAH11",
                                "DCN","LAMA2","RERGL","ABCG2","EMCN","LEF1","FLT1","SLC6A12","CFH","SLCO1B3",
                                "RBFOX3","SNAP25","SYT1","GRIN1" # Neurons


),
Type=c(rep("Astr",8),
       rep("OPC",4),
       rep("Oli",11),
       rep("I",1),
       rep("Micro",8),
       rep("Ependymal",2),
       rep("Vascular",10),
       rep("Neurons",4)

)))



myBubblePlot <- function(tmp_obj,tmp_Marker,tmp_Scale) {
  tmp_Plist<-list()
  tmp_width<-c()
  tmp_df<-DotPlot(object = tmp_obj, features = unique(tmp_Marker$Gname),group.by  = tmp_Idents,scale = tmp_Scale)
  for (tmp_n in 1:length(unique(tmp_Marker$Type))) {
    # tmp_n<-1
    tmpG<-tmp_Marker[tmp_Marker$Type%in%unique(tmp_Marker$Type)[tmp_n],]$Gname %>% unique() %>%
      intersect(.,rownames(tmp_obj))
    tmp_width<-c(tmp_width,length(tmpG))
    pdata<-DotPlot(object = tmp_obj, features = tmpG,group.by  = tmp_Idents,scale = tmp_Scale) %>%
      {left_join(.$data,tmp_Marker,by=c("features.plot"="Gname"))} %>% .[.$Type%in%unique(tmp_Marker$Type)[tmp_n],] %>%
      {.$features.plot<-factor(.$features.plot,levels = unique(tmp_Marker$Gname));.}
    if (tmp_n==1) {
      p<-ggplot(pdata,aes(x = features.plot,y=id,size=pct.exp,colour=avg.exp.scaled,ylab=''))+
        geom_point()+
        scale_size_continuous(range=c(1,4),limits = c(0,100))+
        scale_color_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","grey90","white")),
                              limits = c(min(tmp_df$data$avg.exp.scaled),max(tmp_df$data$avg.exp.scaled)))+
        facet_wrap(.~Type)+
        guides(colour = guide_colourbar(title = "Average Expression"))+
        guides(size = guide_legend(title = "Percent Expressed"))+
        theme_classic()+
        theme(panel.grid = element_blank(),panel.background = element_blank(),
              # panel.spacing = unit(0.1, "lines"),
              axis.title = element_blank(),
              axis.text = element_text(size = tmp_label,colour = "black"),
              axis.text.x = element_text(angle = 45,hjust = 1),
              legend.text = element_text (size = tmp_label),
              legend.title = element_text (size = tmp_label),
              strip.background = element_rect(color="grey90",fill = "grey90"),
              strip.text = element_text(size = tmp_label,colour = "black"),
              strip.placement = "outside",
              strip.switch.pad.grid = unit(2, "mm"))
      tmp_Plist[[tmp_n]]<-p
    } else {
      p<-ggplot(pdata,aes(x = features.plot,y=id,size=pct.exp,colour=avg.exp.scaled,ylab=''))+
        geom_point()+
        scale_size_continuous(range=c(1,4),limits = c(0,100))+
        scale_color_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","grey90","white")),
                              limits = c(min(tmp_df$data$avg.exp.scaled),max(tmp_df$data$avg.exp.scaled)))+
        facet_wrap(.~Type)+
        guides(colour = guide_colourbar(title = "Average Expression"))+
        guides(size = guide_legend(title = "Percent Expressed"))+
        theme_classic()+
        theme(panel.grid = element_blank(),panel.background = element_blank(),
              # panel.spacing = unit(0.1, "lines"),
              axis.title = element_blank(),
              axis.title.y = element_blank(),axis.line.y = element_blank(),
              axis.text.y = element_blank(),axis.ticks.y = element_blank(),
              axis.text = element_text(size = tmp_label,colour = "black"),
              axis.text.x = element_text(angle = 45,hjust = 1),
              legend.text = element_text (size = tmp_label),
              legend.title = element_text (size = tmp_label),
              strip.background = element_rect(color="grey90",fill = "grey90"),
              strip.text = element_text(size = tmp_label,colour = "black"),
              strip.placement = "outside",
              strip.switch.pad.grid = unit(2, "mm"))
      tmp_Plist[[tmp_n]]<-p
    }
  }
  p_out<-Reduce("+",tmp_Plist)+plot_layout(ncol = length(tmp_Plist),
                                           widths = tmp_width,guides = "collect")
  return(p_out)
}
(p2<-tmp_obj %>%
    # {.$Major<-factor(.$Major,levels = rev(levels(.$Major)));.} %>%
    myBubblePlot(tmp_obj=.,tmp_Marker=tmp_Marker,tmp_Scale=T))

(pp1<-(as.ggplot(p1)|as.ggplot(p2))+plot_layout(widths = c(1,2.5)))
ggsave2("results/figure/res.0.2_marker.pdf",plot = pp1,width = 20,height = 7)
ggsave2("results/figure/res.0.2_marker.png",plot = pp1,width = 20,height = 7)

# 6. 粗定义细胞类型-------------------
tmp_obj$Raw_type<-as.character(tmp_obj$RNA_snn_res.0.2)

tmp_obj@meta.data %<>%
  mutate(Raw_type=case_when(RNA_snn_res.0.2%in%c("2","4","9","10","25")~"Astro",
                            RNA_snn_res.0.2%in%c("5","8","20")~"OPC",
                            RNA_snn_res.0.2%in%c("0")~"Oli",
                            RNA_snn_res.0.2%in%c("3","21")~"Micro",
                            RNA_snn_res.0.2%in%c("13","23")~"Ependymal",
                            RNA_snn_res.0.2%in%c("12")~"Vasular",
                            RNA_snn_res.0.2%in%c("1","6","7","11","14","15","16","17","18","19","22","24")~"Neurons"
                            ))

(p1<-DimPlot(tmp_obj,group.by = "Raw_type",label = T))
qsave(tmp_obj,"data/processed_data/Raw_data.raw_type.qs")

# 7. DoubletFinder------------------------------------
tmp_obj <- qread("data/processed_data/Raw_data.raw_type.qs")

raw_type <- "Raw_type"
tmp_Idents<-"RNA_snn_res.0.2"
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
  seu_obj <- NormalizeData(subobj)
  seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
  seu_obj <- ScaleData(seu_obj)
  seu_obj <- RunPCA(seu_obj)
  seu_obj <- RunUMAP(seu_obj, dims = 1:10)

  ## pK Identification (no ground-truth)
  sweep.res.list <- paramSweep_v3(seu_obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

  ## Homotypic Doublet Proportion Estimate
  annotations <- tmp_meta[tmp_meta$orig.ident%in%i,][[raw_type]]
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(DoubletRate*nrow(seu_obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  ## Run DoubletFinder with varying classification stringencies
  seu_obj <- doubletFinder_v3(seu_obj, PCs = 1:10, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)

  p<-DimPlot(seu_obj, reduction = "umap", group.by = grep("^DF",colnames(seu_obj@meta.data),value = T))+labs(title = i)
  tmp_Plist[[i]]<-p

  DL<-seu_obj@meta.data[,c(grep("^DF",colnames(seu_obj@meta.data),"orig.ident"),value = T),drop=F] %>%
    {colnames(.)<-c("if_DL","orig.ident");.}

  all_DL<-rbind(all_DL,DL)
  summary<-rbind(summary,data.frame(sample=i,nCells=nrow(seu_obj@meta.data),
                                    Raw_Doublet_rate=DoubletRate,Raw_Doublet=nExp_poi,
                                    Homotypic_rate=homotypic.prop,Hetertypic_Doublet=nExp_poi.adj))

}
qsavem(all_DL, summary, tmp_Plist, file = "data/processed_data/DoubletFinder.result.qs")

Reduce("+",tmp_Plist)+plot_layout(ncol = 4,guides = "collect")

tmp_obj@meta.data %<>% {.$DoubletFinder<-"Singlet";.} %<>%
  {.[rownames(all_DL[all_DL$if_DL%in%"Doublet",]),]$DoubletFinder<-"Doublet";.}

qsave(tmp_obj,"data/processed_data/Raw_data.DoubletFinder.qs")

# 8. 第一轮质控--------------------

tmp_obj <- qread('data/processed_data/Raw_data.DoubletFinder.qs')
alldata<-tmp_obj
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^MT-")

# step1
data.filt <- alldata %>% subset(.,subset=DoubletFinder%in%"Singlet")

# step2
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.mt>10) & orig.ident == "l_4D")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.mt>10) & orig.ident == "l_48y")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.mt>10) & orig.ident == "l_50y")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.mt>10) & orig.ident == "l_56y")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.mt>10) & orig.ident == "l_60y")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.mt>10) & orig.ident == "l_64f")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.mt>10) & orig.ident == "l_64m")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.mt>10) & orig.ident == "l_68y")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.mt>10) & orig.ident == "l_32y")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))

VlnPlot(data.filt,features = c('nFeature_RNA','nCount_RNA','percent.mt'),group.by = 'orig.ident')

qsave(data.filt,file = "data/processed_data/Data_after_filtered.v1.qs")

# 9. 单细胞整合------------------
# — 9.1 CCA 整合(dims50)----------
data.filt <- qread("data/processed_data/Data_after_filtered.v1.qs")
tmp_meta<-data.filt@meta.data %>%
  .[,setdiff(colnames(.),grep("RNA_snn_res.",colnames(.),value = T))] %>%
  dplyr::select(.,-seurat_clusters)

tmp_obj<-data.filt %>%
  {.@assays[["RNA"]]@counts} %>%
  CreateSeuratObject(.,min.features=200,min.cell=3,meta.data = tmp_meta)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

Batch_name <- "orig.ident"
tmp_obj.list <- SplitObject(tmp_obj, split.by = Batch_name)
tmp_obj.list <- lapply(X = tmp_obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = tmp_obj.list,nfeatures = 2000)
seurat.anchors <- FindIntegrationAnchors(object.list = tmp_obj.list,anchor.features = features)
tmp_Sobj <- IntegrateData(anchorset = seurat.anchors)
DefaultAssay(tmp_Sobj) <- "RNA"
tmp_Sobj <- ScaleData(tmp_Sobj,verbose = T)
tmp_Sobj <- CellCycleScoring(tmp_Sobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(tmp_Sobj) <- "integrated"
tmp_Sobj <- ScaleData(tmp_Sobj,verbose = T)
tmp_Sobj_pca <- RunPCA(tmp_Sobj, npcs = 100, verbose = T)
qsave(tmp_Sobj_pca,"data/processed_data/tmp_Sobj.after_PCA.v1.qs")

tmp_Sobj <- RunUMAP(tmp_Sobj_pca, seed.use = 42, reduction = "pca", dims = 1:70) %>%
  FindNeighbors(., reduction = "pca",dims = 1:70) %>%
  FindClusters(., resolution = seq(0.1,1,0.1))
qsave(tmp_Sobj,"data/processed_data/tmp_Sobj.after_FindClusters.v1.qs")

tmp_Plist<-list()
for (i in grep("integrated_snn_res",colnames(tmp_Sobj@meta.data),value = T)[1:10]) {
  title<-gsub(".*res.","",i) %>% paste("Resolution ",.,sep = "")
  p<-DimPlot(tmp_Sobj, reduction = "umap", label = TRUE,group.by = i,raster = T) + labs(title = title) + NoLegend()
  tmp_Plist[[as.character(i)]]<-p
}
(p1<-Reduce("+",tmp_Plist)+plot_layout(ncol = 5))
ggsave(paste("results/figure/tmp_Sobj.Resolution.Umap.1.v1.pdf",sep = ""),p1,width = 18,height = 9)

# tmp_Plist<-list()
# for (i in grep("integrated_snn_res",colnames(tmp_Sobj@meta.data),value = T)[11:20]) {
#   title<-gsub(".*res.","",i) %>% paste("Resolution ",.,sep = "")
#   p<-DimPlot(tmp_Sobj, reduction = "umap", label = TRUE,group.by = i,raster = T) + labs(title = title) + NoLegend()
#   tmp_Plist[[as.character(i)]]<-p
# }
# (p2<-Reduce("+",tmp_Plist)+plot_layout(ncol = 5))
# ggsave(paste("results/figure/tmp_Sobj.Resolution.Umap.2.v1.pdf",sep = ""),p2,width = 18,height = 9)

# — 9.2 RPCA 整合-------------------
data.filt <- qread("data/processed_data/Data_after_filtered.v1.qs")

# 准备 Meta data 和创建 Seurat 对象
tmp_meta <- data.filt@meta.data %>%
  .[,setdiff(colnames(.), grep("RNA_snn_res.", colnames(.), value = T))] %>%
  dplyr::select(., -seurat_clusters)

tmp_obj <- data.filt %>%
  {.@assays[["RNA"]]@counts} %>%
  CreateSeuratObject(., min.features = 200, min.cell = 3, meta.data = tmp_meta)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
Batch_name <- "orig.ident"

# 拆分对象
tmp_obj.list <- SplitObject(tmp_obj, split.by = Batch_name)

# 标准化和寻找高变基因
tmp_obj.list <- lapply(X = tmp_obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# 选择整合特征
features <- SelectIntegrationFeatures(object.list = tmp_obj.list, nfeatures = 2000)

# RPCA 核心步骤：整合前需对每个对象单独 Scale 和 PCA
tmp_obj.list <- lapply(X = tmp_obj.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# 使用 RPCA 寻找锚点
seurat.anchors <- FindIntegrationAnchors(object.list = tmp_obj.list,
                                         anchor.features = features,
                                         reduction = "rpca")

tmp_Sobj <- IntegrateData(anchorset = seurat.anchors)

# 切换回 RNA assay 进行细胞周期评分
DefaultAssay(tmp_Sobj) <- "RNA"
tmp_Sobj <- ScaleData(tmp_Sobj, verbose = T)
tmp_Sobj <- CellCycleScoring(tmp_Sobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# 切换到 integrated assay 进行后续降维
DefaultAssay(tmp_Sobj) <- "integrated"
tmp_Sobj <- ScaleData(tmp_Sobj, verbose = T)
tmp_Sobj_pca <- RunPCA(tmp_Sobj, npcs = 100, verbose = T)
qsave(tmp_Sobj_pca, "data/processed_data/tmp_Sobj.rpca.after_PCA.v1.qs")

# UMAP 和 聚类
tmp_Sobj <- RunUMAP(tmp_Sobj_pca, seed.use = 42, reduction = "pca", dims = 1:70) %>%
  FindNeighbors(., reduction = "pca",dims = 1:70) %>%
  FindClusters(., resolution = seq(0.1, 2, 0.1))
qsave(tmp_Sobj, "data/processed_data/tmp_Sobj.rpca.after_FindClusters.v1.qs")

# 绘图部分 1
tmp_Plist <- list()
for (i in grep("integrated_snn_res", colnames(tmp_Sobj@meta.data), value = T)[1:10]) {
  title <- gsub(".*res.", "", i) %>% paste("Resolution ", ., sep = "")
  p <- DimPlot(tmp_Sobj, reduction = "umap", label = TRUE, group.by = i, raster = T) +
    labs(title = title) + NoLegend()
  tmp_Plist[[as.character(i)]] <- p
}
(p1 <- Reduce("+", tmp_Plist) + plot_layout(ncol = 5))
ggsave("results/figure/tmp_Sobj.rpca.Resolution.Umap.1.v1.pdf", p1, width = 18, height = 9)

# 绘图部分 2
tmp_Plist <- list()
for (i in grep("integrated_snn_res", colnames(tmp_Sobj@meta.data), value = T)[11:20]) {
  title <- gsub(".*res.", "", i) %>% paste("Resolution ", ., sep = "")
  p <- DimPlot(tmp_Sobj, reduction = "umap", label = TRUE, group.by = i, raster = T) +
    labs(title = title) + NoLegend()
  tmp_Plist[[as.character(i)]] <- p
}
(p2 <- Reduce("+", tmp_Plist) + plot_layout(ncol = 5))

# — 9.3 Harmony 整合--------------
library(harmony)

data.filt <- qread("data/processed_data/Data_after_filtered.v1.qs")

# 准备 Meta data
tmp_meta <- data.filt@meta.data %>%
  .[,setdiff(colnames(.), grep("RNA_snn_res.", colnames(.), value = T))] %>%
  dplyr::select(., -seurat_clusters)

# 创建 Seurat 对象
tmp_obj <- data.filt %>%
  {.@assays[["RNA"]]@counts} %>%
  CreateSeuratObject(., min.features = 200, min.cell = 3, meta.data = tmp_meta)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
Batch_name <- "orig.ident"

# 1. 筛选高质量整合特征
tmp_obj.list <- SplitObject(tmp_obj, split.by = Batch_name)
tmp_obj.list <- lapply(X = tmp_obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = tmp_obj.list, nfeatures = 2000)

# 2. 数据处理与 PCA
VariableFeatures(tmp_obj) <- features
tmp_obj <- NormalizeData(tmp_obj)
tmp_obj <- ScaleData(tmp_obj, verbose = T)

# 细胞周期评分
tmp_obj <- CellCycleScoring(tmp_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# 再次 Scale (基于评分后的数据) 并运行 PCA
tmp_obj <- ScaleData(tmp_obj, verbose = T)
tmp_Sobj_pca <- RunPCA(tmp_obj, npcs = 30, verbose = T)

# 3. 运行 Harmony
tmp_Sobj_pca <- RunHarmony(tmp_Sobj_pca, group.by.vars = Batch_name)

# 保存 PCA 结果 (文件名包含 harmony)
qsave(tmp_Sobj_pca, "data/processed_data/tmp_Sobj.harmony.after_PCA.v1.qs")

# 4. UMAP 和 聚类 (使用 harmony reduction)
tmp_Sobj <- RunUMAP(tmp_Sobj_pca, reduction = "harmony", dims = 1:30, seed.use = 42) %>%
  FindNeighbors(., reduction = "harmony", dims = 1:30) %>%
  FindClusters(., resolution = seq(0.1, 2, 0.1))

# 保存聚类结果 (文件名包含 harmony)
qsave(tmp_Sobj, "data/processed_data/tmp_Sobj.harmony.after_FindClusters.v1.qs")

# 5. 绘图 (注意：Harmony 模式下聚类列名通常是 RNA_snn_res 或 snn_res)
# 这里使用 grep 匹配 "snn_res" 以兼容可能的命名差异
cluster_cols <- grep("_snn_res", colnames(tmp_Sobj@meta.data), value = T)

tmp_Plist <- list()
for (i in cluster_cols[1:10]) {
  title <- gsub(".*res.", "", i) %>% paste("Resolution ", ., sep = "")
  p <- DimPlot(tmp_Sobj, reduction = "umap", label = TRUE, group.by = i, raster = T) +
    labs(title = title) + NoLegend()
  tmp_Plist[[as.character(i)]] <- p
}
(p1 <- Reduce("+", tmp_Plist) + plot_layout(ncol = 5))
# 保存图片 1 (文件名包含 harmony)
ggsave("results/figure/tmp_Sobj.harmony.Resolution.Umap.1.v1.pdf", p1, width = 18, height = 9)

if(length(cluster_cols) > 10) {
  tmp_Plist <- list()
  for (i in cluster_cols[11:min(20, length(cluster_cols))]) {
    title <- gsub(".*res.", "", i) %>% paste("Resolution ", ., sep = "")
    p <- DimPlot(tmp_Sobj, reduction = "umap", label = TRUE, group.by = i, raster = T) +
      labs(title = title) + NoLegend()
    tmp_Plist[[as.character(i)]] <- p
  }
  (p2 <- Reduce("+", tmp_Plist) + plot_layout(ncol = 5))
  # 保存图片 2 (文件名包含 harmony)
  ggsave("results/figure/tmp_Sobj.harmony.Resolution.Umap.2.v1.pdf", p2, width = 18, height = 9)
}

# 10. 整合后marker dotplot------------------
tmp_obj <- qread('data/processed_data/tmp_Sobj.rpca.after_FindClusters.v1.qs')
DefaultAssay(tmp_obj) <- 'RNA'
tmp_Idents<-"integrated_snn_res.0.5"
tmp_label=13

# 调整顺序
current_idents <- tmp_obj@meta.data[[tmp_Idents]]
proper_order <- sort(unique(as.numeric(as.character(current_idents))))
tmp_obj@meta.data[[tmp_Idents]] <- factor(current_idents, levels = proper_order)

# 定义各群marker
(tmp_Marker<-data.frame(Gname=c("AQP4","ITGB4","ETNPPL","SLC1A2","RFX4","GLI3","GRAMD2B","GFAP",
                                "VCAN","C1QL1","SUSD5","COL9A1",
                                "MOG","MOBP","SLC5A11","PICK1","HAPLN2","CNDP1","MYRF","CPM","TF","TMC6","CDK18",
                                "PTPRC", # Immune
                                "CSF1R","MRC1","IL18","C3","SLC2A5","BIN2","DOCK2","ENTPD1",
                                "CFAP54","DNAH11",
                                "DCN","LAMA2","RERGL","ABCG2","EMCN","LEF1","FLT1","SLC6A12","CFH","SLCO1B3",
                                "RBFOX3","SNAP25","SYT1","GRIN1", # Neurons
                                "CD247" # T cells
                                
                                
),
Type=c(rep("Astr",8),
       rep("OPC",4),
       rep("Oli",11),
       rep("I",1),
       rep("Micro",8),
       rep("Ependymal",2),
       rep("Vascular",10),
       rep("Neurons",4),
       rep("T cells",1)
       
)))



myBubblePlot <- function(tmp_obj,tmp_Marker,tmp_Scale) {
  tmp_Plist<-list()
  tmp_width<-c()
  tmp_df<-DotPlot(object = tmp_obj, features = unique(tmp_Marker$Gname),group.by  = tmp_Idents,scale = tmp_Scale)
  for (tmp_n in 1:length(unique(tmp_Marker$Type))) {
    # tmp_n<-1
    tmpG<-tmp_Marker[tmp_Marker$Type%in%unique(tmp_Marker$Type)[tmp_n],]$Gname %>% unique() %>%
      intersect(.,rownames(tmp_obj))
    tmp_width<-c(tmp_width,length(tmpG))
    pdata<-DotPlot(object = tmp_obj, features = tmpG,group.by  = tmp_Idents,scale = tmp_Scale) %>%
      {left_join(.$data,tmp_Marker,by=c("features.plot"="Gname"))} %>% .[.$Type%in%unique(tmp_Marker$Type)[tmp_n],] %>%
      {.$features.plot<-factor(.$features.plot,levels = unique(tmp_Marker$Gname));.}
    if (tmp_n==1) {
      p<-ggplot(pdata,aes(x = features.plot,y=id,size=pct.exp,colour=avg.exp.scaled,ylab=''))+
        geom_point()+
        scale_size_continuous(range=c(1,4),limits = c(0,100))+
        scale_color_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","grey90","white")),
                              limits = c(min(tmp_df$data$avg.exp.scaled),max(tmp_df$data$avg.exp.scaled)))+
        facet_wrap(.~Type)+
        guides(colour = guide_colourbar(title = "Average Expression"))+
        guides(size = guide_legend(title = "Percent Expressed"))+
        theme_classic()+
        theme(panel.grid = element_blank(),panel.background = element_blank(),
              # panel.spacing = unit(0.1, "lines"),
              axis.title = element_blank(),
              axis.text = element_text(size = tmp_label,colour = "black"),
              axis.text.x = element_text(angle = 45,hjust = 1),
              legend.text = element_text (size = tmp_label),
              legend.title = element_text (size = tmp_label),
              strip.background = element_rect(color="grey90",fill = "grey90"),
              strip.text = element_text(size = tmp_label,colour = "black"),
              strip.placement = "outside",
              strip.switch.pad.grid = unit(2, "mm"))
      tmp_Plist[[tmp_n]]<-p
    } else {
      p<-ggplot(pdata,aes(x = features.plot,y=id,size=pct.exp,colour=avg.exp.scaled,ylab=''))+
        geom_point()+
        scale_size_continuous(range=c(1,4),limits = c(0,100))+
        scale_color_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","grey90","white")),
                              limits = c(min(tmp_df$data$avg.exp.scaled),max(tmp_df$data$avg.exp.scaled)))+
        facet_wrap(.~Type)+
        guides(colour = guide_colourbar(title = "Average Expression"))+
        guides(size = guide_legend(title = "Percent Expressed"))+
        theme_classic()+
        theme(panel.grid = element_blank(),panel.background = element_blank(),
              # panel.spacing = unit(0.1, "lines"),
              axis.title = element_blank(),
              axis.title.y = element_blank(),axis.line.y = element_blank(),
              axis.text.y = element_blank(),axis.ticks.y = element_blank(),
              axis.text = element_text(size = tmp_label,colour = "black"),
              axis.text.x = element_text(angle = 45,hjust = 1),
              legend.text = element_text (size = tmp_label),
              legend.title = element_text (size = tmp_label),
              strip.background = element_rect(color="grey90",fill = "grey90"),
              strip.text = element_text(size = tmp_label,colour = "black"),
              strip.placement = "outside",
              strip.switch.pad.grid = unit(2, "mm"))
      tmp_Plist[[tmp_n]]<-p
    }
  }
  p_out<-Reduce("+",tmp_Plist)+plot_layout(ncol = length(tmp_Plist),
                                           widths = tmp_width,guides = "collect")
  return(p_out)
}
(p1<-DimPlot(tmp_obj, reduction = "umap",label.size = 5,
             group.by = tmp_Idents,
             # cols = Major_color[levels(tmp_Sobj$Major),]$Color,
             raster = F,label = T, pt.size = 0.1,repel = T)+
    guides(color = guide_legend(ncol = 1,override.aes = list(size = 2.5)))+
    tidydr::theme_dr(xlength = 0.2,ylength = 0.2,
                     arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
    theme(panel.grid = element_blank(),plot.title = element_blank(),
          legend.position = "none",
          axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
          axis.title.y.right = element_text(size = tmp_label,hjust = 0.5,face = "bold")))

(p2<-tmp_obj %>%
    # {.$Major<-factor(.$Major,levels = rev(levels(.$Major)));.} %>%
    myBubblePlot(tmp_obj=.,tmp_Marker=tmp_Marker,tmp_Scale=T))

(pp1<-(as.ggplot(p1)|as.ggplot(p2))+plot_layout(widths = c(1,2.5)))
ggsave2("results/figure/rpca_integrated_res.0.5_marker.pdf",plot = pp1,width = 20,height = 7)
ggsave2("results/figure/rpca_integrated_res.0.5_marker.png",plot = pp1,width = 20,height = 7)

# 11. 细胞类型定义v1-------------------
tmp_obj$CellType<-as.character(tmp_obj$integrated_snn_res.0.5)

tmp_obj@meta.data %<>%
  mutate(CellType=case_when(integrated_snn_res.0.5%in%c("4","5","8","10","19")~"Astro",
                            integrated_snn_res.0.5%in%c("6","7","14","27")~"OPC",
                            integrated_snn_res.0.5%in%c("1","2","11")~"Oli",
                            integrated_snn_res.0.5%in%c("3","18")~"Micro",
                            integrated_snn_res.0.5%in%c("36")~"T cells",
                            integrated_snn_res.0.5%in%c("33")~"Ependymal",
                            integrated_snn_res.0.5%in%c("23","30","31")~"Vasular",
                            integrated_snn_res.0.5%in%c("0","12","13","15","17","20","21","22","24","25","26","28","29","32","34","35","37")~"Neurons",
                            integrated_snn_res.0.5%in%c("9")~"Unknown1",
                            integrated_snn_res.0.5%in%c("16")~"Unknown2",
  ))

(p1<-DimPlot(tmp_obj,group.by = "CellType",label = T))
qsave(tmp_obj,"data/processed_data/Seuobj.CellType.v1.qs")
