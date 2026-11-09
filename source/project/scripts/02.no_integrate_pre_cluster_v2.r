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
# 1.分大群 ---------------
# — 1.1 ----------------------
seu.obj <- readRDS('data/processed_data/04.no_integrate_after_cluster.rds')

(tmp_Marker<-data.frame(Gname=c('AQP4','ALDH1L1','GFAP', #AST
                                'CLIC6', #Choroid.plexus
                                'DCX','STMN1','STMN2','SOX4', #NB
                                'SYT1','SV2B', #GC
                                'RELN','CCK','GAD2', #IN
                                'OLIG1','PDGFRA', #OPC
                                'MBP','PLP1','MAL', #OLG
                                'CSF1R','P2RY12','CTSS', #MG
                                'VWF', #Endo
                                'DNAH9' #Epen
),
Type=c(rep("AST",3),
       rep("Choroid.plexus",1),
       rep("NB",4),
       rep("GC",2),
       rep("IN",3),
       rep("OPC",2),
       rep("OLG",3),
       rep("MG",3),
       rep("Endo",1),
       rep("Epen",1)
)))


tmp_obj <- seu.obj

unique_clusters <- unique(tmp_obj$RNA_snn_res.0.4)
length(unique_clusters)

tmp_obj$RNA_snn_res.0.4<-factor(tmp_obj$RNA_snn_res.0.4,levels = seq(0,length(unique_clusters)-1))
(p1<-DimPlot(tmp_obj,group.by = "RNA_snn_res.0.4",label = T,raster = F)+NoLegend())
tmp_Idents<-"RNA_snn_res.0.4"

myBubblePlot <- function(tmp_obj,tmp_Marker,tmp_Scale,tmp_label) {
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
    myBubblePlot(tmp_obj=.,tmp_Marker=tmp_Marker,tmp_Scale=T,tmp_label=9))

(pp1<-(as.ggplot(p1)|as.ggplot(p2))+plot_layout(widths = c(1,2.5)))

ggsave('results/figure/no_integrate_define_celltype_v1.png',plot = pp1,width = 20,height = 7)
ggsave('results/figure/no_integrate_define_celltype_v1.pdf',plot = pp1,width = 20,height = 7)

# 2.粗定义细胞群 -----------------------------
print("-------------------load data------------------")
seu.obj <- readRDS('data/processed_data/04.no_integrate_after_cluster.rds')
tmp_obj <- seu.obj


tmp_obj$Raw_type<-as.character(tmp_obj$RNA_snn_res.0.4)

tmp_obj@meta.data %<>% 
  mutate(Raw_type=case_when(RNA_snn_res.0.4%in%c("4","5","10","20")~"AST",
                            RNA_snn_res.0.4%in%c("34")~"Choroid.plexus",
                            RNA_snn_res.0.4%in%c("2","3","6","7","12","13","14","15","17","18","21","22","23","25","27","28","30","31","33","35")~"Neuron",
                            RNA_snn_res.0.4%in%c("9","16","29")~"OPC",
                            RNA_snn_res.0.4%in%c("0","1")~"OLG",
                            RNA_snn_res.0.4%in%c("8","11","26")~"MG",
                            RNA_snn_res.0.4%in%c("19")~"Endo",
                            RNA_snn_res.0.4%in%c("24","32")~"Epen"))

(p3<-DimPlot(tmp_obj,group.by = "Raw_type",label = T,raster = F))
ggsave('results/figure/no_integrate_pre_celltype_v2.png',plot = p3,width = 10,height = 7)
ggsave('results/figure/no_integrate_pre_celltype_v2.pdf',plot = p3,width = 10,height = 7)

saveRDS(tmp_obj,file = 'data/processed_data/05.no_integrate_pre_cluster_v2.rds')


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
seu.obj <- readRDS('data/processed_data/05.no_integrate_pre_cluster_v2.rds')
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
save(all_DL,summary,tmp_Plist,file = paste("data/processed_data/06.no_integrate_DoubletFinder.result_v2.RData",sep = ""))


# --------------------------------
load("data/processed_data/06.no_integrate_DoubletFinder.result_v2.RData")

summary

Reduce("+",tmp_Plist)+plot_layout(ncol = 4,guides = "collect")


tmp_obj@meta.data %<>% {.$DoubletFinder<-"Singlet";.} %<>% 
  {.[rownames(all_DL[all_DL$if_DL%in%"Doublet",]),]$DoubletFinder<-"Doublet";.}

saveRDS(tmp_obj,file = 'data/processed_data/07.no_integrate_after_df_v2.rds')

p1 <- DimPlot(tmp_obj, reduction = "umap", group.by = "DoubletFinder",pt.size = 0.2,raster = F)
p2 <- DimPlot(tmp_obj, reduction = "umap", group.by = "RNA_snn_res.0.4",pt.size = 0.2,raster = F,label = T,repel = T)
p3 <- DimPlot(tmp_obj, reduction = "umap", group.by = "Raw_type",pt.size = 0.2,raster = F,label = T)

p_df <- p1|p2|p3
ggsave("results/figure/no_integrate_after_df_v2.png",plot = p_df,width = 25,height = 8)
ggsave("results/figure/no_integrate_after_df_v2.pdf",plot = p_df,width = 25,height = 8)


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

# 1. load data------------------------
tmp_obj <- readRDS('data/processed_data/07.no_integrate_after_df_v2.rds')


alldata<-tmp_obj

alldata[["percent.mito"]] <- PercentageFeatureSet(alldata, pattern = "^MT-")
alldata[["percent.HB"]] <- PercentageFeatureSet(alldata, pattern = "^HB")

tmp_Plist<-list()
for (i in c("nFeature_RNA","nCount_RNA","percent.mito","percent.HB")) {
  p<-VlnPlot(alldata, features = i, group.by = "orig.ident", pt.size = 0) + 
    NoLegend() + 
    theme(axis.title = element_blank())
  tmp_Plist[[i]]<-p
}

(p_before<-Reduce("+",tmp_Plist)+
    plot_layout(ncol = length(tmp_Plist))+
    plot_annotation(title = "Before filtered",theme = theme(plot.title = element_text(size = 20))))


tmp_Summary <- table(alldata$orig.ident) %>% as.data.frame() %>% {.$Step="Raw";.}

# step1
data.filt <- alldata %>% subset(.,subset=DoubletFinder%in%"Singlet")
tmp_Summary<-table(data.filt$orig.ident) %>% as.data.frame() %>% {.$Step="Singlets";.} %>% 
  rbind(tmp_Summary,.)

# step2
tmp_minGenes<-200
data.filt <- data.filt %>% subset(.,cells=setdiff(WhichCells(.),
                                                  WhichCells(.,expression = (nFeature_RNA<=tmp_minGenes))))
tmp_Summary<-table(data.filt$orig.ident) %>% as.data.frame() %>% {.$Step=paste("nGenes > ",tmp_minGenes,sep = "");.} %>% 
  rbind(tmp_Summary,.)

# step3
tmp_maxGenes<-10000
data.filt <- data.filt %>% subset(.,cells=setdiff(WhichCells(.),
                                                  WhichCells(.,expression = (nFeature_RNA>=tmp_maxGenes))))
tmp_Summary<-table(data.filt$orig.ident) %>% as.data.frame() %>% {.$Step=paste("nGenes < ",tmp_maxGenes,sep = "");.} %>% 
  rbind(tmp_Summary,.)

# step4
tmp_MT<-20
data.filt <- data.filt %>% subset(.,cells=setdiff(WhichCells(.),
                                                  WhichCells(.,expression = (percent.mito>=tmp_MT))))
tmp_Summary<-table(data.filt$orig.ident) %>% as.data.frame() %>% {.$Step=paste("MT% < ",tmp_MT,sep = "");.} %>% 
  rbind(tmp_Summary,.)

# step5
tmp_HB<-1
data.filt <- data.filt %>% subset(.,cells=setdiff(WhichCells(.),
                                                  WhichCells(.,expression = (percent.HB>=tmp_HB))))
tmp_Summary<-table(data.filt$orig.ident) %>% as.data.frame() %>% {.$Step=paste("HB% < ",tmp_HB,sep = "");.} %>% 
  rbind(tmp_Summary,.)

# saveRDS(data.filt,file = 'data/processed_data/0507test.rds')
# data.filt <- readRDS('data/processed_data/0507test.rds')

# Before
i<-"nFeature_RNA"
p1<-VlnPlot(alldata, features = i, group.by = "orig.ident", pt.size = 0) + 
  geom_hline(yintercept = tmp_minGenes,lty=4,lwd=0.6,alpha=0.3)+
  geom_hline(yintercept = tmp_maxGenes,lty=4,lwd=0.6,alpha=0.3)+
  labs(title = "nGenes")+
  NoLegend()+ 
  theme(axis.title = element_blank())

i<-"nCount_RNA"
p2<-VlnPlot(alldata, features = i, group.by = "orig.ident", pt.size = 0) + 
  labs(title = "nUMIs")+
  NoLegend() + 
  theme(axis.title = element_blank())

i<-"percent.mito"
p3<-VlnPlot(alldata, features = i, group.by = "orig.ident", pt.size = 0) + 
  geom_hline(yintercept = tmp_MT,lty=4,lwd=0.6,alpha=0.3)+
  labs(title = "Mitochondrial ratio (%)")+
  NoLegend() + 
  theme(axis.title = element_blank())

i<-"percent.HB"
p4<-VlnPlot(alldata, features = i, group.by = "orig.ident", pt.size = 0) + 
  geom_hline(yintercept = tmp_HB,lty=4,lwd=0.6,alpha=0.3)+
  labs(title = "Erythrocytic ratio (%)")+
  NoLegend() + 
  theme(axis.title = element_blank())

(p_before<-(p1+p2+p3+p4)+
    plot_layout(ncol = 4)+
    plot_annotation(title = "Before filtered",theme = theme(plot.title = element_text(size = 20))))

# After
i<-"nFeature_RNA"
p1<-VlnPlot(data.filt, features = i, group.by = "orig.ident", pt.size = 0) + 
  geom_hline(yintercept = tmp_minGenes,lty=4,lwd=0.6,alpha=0.3)+
  geom_hline(yintercept = tmp_maxGenes,lty=4,lwd=0.6,alpha=0.3)+
  labs(title = "nGenes")+
  NoLegend()+ 
  theme(axis.title = element_blank())

i<-"nCount_RNA"
p2<-VlnPlot(data.filt, features = i, group.by = "orig.ident", pt.size = 0) + 
  labs(title = "nUMIs")+
  NoLegend() + 
  theme(axis.title = element_blank())

i<-"percent.mito"
p3<-VlnPlot(data.filt, features = i, group.by = "orig.ident", pt.size = 0) + 
  geom_hline(yintercept = tmp_MT,lty=4,lwd=0.6,alpha=0.3)+
  labs(title = "Mitochondrial ratio (%)")+
  NoLegend() + 
  theme(axis.title = element_blank())

i<-"percent.HB"
p4<-VlnPlot(data.filt, features = i, group.by = "orig.ident", pt.size = 0) + 
  geom_hline(yintercept = tmp_HB,lty=4,lwd=0.6,alpha=0.3)+
  labs(title = "Erythrocytic ratio (%)")+
  NoLegend() + 
  theme(axis.title = element_blank())

(p_after<-(p1+p2+p3+p4)+
    plot_layout(ncol = 4)+
    plot_annotation(title = "After filtered",theme = theme(plot.title = element_text(size = 20))))

(p<-as.ggplot(p_before)/as.ggplot(p_after))
ggsave("results/figure/qc_round2_summary.Violin.Filtered_v2.pdf",p,width = 25,height = 7)


tmp_Summary<-tmp_Summary[tmp_Summary$Step%in%"HB% < 1",] %>% {.$Step<-"nCells after filtered";.} %>% rbind(tmp_Summary,.)
tmp_Summary$Step<-factor(tmp_Summary$Step,levels = rev(unique(tmp_Summary$Step)))
tmp_Plist<-list()
for (i in 1:length(levels(tmp_Summary$Var1))) {
  if (i==1) {
    (pp<-ggplot(tmp_Summary[tmp_Summary$Var1%in%levels(tmp_Summary$Var1)[i],],aes(x = Freq,y=Step))+
       geom_bar(stat= 'identity',color="grey70",fill="grey70")+
       geom_text(aes(x=Freq,label=Freq),size=4,hjust = 1,color="black")+
       # scale_x_continuous(expand = c(0, 0.1))+
       # xlim(0,16)+
       facet_grid(~Var1)+
       scale_color_manual(values = c("#F9423AFF","#FFA500","#8DC6FF","#00A087B2"))+
       scale_fill_manual(values = c("#F9423AFF","#FFA500","#8DC6FF","#00A087B2"))+
       # scale_color_gradientn(colors = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
       # labs(title = "Protein family")+
       theme_bw()+
       theme(axis.title = element_blank(),
             axis.text.x = element_text(size=12,colour = "black",angle = 45,hjust = 1),
             axis.text.y = element_text(size=12,colour = "black"),
             # axis.ticks.y = element_blank(),axis.line.y = element_blank(),
             # legend.position = "none",
             legend.text= element_text(size=12),
             legend.title = element_text(size=12),
             strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
             panel.grid = element_blank(),
             strip.placement = "outside",
             strip.switch.pad.wrap = unit(0.3, "mm"),
             plot.title = element_blank(),
             strip.text = element_text(size =13,face = "bold",colour = "white"))) 
    tmp_Plist[[as.character(i)]]<-pp
  } else {
    (pp<-ggplot(tmp_Summary[tmp_Summary$Var1%in%levels(tmp_Summary$Var1)[i],],aes(x = Freq,y=Step))+
       geom_bar(stat= 'identity',color="grey70",fill="grey70")+
       geom_text(aes(x=Freq,label=Freq),size=4,hjust = 1,color="black")+
       # scale_x_continuous(expand = c(0, 0.1))+
       # xlim(0,16)+
       facet_grid(~Var1)+
       scale_color_manual(values = c("#F9423AFF","#FFA500","#8DC6FF","#00A087B2"))+
       scale_fill_manual(values = c("#F9423AFF","#FFA500","#8DC6FF","#00A087B2"))+
       # scale_color_gradientn(colors = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
       # labs(title = "Protein family")+
       theme_bw()+
       theme(axis.title = element_blank(),
             axis.text.x = element_text(size=12,colour = "black",angle = 45,hjust = 1),
             axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
             # legend.position = "none",
             legend.text= element_text(size=12),
             legend.title = element_text(size=12),
             strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
             panel.grid = element_blank(),
             strip.placement = "outside",
             strip.switch.pad.wrap = unit(0.3, "mm"),
             plot.title = element_blank(),
             strip.text = element_text(size =13,face = "bold",colour = "white"))) 
    tmp_Plist[[as.character(i)]]<-pp
  }
}

(p5<-Reduce("+",tmp_Plist)+plot_layout(ncol = length(tmp_Plist)))
(ppp<-p_before/p_after/p5)
ggsave("results/figure/qc_round2_summary.pdf",ppp,width = 30,height = 9)

save(alldata,data.filt,tmp_Summary,file = "data/saved_data/qc.summary_v2.RData")
saveRDS(data.filt,file = "data/processed_data/08.no_integrate_after_qc_round2_v2.rds")

# cca方法整合 ------------------------
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

library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)

library(future) # 新增
plan("multicore", workers = 5) # 根据你的服务器负载调整 workers 数量
options(future.globals.maxSize = 80 * 1024^3) # 设置 80GB 内存上限


data.filt <- readRDS("data/processed_data/08.no_integrate_after_qc_round2_v2.rds")

data.filt.list <- SplitObject(data.filt, split.by = "sample")
data.filt.list <- lapply(X = data.filt.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = data.filt.list)

# perform integration
data.filt.anchors <- FindIntegrationAnchors(object.list = data.filt.list, anchor.features = features)
data.filt.combined <- IntegrateData(anchorset = data.filt.anchors)
saveRDS(data.filt.combined,file = 'data/processed_data/09.integrate_cca_before_cluster_v2.rds')

DefaultAssay(data.filt.combined) <- "integrated"
data.filt.combined <- ScaleData(data.filt.combined)
data.filt.combined <- RunPCA(data.filt.combined, npcs = 70, verbose = FALSE)
data.filt.combined <- RunUMAP(data.filt.combined, reduction = "pca", dims = 1:70)
data.filt.combined <- FindNeighbors(data.filt.combined, reduction = "pca", dims = 1:70)
data.filt.combined <- FindClusters(data.filt.combined, resolution = seq(0.1,2,0.1))

saveRDS(data.filt.combined,file = 'data/processed_data/10.integrate_cca_after_cluster_v2.rds')

p1 <- DimPlot(data.filt.combined,label = T)
ggsave(Sys.getenv("MG_EXTERNAL_067"),p1,width = 6,height = 8)



