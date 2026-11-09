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
                            RNA_snn_res.0.4%in%c("14","17")~"NB",
                            RNA_snn_res.0.4%in%c("2","3","6","7","13","15","22","23","25","27","28","31","33","35")~"GC",
                            RNA_snn_res.0.4%in%c("12","18","21","30")~"IN",
                            RNA_snn_res.0.4%in%c("9","16","29")~"OPC",
                            RNA_snn_res.0.4%in%c("0","1")~"OLG",
                            RNA_snn_res.0.4%in%c("8","11","26")~"MG",
                            RNA_snn_res.0.4%in%c("19")~"Endo",
                            RNA_snn_res.0.4%in%c("24","32")~"Epen"))

(p3<-DimPlot(tmp_obj,group.by = "Raw_type",label = T,raster = F))
ggsave('results/figure/no_integrate_pre_celltype_v1.png',plot = p3,width = 10,height = 7)
ggsave('results/figure/no_integrate_pre_celltype_v1.pdf',plot = p3,width = 10,height = 7)

saveRDS(tmp_obj,file = 'data/processed_data/05.no_integrate_pre_cluster.rds')
