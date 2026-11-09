# — Celldefine : Major v1 ---------------------------------------------------------------
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

tmp_obj <- readRDS("data/processed_data/10.integrate_cca_after_cluster.rds")

# 重新排列 levels，将其强制转换为数值排序后再转回字符
levels(tmp_obj$integrated_snn_res.0.4) <- sort(as.numeric(levels(tmp_obj$integrated_snn_res.0.4)))


tmp_label=13
(p1<-DimPlot(tmp_obj, reduction = "umap",label.size = 5,
             group.by = "integrated_snn_res.0.4",
             # cols = Major_color[levels(tmp_Sobj$Major),]$Color,
             raster = F,label = T, pt.size = 0.1,repel = T)+
    guides(color = guide_legend(ncol = 1,override.aes = list(size = 2.5))) +
    tidydr::theme_dr(xlength = 0.2,ylength = 0.2,
                     arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
    theme(panel.grid = element_blank(),plot.title = element_blank(),
          legend.position = "none",
          axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
          axis.title.y.right = element_text(size = tmp_label,hjust = 0.5,face = "bold")))


# 扩充后的各群 marker 基因定义（包含兴奋性神经元 Exc）
tmp_Marker <- data.frame(
  Gname = c(
    # AST (星状胶质细胞)
    'AQP4', 'ALDH1L1', 'GFAP', 'S100B', 'SLC1A2', 
    # Choroid.plexus (脉络丛)
    'CLIC6', 'TTR', 'FOLR1', 'ENPP2', 'KRT18', 
    # Exc (兴奋性神经元)
    'SLC17A7', 'CAMK2A', 'GRIN1', 'SATB2', 'NEUROD6',
    # IN (抑制性神经元)
    'GAD1', 'GAD2', 'RELN', 'CCK', 'SST', 
    # OPC (少突前体细胞)
    'OLIG1', 'PDGFRA', 'CSPG4', 'SOX10', 'GPR17', 
    # OLG (成熟少突胶质细胞)
    'MBP', 'PLP1', 'MAL', 'MOG', 'MOBP', 
    # MG (小胶质细胞)
    'CSF1R', 'P2RY12', 'CTSS', 'CX3CR1', 'C1QA', 
    # Endo (内皮细胞)
    'VWF', 'CLDN5', 'PECAM1', 'ABCB1', 'TIE1', 
    # Epen (室管膜细胞)
    'DNAH9', 'FOXJ1', 'CAPS', 'PIFO', 'SNTN'
  ),
  Type = c(
    rep("AST", 5),
    rep("Choroid.plexus", 5),
    rep("Exc", 5), # 新增的兴奋性神经元分类
    rep("IN", 5),
    rep("OPC", 5),
    rep("OLG", 5),
    rep("MG", 5),
    rep("Endo", 5),
    rep("Epen", 5)
  )
)
DefaultAssay(tmp_obj) <- "RNA"
tmp_Idents<-"integrated_snn_res.0.4"

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

ggsave('results/figure/integrate_cca_define_celltype_v1.png',plot = pp1,width = 20,height = 9)
ggsave('results/figure/integrate_cca_define_celltype_v1.pdf',plot = pp1,width = 20,height = 9)


# define--------------------------------
# 2.粗定义细胞群 -----------------------------
print("-------------------load data------------------")
tmp_obj$CellType<-as.character(tmp_obj$integrated_snn_res.0.4)

tmp_obj@meta.data %<>% 
  mutate(CellType=case_when(integrated_snn_res.0.4%in%c("1","2")~"AST",
                            integrated_snn_res.0.4%in%c("16")~"Choroid.plexus",
                            integrated_snn_res.0.4%in%c("3","6","8","9","17","18","20","21","24","25","26")~"Exc",
                            integrated_snn_res.0.4%in%c("4","5","10","15")~"IN",
                            integrated_snn_res.0.4%in%c("12")~"OPC",
                            integrated_snn_res.0.4%in%c("0","19","22")~"OLG",
                            integrated_snn_res.0.4%in%c("13","23")~"MG",
                            integrated_snn_res.0.4%in%c("11","14")~"Endo",
                            integrated_snn_res.0.4%in%c("7")~"Epen"))

(p3<-DimPlot(tmp_obj,group.by = "CellType",label = T,raster = F))
ggsave('results/figure/integrate_cca_CellType_v1.png',plot = p3,width = 10,height = 7)
ggsave('results/figure/integrate_cca_CellType_v1.pdf',plot = p3,width = 10,height = 7)

saveRDS(tmp_obj,file = 'data/processed_data/11.integrate_cca_CellType_v1.rds')


library(qs)
qsave(tmp_obj, 'data/processed_data/11.integrate_cca_CellType_v1.qs', nthreads = 4)
tmp_obj <- qread('data/processed_data/11.integrate_cca_CellType_v1.qs', nthreads = 4)
