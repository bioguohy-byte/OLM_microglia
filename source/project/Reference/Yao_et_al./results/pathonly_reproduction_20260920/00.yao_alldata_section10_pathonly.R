# Path-only executable copy of original Yao script section 10 (original lines 541-661).
# Plotting logic, markers, function body, and parameters are unchanged; only the
# migrated working/input/output paths and minimal package-loading setup are supplied.
.libPaths(c(file.path(Sys.getenv("MG_WORKSPACE_ROOT"), "R_libs"), .libPaths()))
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Yao_et_al."))
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(magrittr)
library(ggplotify)
library(qs)

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
ggsave2("results/pathonly_reproduction_20260920/rpca_integrated_res.0.5_marker.pdf",plot = pp1,width = 20,height = 7)
ggsave2("results/pathonly_reproduction_20260920/rpca_integrated_res.0.5_marker.png",plot = pp1,width = 20,height = 7)
