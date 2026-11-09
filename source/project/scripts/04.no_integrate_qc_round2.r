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
tmp_obj <- readRDS('data/processed_data/07.no_integrate_after_df.rds')


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
ggsave("results/figure/qc_round2_summary.Violin.Filtered.pdf",p,width = 25,height = 7)


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

save(alldata,data.filt,tmp_Summary,file = "data/saved_data/qc.summary.RData")
saveRDS(data.filt,file = "data/processed_data/08.no_integrate_after_qc_round2.rds")
