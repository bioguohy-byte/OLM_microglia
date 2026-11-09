# 20250722 combine 24 samples------------------------------
setwd(Sys.getenv("MG_PROJECT_ROOT"))
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

# # 1. 加载数据 --------------
# # — 1.1 加载Song的样本 ------------------------
# 
# s_18.5ys <- read.delim2(Sys.getenv("MG_EXTERNAL_037"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_40ys <- read.delim2(Sys.getenv("MG_EXTERNAL_038"),row.names = 1)
# s_50ys <- read.delim2(Sys.getenv("MG_EXTERNAL_039"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_60ys <- read.delim2(Sys.getenv("MG_EXTERNAL_040"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_86ys <- read.delim2(Sys.getenv("MG_EXTERNAL_041"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_88ys <- read.delim2(Sys.getenv("MG_EXTERNAL_042"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_0.6ys <- read.delim2(Sys.getenv("MG_EXTERNAL_043"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_89_1ys <- read.delim2(Sys.getenv("MG_EXTERNAL_044"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_89_2ys <- read.delim2(Sys.getenv("MG_EXTERNAL_045"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_73_2ys <- read.delim2(Sys.getenv("MG_EXTERNAL_046"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_73_1ys <- read.delim2(Sys.getenv("MG_EXTERNAL_047"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_77_1ys <- read.delim2(Sys.getenv("MG_EXTERNAL_048"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_77_2ys <- read.delim2(Sys.getenv("MG_EXTERNAL_049"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_1.2ys <- read.delim2(Sys.getenv("MG_EXTERNAL_050"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_79_1ys <- read.delim2(Sys.getenv("MG_EXTERNAL_051"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_79_2ys <- read.delim2(Sys.getenv("MG_EXTERNAL_052"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_80M_1ys <- read.delim2(Sys.getenv("MG_EXTERNAL_053"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_80M_2ys <- read.delim2(Sys.getenv("MG_EXTERNAL_054"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_80F_1ys <- read.delim2(Sys.getenv("MG_EXTERNAL_055"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_80F_2ys <- read.delim2(Sys.getenv("MG_EXTERNAL_056"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# s_3.2ys <- read.delim2(Sys.getenv("MG_EXTERNAL_057"), stringsAsFactors = F,sep = ' ',header = T,row.names = 1)
# 
# 
# # — 1.2 加载Li的样本数据 ------------------
# l_4D <- Read10X_h5(Sys.getenv("MG_EXTERNAL_058"), use.names = T)
# l_32y <- Read10X_h5(Sys.getenv("MG_EXTERNAL_059"), use.names = T)
# l_48y <- Read10X_h5(Sys.getenv("MG_EXTERNAL_060"), use.names = T)
# l_50y <- Read10X_h5(Sys.getenv("MG_EXTERNAL_061"), use.names = T)
# l_56y <- Read10X_h5(Sys.getenv("MG_EXTERNAL_062"), use.names = T)
# l_60y <- Read10X_h5(Sys.getenv("MG_EXTERNAL_063"), use.names = T)
# l_64f <- Read10X_h5(Sys.getenv("MG_EXTERNAL_064"), use.names = T)
# l_64m <- Read10X_h5(Sys.getenv("MG_EXTERNAL_065"), use.names = T)
# l_68y <- Read10X_h5(Sys.getenv("MG_EXTERNAL_066"), use.names = T)
# 
# save.image(file = 'data/processed_data/00.load_24sample_data.RData')

# 2. 创建seurat对象 ----------------
load('data/processed_data/00.load_24sample_data.RData')

s_18.5ys <- CreateSeuratObject(s_18.5ys, project = "s_18.5ys", min.cells = 3, min.features = 200)
s_40ys <- CreateSeuratObject(s_40ys, project = "s_40ys", min.cells = 3, min.features = 200)
s_50ys <- CreateSeuratObject(s_50ys, project = "s_50ys", min.cells = 3, min.features = 200)
s_60ys <- CreateSeuratObject(s_60ys, project = "s_60ys", min.cells = 3, min.features = 200)
s_86ys <- CreateSeuratObject(s_86ys, project = "s_86ys", min.cells = 3, min.features = 200)
s_88ys <- CreateSeuratObject(s_88ys, project = "s_88ys", min.cells = 3, min.features = 200)
s_0.6ys <- CreateSeuratObject(s_0.6ys, project = "s_0.6ys", min.cells = 3, min.features = 200)
s_89_1ys <- CreateSeuratObject(s_89_1ys, project = "s_89_1ys", min.cells = 3, min.features = 200)
s_89_2ys <- CreateSeuratObject(s_89_2ys, project = "s_89_2ys", min.cells = 3, min.features = 200)
s_73_2ys <- CreateSeuratObject(s_73_2ys, project = "s_73_2ys", min.cells = 3, min.features = 200)
s_73_1ys <- CreateSeuratObject(s_73_1ys, project = "s_73_1ys", min.cells = 3, min.features = 200)
s_77_1ys <- CreateSeuratObject(s_77_1ys, project = "s_77_1ys", min.cells = 3, min.features = 200)
s_77_2ys <- CreateSeuratObject(s_77_2ys, project = "s_77_2ys", min.cells = 3, min.features = 200)
s_1.2ys <- CreateSeuratObject(s_1.2ys, project = "s_1.2ys", min.cells = 3, min.features = 200)
s_79_1ys <- CreateSeuratObject(s_79_1ys, project = "s_79_1ys", min.cells = 3, min.features = 200)
s_79_2ys <- CreateSeuratObject(s_79_2ys, project = "s_79_2ys", min.cells = 3, min.features = 200)
s_80M_1ys <- CreateSeuratObject(s_80M_1ys, project = "s_80M_1ys", min.cells = 3, min.features = 200)
s_80M_2ys <- CreateSeuratObject(s_80M_2ys, project = "s_80M_2ys", min.cells = 3, min.features = 200)
s_80F_1ys <- CreateSeuratObject(s_80F_1ys, project = "s_80F_1ys", min.cells = 3, min.features = 200)
s_80F_2ys <- CreateSeuratObject(s_80F_2ys, project = "s_80F_2ys", min.cells = 3, min.features = 200)
s_3.2ys <- CreateSeuratObject(s_3.2ys, project = "s_3.2ys", min.cells = 3, min.features = 200)


l_4D <- CreateSeuratObject(l_4D, project = "l_4D", min.cells = 3, min.features = 200)
l_48y <- CreateSeuratObject(l_48y, project = "l_48y", min.cells = 3, min.features = 200)
l_32y <- CreateSeuratObject(l_32y, project = "l_32y", min.cells = 3, min.features = 200)
l_50y <- CreateSeuratObject(l_50y, project = "l_50y", min.cells = 3, min.features = 200)
l_56y <- CreateSeuratObject(l_56y, project = "l_56y", min.cells = 3, min.features = 200)
l_60y <- CreateSeuratObject(l_60y, project = "l_60y", min.cells = 3, min.features = 200)
l_64f <- CreateSeuratObject(l_64f, project = "l_64f", min.cells = 3, min.features = 200)
l_64m <- CreateSeuratObject(l_64m, project = "l_64m", min.cells = 3, min.features = 200)
l_68y <- CreateSeuratObject(l_68y, project = "l_68y", min.cells = 3, min.features = 200)

alldata <- merge(x=s_18.5ys,
                 y=c(s_40ys,s_50ys,s_60ys,s_86ys,s_88ys,s_0.6ys,s_89_1ys,s_89_2ys,s_73_2ys,s_73_1ys,s_77_1ys,s_77_2ys,s_79_1ys,s_79_2ys,s_1.2ys,s_80M_1ys,s_80M_2ys,s_80F_1ys,s_80F_2ys,s_3.2ys,
                     l_4D,l_32y,l_48y,l_50y,l_56y,l_60y,l_64f,l_64m,l_68y), 
                 add.cell.ids=c("s_18.5ys","s_40ys","s_50ys","s_60ys","s_86ys",
                                "s_88ys","s_0.6ys","s_89_1ys","s_89_2ys","s_73_2ys",
                                "s_73_1ys","s_77_1ys","s_77_2ys","s_79_1ys","s_79_2ys",
                                "s_1.2ys","s_80M_1ys","s_80M_2ys","s_80F_1ys","s_80F_2ys","s_3.2ys",
                                "l_4D","l_32y","l_48y","l_50y","l_56y",
                                "l_60y","l_64f","l_64m","l_68y"
                 ))

saveRDS(alldata,file = 'data/processed_data/01.24samples_before_filter_alldata.rds')


# 3. 质控QC -------------------------
alldata <- readRDS('data/processed_data/01.24samples_before_filter_alldata.rds')

# — 3.1 计算线粒体、核糖体、血红蛋白百分比------------------
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^MT-")
alldata[["percent.rb"]] <- PercentageFeatureSet(alldata, pattern = "^RP[SL]")
alldata[["percent.HB"]] <- PercentageFeatureSet(alldata, pattern = "^HB")
(p0 <- VlnPlot(alldata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb","percent.HB"), ncol = 5,pt.size = 0,raster=FALSE))
ggsave("results/figure/qc_before.png",p0,width = 30,height = 8)
ggsave("results/figure/qc_before.pdf",p0,width = 30,height = 8)

# — 3.2 针对具体样本调整阈值----------------
del_all <- WhichCells(alldata, expression = (nFeature_RNA > 4000 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_18.5ys")
data.filt <- subset(alldata, cells = setdiff(WhichCells(alldata), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 3000 | nFeature_RNA < 200 | percent.mt > 20) & orig.ident == "s_40ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 5500 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_50ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 5500 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_60ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 5500 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_86ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 5000 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_88ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 4000 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_0.6ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 4500 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_89_1ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 4500 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_89_2ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 2800 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_73_2ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 3000 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_73_1ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 2500 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_77_1ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 2500 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_77_2ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 5000 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_1.2ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 5000 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_79_1ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 5000 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_79_2ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 2500 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_80M_1ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 2500 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_80M_2ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 5000 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_80F_1ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 5000 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_80F_2ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))
del_all <- WhichCells(data.filt, expression = (nFeature_RNA > 5500 | nFeature_RNA < 200 | percent.mt > 5) & orig.ident == "s_3.2ys")
data.filt <- subset(data.filt, cells = setdiff(WhichCells(data.filt), del_all))



del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.mt>5) & orig.ident == "l_4D")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.mt>5) & orig.ident == "l_48y")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>8300 | nFeature_RNA<200 | percent.mt>10) & orig.ident == "l_50y")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>8300 | nFeature_RNA<200 | percent.mt>10) & orig.ident == "l_56y")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.mt>10) & orig.ident == "l_60y")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.HB>1 | percent.mt>15) & orig.ident == "l_64f")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.HB>1 | percent.mt>20) & orig.ident == "l_64m")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>8600 | nFeature_RNA<200 | percent.mt>10) & orig.ident == "l_68y")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))
del_all <- WhichCells(data.filt,expression = (nFeature_RNA>7500 | nFeature_RNA<200 | percent.mt>20) & orig.ident == "l_32y")
data.filt <- subset(data.filt,cells=setdiff(WhichCells(data.filt),del_all))

data.filt$sample <- data.filt$orig.ident
saveRDS(data.filt,file = 'data/processed_data/02.24samples_after_filter_data.filt.rds')

# — 3.3 质控前后作图对比--------------------
p1 <- VlnPlot(alldata, features = "nFeature_RNA", pt.size = 0) + NoLegend() + theme()
p2 <- VlnPlot(alldata, features = "nCount_RNA", pt.size = 0) + NoLegend()+ theme()
p3 <- VlnPlot(alldata, features = "percent.mt", pt.size = 0) + NoLegend()+ theme()
p4 <- VlnPlot(alldata, features = "percent.rb", pt.size = 0) + NoLegend()+ theme()
p5 <- VlnPlot(alldata, features = "percent.HB", pt.size = 0)+ NoLegend()+ theme()
p6 <- VlnPlot(data.filt, features = "nFeature_RNA", pt.size = 0) + NoLegend() + theme()
p7 <- VlnPlot(data.filt, features = "nCount_RNA", pt.size = 0) + NoLegend()+ theme()
p8 <- VlnPlot(data.filt, features = "percent.mt", pt.size = 0) + NoLegend()+ theme()
p9 <- VlnPlot(data.filt, features = "percent.rb", pt.size = 0) + NoLegend()+ theme()
p10 <- VlnPlot(data.filt, features = "percent.HB", pt.size = 0)+ NoLegend()+ theme()
p<-(p1|p2|p3|p4|p5)/(p6|p7|p8|p9|p10)
ggsave("results/figure/qc_round1_summary.png",plot = p,width = 40,height = 9)
ggsave("results/figure/qc_round1_summary.pdf",plot = p,width = 40,height = 9)