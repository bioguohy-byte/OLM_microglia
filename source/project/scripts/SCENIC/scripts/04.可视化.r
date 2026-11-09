# 4.1 SCENIC依赖包安装---------------
## SCENIC需要一些依赖包，先安装好
# BiocManager::install(c("AUCell", "RcisTarget"))
# BiocManager::install(c("GENIE3"))
# BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# BiocManager::install(c("DT", "NMF", "pheatmap", "R2HTML", "Rtsne"))
# BiocManager::install(c("doMC", "doRNG"))
# devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
# devtools::install_github("aertslab/SCENIC")

library(SCENIC)
packageVersion("SCENIC")
#[1] ‘1.2.4’


# 4.2 提取out_SCENIC.loom文件的信息-------------
##可视化
rm(list=ls())
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)

# 多线程运行
library(future)
plan(multisession, workers = 30)
nbrOfWorkers()
options(future.globals.maxSize = 1999999 * 1924^2)
# library(scRNAseq)

## 提取 out_SCENIC.loom 信息
#inputDir='./outputs/'
#scenicLoomPath=file.path(inputDir,'out_SCENIC.loom')
library(SCENIC)
loom <- open_loom(Sys.getenv("MG_EXTERNAL_069")) 
# 24samples subset出来的小胶质
load(Sys.getenv("MG_EXTERNAL_070"))

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)

# 4.3 得到可视化的文件--------------
sub_regulonAUC <- regulonAUC[,match(colnames(MG.combined),colnames(regulonAUC))]
dim(sub_regulonAUC)
MG.combined 
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(MG.combined))
#[1] TRUE

cellClusters <- data.frame(row.names = colnames(MG.combined), 
                           sub_clusters = as.character(MG.combined$sub_clusters2))

### 这里我把带有cellTypes的全都注释掉了，因为在我们的MG.combined中并没有sub_celltype这一meta.data，我们只有已经分好的sub_clusters
### 所以只用这个sub_clusters就可以了其实
# cellTypes <- data.frame(row.names = colnames(MG.combined), 
#                         celltype = MG.combined$sub_celltype)
# head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 
# save(sub_regulonAUC,cellTypes,cellClusters,MG.combined,
#      file = 'for_rss_and_visual.Rdata')

save(sub_regulonAUC,cellClusters,MG.combined,
     file = Sys.getenv("MG_EXTERNAL_071"))

# 4.4 TF活性均值-----------------------
load(Sys.getenv("MG_EXTERNAL_071"))
sub_regulonAUC[1:4,1:2] 
dim(sub_regulonAUC)
# [1]   339 12578

# Split the cells by cluster:
selectedResolution <- "sub_clusters" # select resolution
cellsPerGroup <- split(rownames(cellClusters), 
                       cellClusters[,selectedResolution]) 

# 去除extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)
#[1]   339 12578（339个转录因子 12578个细胞）

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

# Scale expression. 
# Scale函数是对列进行归一化，所以要把regulonActivity_byGroup转置成细胞为行，基因为列
# 参考：https://www.jianshu.com/p/115d07af3029
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
# 同一个regulon在不同cluster的scale处理
dim(regulonActivity_byGroup_Scaled)
#[1] 339 6
regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

# 热图查看TF分布
pheatmap(regulonActivity_byGroup_Scaled)

# 4.5 rss查看特异的TF转录因子----------------
library(openxlsx)
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellClusters[colnames(sub_regulonAUC), 
                                           selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss)
p1 <- rssPlot$plot

# 保存图片
ggsave(Sys.getenv("MG_EXTERNAL_072"),p1,width = 7,height = 12)
ggsave(Sys.getenv("MG_EXTERNAL_073"),p1,width = 7,height = 12)


library(dplyr) 
library(ggplot2)
rss=regulonActivity_byGroup_Scaled
head(rss)
library(dplyr) 
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster =   colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(cluster = top5$cluster) 
n = rss[top5$path,] 
rownames(rowcn) = rownames(n)
write.xlsx(df,Sys.getenv("MG_EXTERNAL_074"))

# 绘图
# 为这些注释分配颜色
annotation_colors <- list(cluster = c(M0 = "#D43F3AFF", 
                                      M1 = "#EEA236FF", 
                                      M2 = "#5CB85CFF",
                                      M3 = "#46B8DAFF", 
                                      M4 = "#357EBDFF", 
                                      M5 = "#9632B8FF" ))

(p2 <- pheatmap(n,
                colorRampPalette(brewer.pal(8, "Blues"))(100),
                annotation_row = rowcn,
                annotation_colors = annotation_colors,
                show_rownames = T,
                cluster_cols = F,
                cluster_rows = F,
                main = "MG Clusters SCENIC"))
# 加了这个cluster_cols=F之后，横坐标的顺序才算正常了
p2

# 保存图片
ggsave(Sys.getenv("MG_EXTERNAL_075"),p2,width = 8,height = 10)
ggsave(Sys.getenv("MG_EXTERNAL_076"),p2,width = 8,height = 10)


# AUCell score 绘图------------------------
# 颜色设定
pal <- viridis(n = 10, option = "D" )
tmp.regulon <- rownames(rowcn)
tmp.SeuObj <- MG.combined
# 
# feature <- gsub("\\(\\+\\)", "", tmp.regulon)
# 这里我暂时只用单个基因和regulon，才能跑通
feature <- 'NKX6-2'
tmp.regulon <- 'NKX6-2(+)'

plot.data <- FeaturePlot(tmp.SeuObj, features = feature, order = T, cols = pal, raster = F)$data
# plot.data$ZNF573 %>% max()
fdata <- FetchData(object = tmp.SeuObj, vars = feature)
# fdata$ZNF573 %>% max()
plot.data$feature <- fdata[match(rownames(plot.data), rownames(fdata)), ]

plot1 <- ggplot(plot.data, aes(UMAP_1, UMAP_2, color = feature)) + 
  geom_point(size = 0.001) +
  scale_color_gradientn(colors = pal) +
  theme_classic() +
  theme(
    # 删除坐标轴文本、刻度、标题和图例
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) + 
  NoAxes() +
  ggtitle(paste0(gsub("\\(\\+\\)", "", tmp.regulon), " gene"," expression ")) + 
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) # 居中标题

plot1

### auc of tfs
tmp.df <- data.frame(tmp.regulon = sub_regulonAUC@assays@data@listData[["AUC"]][tmp.regulon, ]) #这里的regulonAUC我改成sub_regulonAUC
# identical(Cells(tmp.SeuObj), rownames(tmp.df))
tmp.SeuObj@meta.data$tmp.regulon <- tmp.df$tmp.regulon
plot.data <- FeaturePlot(tmp.SeuObj, features = "tmp.regulon", order = T, cols = pal, raster = F)$data
# plot.data$tmp.regulon %>% max()
plot.data$tmp.regulon <- tmp.df[match(rownames(plot.data), rownames(tmp.df)), ]
# plot.data$tmp.regulon %>% max()

plot2 <- ggplot(plot.data, aes(UMAP_1, UMAP_2, color = tmp.regulon)) +
  geom_point(size = 0.001) +
  scale_color_gradientn(colors = pal) +
  theme_classic() +
  theme(
    # 删除坐标轴文本、刻度、标题和图例
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) + 
  NoAxes() +
  ggtitle(paste0("AUCell score of the ", gsub("\\(\\+\\)", "", tmp.regulon), " regulon")) + 
  theme(plot.title = element_text(hjust = 0.5,face = "bold"))
plot <- plot1 + plot2

plot
ggsave(paste0(fwd, "/", gsub("\\(\\+\\)", "", tmp.regulon), ".png"), plot, width = 8, height = 4)


