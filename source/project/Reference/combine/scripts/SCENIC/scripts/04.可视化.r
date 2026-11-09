annotation_colors <- list(
  cluster = c("MG0" = "#4DAF4A", "MG1" = "#A65628", "MG2" = "#377EB8", 
              "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG5" = "#8B0000")
)
# 4.2 提取out_SCENIC.loom文件的信息-------------
##可视化
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
library(qs)

# 加载数据------------
loom <- open_loom(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/SCENIC/results/out_SCENIC.loom")) 
tmp_obj <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))

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
sub_regulonAUC <- regulonAUC[,match(colnames(tmp_obj),colnames(regulonAUC))]
dim(sub_regulonAUC)
tmp_obj 
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(tmp_obj))
#[1] TRUE

cellClusters <- data.frame(row.names = colnames(tmp_obj), 
                           sub_clusters = as.character(tmp_obj$CellType1))

### 这里我把带有cellTypes的全都注释掉了，因为在我们的tmp_obj中并没有sub_celltype这一meta.data，我们只有已经分好的sub_clusters
### 所以只用这个sub_clusters就可以了其实
# cellTypes <- data.frame(row.names = colnames(tmp_obj), 
#                         celltype = tmp_obj$sub_celltype)
# head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 
# save(sub_regulonAUC,cellTypes,cellClusters,tmp_obj,
#      file = 'for_rss_and_visual.Rdata')

qsavem(sub_regulonAUC, cellClusters, tmp_obj, file = file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/SCENIC/results/for_rss_and_visual.qsm"))

# 可视化2---------------------
# Fig.1 SCENIC--------------
# TF活性均值
qload(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/SCENIC/results/for_rss_and_visual.qsm"))
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

#  rss查看特异的TF转录因子
library(openxlsx)
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellClusters[colnames(sub_regulonAUC), 
                                           selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss)
p1 <- rssPlot$plot











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

top5 <- top5[!duplicated(top5$path), ]
rowcn = data.frame(cluster = top5$cluster) 
n = rss[top5$path,] 
rownames(rowcn) = top5$path  # 去重后直接赋值即可
rownames(n) = top5$path

write.xlsx(df,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/SCENIC/results/SCENIC.xlsx"),overwrite = T)

# Define paths
base_path <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/SCENIC/results/figure/")
pdf_file <- paste0(base_path, "SCENIC_heatmap.pdf")
png_file <- paste0(base_path, "SCENIC_heatmap.png")

# Generate and save directly via pheatmap (Example for PDF)
p2 <- pheatmap(n,
               color = colorRampPalette(rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", 
                                              "#ffe6b7", "#aadce0", "#72bcd5", "#528fad", 
                                              "#376795", "#1e466e")))(100),
               annotation_row = rowcn,
               annotation_colors = annotation_colors,
               show_rownames = TRUE,
               cluster_cols = FALSE,
               cluster_rows = FALSE,
               main = "MG Clusters SCENIC",
               filename = pdf_file, # Adding this directly saves the plot
               width = 6, 
               height = 12)

# 1. Save as PDF
pdf(file = pdf_file, width = 5, height = 8, bg = "white")
pheatmap(n,
         color = colorRampPalette(rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", 
                                        "#ffe6b7", "#aadce0", "#72bcd5", "#528fad", 
                                        "#376795", "#1e466e")))(100),
         annotation_row = rowcn,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         main = "MG Clusters SCENIC")
dev.off() # Close the device and save

# 2. Save as PNG (Using res = 300 for high quality)
png(file = png_file, width = 6, height = 12, units = "in", res = 300, bg = "white")
pheatmap(n,
         color = colorRampPalette(rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", 
                                        "#ffe6b7", "#aadce0", "#72bcd5", "#528fad", 
                                        "#376795", "#1e466e")))(100),
         annotation_row = rowcn,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         main = "MG Clusters SCENIC")
dev.off()


# Fig.1 AUCell Score-----------
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/SCENIC"))
# 1. library
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
library(viridis)

# 2.load MG data
tmp_obj <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))

# 3.load TF data
# TF活性均值
qload(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/SCENIC/results/for_rss_and_visual.qsm"))
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


#  rss查看特异的TF转录因子
library(openxlsx)
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellClusters[colnames(sub_regulonAUC), 
                                           selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss)
p1 <- rssPlot$plot

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

# 4.开始分基因绘图
# 颜色设定
pal <- viridis(n = 10, option = "D" )
tmp.regulon <- rownames(rowcn)
tmp.SeuObj <- tmp_obj

gene_list <- c("NFIX", "SOX10", "SOX2", "SREBF1", "ZKSCAN1")
for (i in gene_list) {
  j <- paste0(i,"(+)")
  # 4.开始分基因绘图
  # 颜色设定
  pal <- viridis(n = 10, option = "D" )
  tmp.SeuObj <- tmp_obj
  
  # 这里我暂时只用单个基因和regulon，才能跑通
  feature <- i
  tmp.regulon <- j
  
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
      # legend.title = element_blank(),
      # legend.position = "none",
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
      # legend.title = element_blank(),
      # legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) + 
    NoAxes() +
    ggtitle(paste0("AUCell Score of the ", gsub("\\(\\+\\)", "", tmp.regulon), " regulon")) +
    theme(plot.title = element_text(hjust = 0.5,face = "bold"))
  plot <- plot1 + plot2
  
  plot
  
  filename<-paste0("results/figure/","Fig.1_",'AUCell_', gsub("\\(\\+\\)", "", tmp.regulon),".v2",".pdf")
  ggsave(filename = filename,plot,width = 9,height = 3.5) # 加title的时候是3.8 不加title是3.
  
  filename<-paste0("results/figure/","Fig.1_",'AUCell_', gsub("\\(\\+\\)", "", tmp.regulon),".v2",".png")
  ggsave(filename = filename,plot,width = 9,height = 3.5) 
}


