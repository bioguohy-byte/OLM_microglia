# SCENIC-------------------
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
library(openxlsx)
library(ggplot2)
library(pheatmap)
library(viridis)
library(patchwork)

annotation_colors <- list(
  cluster = c("MG0" = "#4DAF4A", "MG1" = "#A65628", "MG2" = "#377EB8", 
              "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG5" = "#8B0000")
)

loom <- open_loom('./results/out_SCENIC.loom') 
tmp_obj <- qread('../../data/processed_data/mg_cca_meta_v2.qs')

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)  
close_loom(loom)

sub_regulonAUC <- regulonAUC[, match(colnames(tmp_obj), colnames(regulonAUC))]

cellClusters <- data.frame(row.names = colnames(tmp_obj), 
                           sub_clusters = as.character(tmp_obj$CellType1))

qsavem(sub_regulonAUC, cellClusters, tmp_obj, file = "./results/for_rss_and_visual.qsm")

selectedResolution <- "sub_clusters"
cellsPerGroup <- split(rownames(cellClusters), cellClusters[,selectedResolution]) 

sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 

regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 

regulonActivity_byGroup_Scaled <- na.omit(regulonActivity_byGroup_Scaled)

# Fig.01 ---------------
pheatmap(regulonActivity_byGroup_Scaled)

rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellClusters[colnames(sub_regulonAUC), selectedResolution]) 
rss <- na.omit(rss) 
rssPlot <- plotRSS(rss)

# Fig.02 ---------------
p1 <- rssPlot$plot

rss_scaled <- regulonActivity_byGroup_Scaled

df <- do.call(rbind,
              lapply(1:ncol(rss_scaled), function(i){
                dat <- data.frame(
                  path  = rownames(rss_scaled),
                  cluster = colnames(rss_scaled)[i],
                  sd.1 = rss_scaled[,i],
                  sd.2 = apply(rss_scaled[,-i], 1, median)  
                )
              }))
df$fc <- df$sd.1 - df$sd.2

top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
top5 <- top5[!duplicated(top5$path), ]

rowcn <- data.frame(cluster = top5$cluster) 
n <- rss_scaled[top5$path,] 
rownames(rowcn) <- top5$path  
rownames(n) <- top5$path

write.xlsx(df, './results/SCENIC.xlsx', overwrite = T)

base_path <- "./results/figure/"
pdf_file <- paste0(base_path, "SCENIC_heatmap.pdf")
png_file <- paste0(base_path, "SCENIC_heatmap.png")

# Fig.03 ---------------
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
               filename = pdf_file, 
               width = 6, 
               height = 12)

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
dev.off() 

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

gene_list <- c("NFIX", "SOX10", "SOX2", "SREBF1", "ZKSCAN1")
pal <- viridis(n = 10, option = "D" )

# Fig.04 ---------------
for (i in gene_list) {
  j <- paste0(i,"(+)")
  
  tmp.SeuObj <- tmp_obj
  feature <- i
  tmp.regulon <- j
  
  plot.data <- FeaturePlot(tmp.SeuObj, features = feature, order = T, cols = pal, raster = F)$data
  fdata <- FetchData(object = tmp.SeuObj, vars = feature)
  plot.data$feature <- fdata[match(rownames(plot.data), rownames(fdata)), ]
  
  plot1 <- ggplot(plot.data, aes(UMAP_1, UMAP_2, color = feature)) + 
    geom_point(size = 0.001) +
    scale_color_gradientn(colors = pal) +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) + 
    NoAxes() +
    ggtitle(paste0(gsub("\\(\\+\\)", "", tmp.regulon), " gene expression"))
  
  tmp.df <- data.frame(tmp.regulon = sub_regulonAUC@assays@data@listData[["AUC"]][tmp.regulon, ]) 
  tmp.SeuObj@meta.data$tmp.regulon <- tmp.df$tmp.regulon
  
  plot.data_auc <- FeaturePlot(tmp.SeuObj, features = "tmp.regulon", order = T, cols = pal, raster = F)$data
  plot.data_auc$tmp.regulon <- tmp.df[match(rownames(plot.data_auc), rownames(tmp.df)), ]
  
  plot2 <- ggplot(plot.data_auc, aes(UMAP_1, UMAP_2, color = tmp.regulon)) +
    geom_point(size = 0.001) +
    scale_color_gradientn(colors = pal) +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) + 
    NoAxes() +
    ggtitle(paste0("AUCell Score of the ", gsub("\\(\\+\\)", "", tmp.regulon), " regulon"))
  
  plot_final <- plot1 + plot2
  
  filename_pdf <- paste0(base_path, "Fig.4_AUCell_", gsub("\\(\\+\\)", "", tmp.regulon), ".v2.pdf")
  ggsave(filename = filename_pdf, plot_final, width = 9, height = 3.5) 
  
  filename_png <- paste0(base_path, "Fig.4_AUCell_", gsub("\\(\\+\\)", "", tmp.regulon), ".v2.png")
  ggsave(filename = filename_png, plot_final, width = 9, height = 3.5) 
}