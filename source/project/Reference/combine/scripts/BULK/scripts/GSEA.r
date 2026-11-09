rm(list = ls())  
options(stringsAsFactors = F)
# library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggplot2)
library(GseaVis)


# 物种设置----------------------
organism = 'hsa'    #  人类'hsa' 小鼠'mmu'   
OrgDb = 'org.Hs.eg.db'#人类"org.Hs.eg.db" 小鼠"org.Mm.eg.db"

# 处理deseq2差异结果---------------
res_df<-read.csv(Sys.getenv("MG_EXTERNAL_009"),row.names = 1)
need_DEG <- res_df
need_DEG <- need_DEG[,c(2,5)] 
colnames(need_DEG) <- c('log2FoldChange','pvalue')
need_DEG_clean <- need_DEG[!is.na(need_DEG$log2FoldChange) & !is.na(need_DEG$pvalue), ]
need_DEG_clean$SYMBOL <- rownames(need_DEG_clean)

# 创建gsea分析的geneList（包含从大到小排列的log2FoldChange和ENTREZID信息）----------
#转化id  
df <- bitr(rownames(need_DEG_clean), 
           fromType = "SYMBOL",
           toType =  "ENTREZID",
           OrgDb = OrgDb) #人数据库org.Hs.eg.db 小鼠org.Mm.eg.db
need_DEG_clean <- merge(need_DEG_clean, df, by='SYMBOL')  #按照SYMBOL合并注释信息
geneList <- need_DEG_clean$log2FoldChange
names(geneList) <- need_DEG_clean$ENTREZID
geneList <- sort(geneList, decreasing = T)   #从大到小排序


##### gsea富集 ####
KEGG_kk_entrez <- gseKEGG(geneList     = geneList,
                          organism     = organism, #人hsa 鼠mmu
                          pvalueCutoff = 0.25)   
KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez, 
                             OrgDb=OrgDb,
                             keyType='ENTREZID')#转化id             
saveRDS(KEGG_kk,file = Sys.getenv("MG_EXTERNAL_010"))
KEGG_kk <- readRDS(Sys.getenv("MG_EXTERNAL_010"))
write.csv(KEGG_kk@result,file = 'results/table/KEGG_gsea_all.xlsx')

GO_kk_entrez <- gseGO(geneList     = geneList,
                      ont          = "BP",  # "BP"、"MF"和"CC"或"ALL"
                      OrgDb        = OrgDb,#人类org.Hs.eg.db 鼠org.Mm.eg.db
                      keyType      = "ENTREZID",
                      pvalueCutoff = 0.25)   
GO_kk <- DOSE::setReadable(GO_kk_entrez, 
                           OrgDb=OrgDb,
                           keyType='ENTREZID')#转化id 
saveRDS(GO_kk,file = Sys.getenv("MG_EXTERNAL_011"))
write.csv(KEGG_kk@result,file = 'results/table/GO_gsea_all.xlsx')

# KEGG------------------
kegg_gsea <- KEGG_kk
kegg_gsea_entrez <- KEGG_kk_entrez

#一般认为|NES|>1，NOM pvalue<0.05，FDR（padj）<0.25的通路是显著富集的
kegg_gsea_cut <- kegg_gsea[kegg_gsea$pvalue<0.05 & kegg_gsea$p.adjust<0.25 & abs(kegg_gsea$NES)>1]
kegg_gsea_cut_down <- kegg_gsea_cut[kegg_gsea_cut$NES < 0,]
kegg_gsea_cut_up <- kegg_gsea_cut[kegg_gsea_cut$NES > 0,]

write.csv(kegg_gsea_cut,file = 'results/table/KEGG_gsea_filter.xlsx')
write.csv(kegg_gsea_cut_down,file = 'results/table/KEGG_gsea_down.xlsx')
write.csv(kegg_gsea_cut_up,file = 'results/table/KEGG_gsea_up.xlsx')

#选择展现NES前12个通路 
kegg_down_gsea <- kegg_gsea_cut_down[tail(order(kegg_gsea_cut_down$NES,decreasing = T),12),]
kegg_up_gsea <- kegg_gsea_cut_up[head(order(kegg_gsea_cut_up$NES,decreasing = T),12),]
kegg_diff_gsea <- kegg_gsea_cut[head(order(abs(kegg_gsea_cut$NES),decreasing = T),12),]

# # 上调合并GSEA图
# (kegg_gsea_up_plot <- gseaNb(
#   object = kegg_gsea,
#   geneSetID = kegg_up_gsea$ID,
#   curveCol= jjAnno::useMyCol('stallion', length(kegg_up_gsea$ID)),
#   # curveCol= col1,
#   subPlot = 2,
#   addPval = T,
#   pvalX = 1,
#   termWidth = 100,
#   pvalY = 1)
# )
# 
# # 下调合并GSEA图
# (kegg_gsea_down_plot <- gseaNb(
#   object = kegg_gsea,
#   geneSetID = kegg_down_gsea$ID,
#   curveCol= jjAnno::useMyCol('stallion', length(kegg_down_gsea$ID)),
#   # curveCol= col1,
#   subPlot = 2,
#   addPval = T,
#   termWidth = 100,
#   pvalX = 0.05,pvalY = 0.05,
#   legend.position = c(0.8,0.8),
#   base_size = 15)
# )

# up分别画图
lapply(kegg_up_gsea$ID, function(x){
  gseaNb(object = kegg_gsea,
         geneSetID = x,
         addPval = T,
         pvalX = 0.7,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,
         subPlot = 2)
}) -> kegg_up_gseaList

# combine
(kegg_gsea_up_sep_plot <- cowplot::plot_grid(plotlist = kegg_up_gseaList,ncol = 3,align = 'hv'))

# save
pdf(Sys.getenv("MG_EXTERNAL_012"),width = 15,height = 12)
kegg_gsea_up_sep_plot
dev.off()


# down分别画图
lapply(kegg_down_gsea$ID, function(x){
  gseaNb(object = kegg_gsea,
         geneSetID = x,
         addPval = T,
         pvalX = 0.7,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,
         subPlot = 2)
}) -> kegg_down_gseaList

# combine
(kegg_gsea_down_sep_plot <- cowplot::plot_grid(plotlist = kegg_down_gseaList,ncol = 3,align = 'hv'))

# save
pdf(Sys.getenv("MG_EXTERNAL_013"),width = 15,height = 12)
kegg_gsea_down_sep_plot
dev.off()



# GO------------------
go_gsea <- GO_kk
go_gsea_entrez <- GO_kk_entrez

#一般认为|NES|>1，NOM pvalue<0.05，FDR（padj）<0.25的通路是显著富集的
go_gsea_cut <- go_gsea[go_gsea$pvalue<0.05 & go_gsea$p.adjust<0.25 & abs(go_gsea$NES)>1]
go_gsea_cut_down <- go_gsea_cut[go_gsea_cut$NES < 0,]
go_gsea_cut_up <- go_gsea_cut[go_gsea_cut$NES > 0,]

write.csv(go_gsea_cut,file = 'results/table/go_gsea_filter.xlsx')
write.csv(go_gsea_cut_down,file = 'results/table/go_gsea_down.xlsx')
write.csv(go_gsea_cut_up,file = 'results/table/go_gsea_up.xlsx')

#选择展现NES前12个通路 
go_down_gsea <- go_gsea_cut_down[tail(order(go_gsea_cut_down$NES,decreasing = T),12),]
go_up_gsea <- go_gsea_cut_up[head(order(go_gsea_cut_up$NES,decreasing = T),12),]
go_diff_gsea <- go_gsea_cut[head(order(abs(go_gsea_cut$NES),decreasing = T),12),]

# # 上调合并GSEA图
# (go_gsea_up_plot <- gseaNb(
#   object = go_gsea,
#   geneSetID = go_up_gsea$ID,
#   curveCol= jjAnno::useMyCol('stallion', length(go_up_gsea$ID)),
#   # curveCol= col1,
#   subPlot = 2,
#   addPval = T,
#   pvalX = 1,
#   termWidth = 100,
#   pvalY = 1)
# )
# 
# # 下调合并GSEA图
# (go_gsea_down_plot <- gseaNb(
#   object = go_gsea,
#   geneSetID = go_down_gsea$ID,
#   curveCol= jjAnno::useMyCol('stallion', length(go_down_gsea$ID)),
#   # curveCol= col1,
#   subPlot = 2,
#   addPval = T,
#   pvalX = 1,
#   termWidth = 100,
#   pvalY = 1)
# )

# up分别画图
lapply(go_up_gsea$ID, function(x){
  gseaNb(object = go_gsea,
         geneSetID = x,
         addPval = T,
         pvalX = 0.6,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,
         subPlot = 2)
}) -> go_up_gseaList

# combine
(go_gsea_up_sep_plot <- cowplot::plot_grid(plotlist = go_up_gseaList,ncol = 3,align = 'hv'))

# save
pdf(Sys.getenv("MG_EXTERNAL_014"),width = 15,height = 12)
go_gsea_up_sep_plot
dev.off()


# down分别画图
lapply(go_down_gsea$ID, function(x){
  gseaNb(object = go_gsea,
         geneSetID = x,
         addPval = T,
         pvalX = 0.6,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,
         subPlot = 2)
}) -> go_down_gseaList

# combine
(go_gsea_down_sep_plot <- cowplot::plot_grid(plotlist = go_down_gseaList,ncol = 3,align = 'hv'))

# save
pdf(Sys.getenv("MG_EXTERNAL_015"),width = 15,height = 12)
go_gsea_down_sep_plot
dev.off()

##### custom gene set #######
# ==========================================
# 使用自定义基因集进行 GSEA 分析
# ==========================================

# 1. 假设这是你的自定义基因集（SYMBOL 格式）
# 请把这里的基因替换为你自己的目标基因
my_custom_genes <- c(
  "IL1RAPL1", "PCDH9", "CTNNA3", "ST18", "RNF220", "SLC24A2", "PPP2R2B", "FRMD5", "MAGI2", "TMTC2",
  "NCAM2", "TTLL7", "EDIL3", "PLP1", "UNC5C", "PDE1C", "SLC44A1", "PEX5L", "SPOCK3", "NKAIN2",
  "PRUNE2", "MAP7", "C10orf90", "PTPRD", "SPOCK1", "AC026316.5", "PHLPP1", "ZNF536", "MOBP", "ANK3",
  "PTPRK", "KLHL32", "NCAM1", "CADM2", "TRIM2", "KIRREL3", "CLMN", "TF", "DNAJC6", "KCNH8",
  "CNTNAP4", "CDH20", "LIMCH1", "SH3GL3", "NRXN3", "GPM6B", "PALM2", "NLGN1", "FBXL7", "GRM3",
  "AGAP1", "DNM3", "AKAP6", "ENPP2", "CREB5", "TMEM144", "ERBB4", "DLG2", "NTM", "APBB2",
  "ENOX1", "TMEFF2", "PLEKHH1", "SLC5A11", "DLC1", "SLAIN1", "ANO4", "NFASC", "BCAS1", "GRID1",
  "PCSK6", "SLCO3A1", "DSCAML1", "COBL", "AK5", "CERCAM", "UGT8", "POLR2F", "SCD5", "OTUD7A",
  "COL4A5", "LINC01608", "ABCA2", "PDE4B", "HHIP", "SLCO1A2", "SHROOM4", "DAAM2", "MBP", "SYNJ2",
  "NPAS3", "ATP8A1", "CDC42BPA", "SIK3", "DOCK1", "ABCA8", "LINC00609", "RAPGEF5", "PLCL1", "GNAO1",
  "CNTN2", "NALCN", "PXK", "CNDP1", "TJP1", "CNP", "ANK2", "USP54", "SH3D19", "PIP4K2A", "DOCK9",
  "FGFR2", "FA2H", "ARAP2", "HIP1", "PHACTR3", "ATP10B", "KCNMB4", "RFTN2", "KIF6", "DOCK5",
  "CACNA2D1", "AC012593.1", "HS3ST5", "MAN2A1", "PKP4", "FAM107B", "AMER2", "DIP2C", "COLGALT2", "PHLDB1",
  "HECW2", "PIEZO2", "KANK1", "HDAC2-AS2", "LINC00639", "ANLN", "CTNND2", "QDPR", "ARFGEF3", "LINC01170",
  "FOLH1", "ANKS1B", "APLP1", "MYO1D", "PTGDS", "PDE1A", "SEPT7", "FUT8", "FAM13C", "MEGF10",
  "KCTD8", "LAMA2", "PHYHIPL", "RNF144A", "LANCL1", "ERBIN", "GALNT13", "GPRC5B", "LINC00320", "CD22",
  "MYRF", "CLASP2", "CDK18", "NEO1", "TMEM165", "CDH19", "DLG1", "DMD", "MOG", "JAM3", "LRP2",
  "ZNF704", "IGSF11", "CLDND1", "SCD", "LRP1B", "PPM1H", "TSPAN5", "MAPT", "CNTNAP2", "SEMA6A",
  "FRYL", "MYO6", "LPGAT1", "LSAMP", "PPP1R9A", "MAG", "EFHD1", "AGPAT4", "AATK", "DOCK3",
  "RASGRF2", "KLHL4", "ITGA2", "GAB1", "SASH1", "LMCD1-AS1", "LIFR", "CADM1"
)

# 2. 将你的自定义基因集转换成 ENTREZID (为了和你之前的 geneList 匹配)
my_custom_df <- bitr(my_custom_genes, 
                     fromType = "SYMBOL",
                     toType =  "ENTREZID",
                     OrgDb = OrgDb)

# 3. 构建 TERM2GENE 数据框
# 假设我们把这个自定义通路命名为 "My_Custom_Pathway"
custom_term2gene <- data.frame(
  term = rep("My_Custom_Pathway", nrow(my_custom_df)), # 通路名称
  gene = my_custom_df$ENTREZID                         # 对应的ENTREZID
)

# 4. 运行自定义 GSEA 分析
# 使用核心函数 GSEA()
custom_gsea_res <- GSEA(geneList = geneList,           # 你之前排好序的 geneList
                        TERM2GENE = custom_term2gene,  # 你的自定义基因集
                        pvalueCutoff = 0.25,           # p值阈值
                        pAdjustMethod = "BH")

# 将结果转换为 SYMBOL 方便阅读（可选）
custom_gsea_res_readable <- DOSE::setReadable(custom_gsea_res, 
                                              OrgDb = OrgDb, 
                                              keyType = 'ENTREZID')

# 查看结果表格
View(custom_gsea_res_readable@result)
write.csv(custom_gsea_res_readable@result, file = 'results/table/Custom_gsea_res.csv')

# 5. 使用 GseaVis 画图
# 画图逻辑和你前面画 KEGG/GO 是一样的
library(GseaVis)

# 画图并保存
gseaNb(object = custom_gsea_res,
       geneSetID = "My_Custom_Pathway", # 这里填你自定义的 term 名称
       addPval = T,
       pvalX = 0.6, pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,
       subPlot = 2)
dev.off()



##### heatmap #####
# 选择前50个显著差异基因
select_genes <- c(
  "FAM13C", "IL1RAPL1", "ITGA2", "ATP8A1", "HS3ST5", "KIF6", 
  "SPOCK3", "SEMA6A", "AK5", "PCSK6", "CD22", "PCDH9", 
  "PEX5L", "RFTN2", "DAAM2", "NTM", "GAB1", "ABCA2", 
  "SPOCK1", "DOCK9", "RNF144A", "PLEKHH1", "SYNJ2", 
  "FAM107B", "AGAP1", "GPM6B", "TTLL7", "TRIM2", "CLMN", 
  "DMD", "KANK1", "DNM3"
)
mat <- assay(vsd)[select_genes, ]

annotation_colors <- list(
  group = c(
    "control" = "#0072B2", 
    "treat" = "#D55E00" 
  )
)

# 赋值给 p
p <- pheatmap(mat,
              scale = "row",
              show_rownames = TRUE,
              show_colnames = TRUE,
              cluster_rows = TRUE,
              cluster_cols = FALSE,
              annotation_col = as.data.frame(col_data["group"]),
              annotation_colors = annotation_colors,
              main = "Heatmap of Significant DEGs",
              fontsize_row = 8,
              color = colorRampPalette(c("#2A788E", "#FFFFFF", "#D8576B"))(50),
              breaks = seq(-2, 2, length.out = 51)
)

p

ggsave(Sys.getenv("MG_EXTERNAL_016"),p,width = 7,height = 6)
