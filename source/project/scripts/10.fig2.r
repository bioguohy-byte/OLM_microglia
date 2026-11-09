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
library(harmony)
library(qs)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(enrichplot)
library(GseaVis)

# 多线程
library(future)
plan("multisession", workers = 6)
options(future.globals.maxSize = 100 * 1024^3)

# 加载数据--------
tmp_obj <- qread('data/processed_data/14.mg_combined_meta.qs')

# 小胶质myelination基因集热图heatmap--------------
gene_list <- c("ST18", "SLC44A1", "PHLPP1", "PTPRD", "NCAM2", "PLP1", "NCAM1", "KCNMB4",
               "GNAO1", "UGT8", "MOBP", "MYRF", "BCAS1", "MBP", "MOG", "OPALIN",
               "QKI", "CLDN11", "FA2H", "DLG1", "SH3TC2", "PTN")
Exp_av <- AverageExpression(tmp_obj,
                            group.by = "sub_clusters",
                            assays = "RNA")

Exp_av <- Exp_av[[1]]
Exp_av <- as.data.frame(Exp_av)
Exp_av$gene <- rownames(Exp_av)
result <- Exp_av %>%
  filter(gene %in% gene_list)

result <- result[match(gene_list, result$gene),]
#现在的result中有一列是gene列，所以不能直接放在pheatmap中进行绘图，#因为里面不全是数字了，所以要删除那一列
result <- subset(result, select = -gene)

pdf("results/figure/MG_myelination_gene_heatmap.pdf", width = 3, height = 3)
pheatmap(result, cluster_cols = FALSE, cluster_rows = FALSE,scale = "row",show_rownames = T,legend = F)
dev.off()

# cluster_cols = FALSE: 这个参数决定是否对列进行聚类。如果设置为FALSE，那么列就不会被聚类。
# cluster_rows = FALSE: 这个参数决定是否对行进行聚类。如果设置为FALSE，那么行就不会被聚类。
# scale = "row": 这个参数决定是否对数据进行缩放，以及如何进行缩放。如果设置为"row"，那么每一行的数据都会被缩放，使得每一行的平均值为0，标准差为1。这样可以帮助你更好地比较不同基因的表达模式。

# 小胶质其他功能基因集热图heatmap--------------
gene_list <- c("CTNNA3", "CNTNAP4", "SPOCK1", "CADM2", "TJP1", "CNTN2", "CTNND2",
               "RNF220", "C10orf90", "FBXL7", "TRIM2", "PPM1H", "RNF144A", 
               "UNC5C","TTLL7", "MAGI2", "NTM", "NFASC", "ENPP2", 
               "DNAJC6", "DLC1", "LPGAT1")

Exp_av <- AverageExpression(tmp_obj,
                            group.by = "sub_clusters",
                            assays = "RNA")

Exp_av <- Exp_av[[1]]
Exp_av <- as.data.frame(Exp_av)
Exp_av$gene <- rownames(Exp_av)
result <- Exp_av %>%
  filter(gene %in% gene_list)

result <- result[match(gene_list, result$gene),]
#现在的result中有一列是gene列，所以不能直接放在pheatmap中进行绘图，#因为里面不全是数字了，所以要删除那一列
result <- subset(result, select = -gene)

pdf("results/figure/MG_other_function_gene_heatmap.pdf", width = 3, height = 3)
pheatmap(result, cluster_cols = FALSE, cluster_rows = FALSE,scale = "row",show_rownames = T,legend = F)
dev.off()

# cluster_cols = FALSE: 这个参数决定是否对列进行聚类。如果设置为FALSE，那么列就不会被聚类。
# cluster_rows = FALSE: 这个参数决定是否对行进行聚类。如果设置为FALSE，那么行就不会被聚类。
# scale = "row": 这个参数决定是否对数据进行缩放，以及如何进行缩放。如果设置为"row"，那么每一行的数据都会被缩放，使得每一行的平均值为0，标准差为1。这样可以帮助你更好地比较不同基因的表达模式。

# 针对M0的GSEA作图--------------
dir_data <- "scripts/GSEA/data/"
dir_res  <- "scripts/GSEA/results/"
if(!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE)

# 2. 基因集 (GMT) 准备与合并

# 2.1 处理自定义基因集 (从 Excel 到 GMT)
# 注意：Excel 第三列起为基因，第二列为描述
gene_sets_raw <- read.xlsx(paste0(dir_data, 'my_gene_set.xlsx'), colNames = FALSE, na.strings = "")
write.table(gene_sets_raw, paste0(dir_data, 'my_genesets.gmt'), 
            sep = "\t", row.names = F, col.names = F, quote = F)

custom_gmt <- read.gmt(paste0(dir_data, 'my_genesets.gmt'))
custom_gmt$gene <- toupper(custom_gmt$gene) # 强制转大写以匹配 DEG

# 2.2 读取并筛选标准 GO 基因集
# 筛选你关注的髓鞘及少突细胞相关条目
gmt_go <- read.gmt(paste0(dir_data, "c5.go.v2023.2.Hs.symbols.gmt"))
target_go_terms <- c(
  'GOBP_POSITIVE_REGULATION_OF_MYELINATION', 'GOBP_REGULATION_OF_MYELINATION',
  'GOBP_OLIGODENDROCYTE_DIFFERENTIATION', 'GOBP_CELL_JUNCTION_ORGANIZATION',
  'GOBP_ADAPTIVE_IMMUNE_RESPONSE', 'GOBP_FATTY_ACID_TRANSPORT',
  'GOBP_MYELIN_ASSEMBLY', 'GOBP_OLIGODENDROCYTE_DEVELOPMENT',
  'GOBP_LIPID_CATABOLIC_PROCESS'
)

filtered_go <- gmt_go %>%
  filter(term %in% target_go_terms) %>%
  mutate(term = gsub("GOBP_", "", term)) %>%
  mutate(term = tolower(gsub("_", " ", term)))

# 2.3 合并自定义与筛选后的 GO 基因集
combined_gmt <- rbind(custom_gmt, filtered_go)

# 3. GSEA 计算 (以 M0 亚群为例) 

# 3.1 准备排序基因列表 (Ranked Gene List)
Idents(tmp_obj) <- 'sub_clusters'
sub_clusters_deg_genes <- FindAllMarkers(tmp_obj,only.pos = F, min.pct = 0.1, logfc.threshold = 0) #差异基因不能只选positive的

# 提取 M0 的 DEG 并按 logFC 降序排列
rt_m0 <- sub_clusters_deg_genes %>%
  filter(cluster == 'M0') %>%
  arrange(desc(avg_log2FC)) %>%
  distinct(gene, .keep_all = TRUE) # 确保基因唯一

# 构建 GSEA 输入向量
logFC_vector <- rt_m0$avg_log2FC
names(logFC_vector) <- toupper(rt_m0$gene) # 确保基因名大写以匹配 GMT

# 3.2 运行 GSEA
# minGSSize 设为 1 是为了兼容较小的自定义基因集
kk <- GSEA(
  geneList = logFC_vector, 
  TERM2GENE = combined_gmt, 
  pvalueCutoff = 1, # 先保留所有结果，后续手动筛选
  pAdjustMethod = 'BH',
  minGSSize = 1,
  maxGSSize = 500,
  scoreType = "pos", # 仅关注正向富集
  eps = 0
)

qsave(kk, file = paste0(dir_data, 'GSEA_kk_M0.qs'))
write.table(as.data.frame(kk), 
            file = paste0(dir_res, "combined_genesets_M0_GSEA_result.txt"), 
            sep = "\t", quote = F, row.names = F)

# 4. 结果可视化 (GseaVis)

# 选择你想要在图中展示的特定条目
selected_ids <- c(
  'oligodendrocyte development',
  'oligodendrocyte differentiation',
  'myelin assembly',
  'cell junction organization',
  'positive regulation of myelination'
)

# 配色方案 (对应 selected_ids 的数量)
colors_set <- c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D23")

# 绘图
p_gsea <- gseaNb(
  object = kk, 
  geneSetID = selected_ids,
  subPlot = 3, 
  curveCol = colors_set,
  addPval = TRUE,
  pvalSize = 5,
  base_size = 15
)

# 保存 PDF
pdf(file = paste0(dir_res, "M0_combined_GSEA_select_plot.pdf"), width = 12, height = 8)
print(p_gsea)
dev.off()

# 绘图
p_gsea_2 <- gseaNb(
  object = kk, 
  geneSetID = selected_ids,
  subPlot = 3, 
  curveCol = colors_set,
  addPval = FALSE,
  pvalSize = 5,
  base_size = 15
)

# 保存 PDF
pdf(file = paste0(dir_res, "M0_combined_GSEA_select_plot_no_pvalue.pdf"), width = 12, height = 8)
print(p_gsea_2)
dev.off()

