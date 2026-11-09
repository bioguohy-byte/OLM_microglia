# 0. 加载包------------------------
library(DESeq2)
library(factoextra)
library(tinyarray)
library(dplyr)
library(tibble)
library(tidyr)
library(ggrepel)
library(ggprism)
library(pheatmap)

# 1. 数据加载---------------------
raw_counts <- read.csv(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/processed_data/All_gene_counts.list"),sep = '\t',header=TRUE,row.names=1)
raw_counts <- raw_counts[,c(-4,-8)]

# 4. 数据预处理与质控----------------
col_data <- read.csv(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/saved_data/sample_info.txt"),sep = '\t',header=TRUE,row.names = 1)
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = col_data,
  design = ~ group)

# 过滤低表达基因
keep <- rowSums(counts(dds) >= 10) >= 4 # 我们的过滤标准（至少4个样本中计数≥10）在统计效能和基因保留之间取得了良好平衡
dds <- dds[keep,]
cat("过滤后保留", nrow(dds), "个基因\n")

# 方差稳定变换
vsd <- vst(dds, blind = FALSE) # vst（方差稳定变换）：计算速度更快，特别适用于大样本量。rlog（正则化对数变换）：计算速度较慢，但对于小样本量可能表现稍好。
saveRDS(vsd,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/processed_data/vsd.rds"))

library(ggrepel)

pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

(pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2) + # 稍微放大一点点看得更清楚
  geom_text_repel(aes(label = name), show.legend = FALSE) + # 添加样本名并排斥重叠
  theme_bw(base_size = 12) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA plot") +
  scale_color_manual(values = c("control" = "blue", "treat" = "red")) + 
  theme(plot.title = element_text(hjust = 0.5)))
ggsave(Sys.getenv("MG_EXTERNAL_005"),width = 4,height = 3.5)

# 5. DESeq2差异分析--------------
# 运行DESeq2分析（包含标准化、离散度估计、模型拟合等所有步骤）
dds <- DESeq(dds)

# 提取分析结果
res <- results(dds, 
               contrast = c("group", "treat", "control"),  # 指定比较：处理组 vs 对照组
               alpha = 0.05,    # 显著性阈值，控制错误发现率(FDR)为5%
               lfcThreshold = 0) # log2FC阈值

# 添加显著性标记
# 注：这部分卡阈值从严格到宽松一共分为4种情况，差异基因数目递减：
# res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1   6+7
# res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 0   32+72
# res_df$pvalue < 0.05 & abs(res_df$log2FoldChange) >= 1 107+71
# res_df$pvalue < 0.05 & abs(res_df$log2FoldChange) >= 0 838+1061

res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 0, 
                             ifelse(res_df$log2FoldChange >= 0, "Up", "Down"), 
                             "Not significant")
res_df <- na.omit(res_df)

# 筛选显著差异表达基因
sig_genes <- subset(res_df[order(res_df$padj), ], padj < 0.05 & abs(log2FoldChange) > 0)
cat("显著差异基因数量:", nrow(sig_genes), "\n")
cat("上调基因数量:", sum(sig_genes$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("下调基因数量:", sum(sig_genes$log2FoldChange < 0, na.rm = TRUE), "\n")

# 保存结果
write.csv(res_df, Sys.getenv("MG_EXTERNAL_006"), row.names = TRUE)
res_df <- read.csv(Sys.getenv("MG_EXTERNAL_006"),row.names = 1)

# 6.结果可视化------------
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.5) +
  scale_color_manual(name = "",
                     values = c("Not significant" = "gray", 
                                "Up" = "brown", "Down" = "steelblue")) +
  xlim(c(-6, 6)) +
  # ylim(c(0, 10)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  # geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  labs(x="log2 (fold change)", y="-log10 (padj)") +
  ggtitle("volcano(treat vs control)") +
  theme_prism(border = TRUE) + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.title = element_blank()
  )

# 添加Top基因标签
p <- volcano_plot + 
  geom_label_repel(data = res_df[head(rownames(sig_genes), 10), ],
                   aes(label = rownames(res_df[head(rownames(sig_genes), 10), ])),
                   box.padding = 0.5, 
                   max.overlaps = Inf)
p
pdf(Sys.getenv("MG_EXTERNAL_007"),width = 7,height = 5)
p
dev.off()


# 选择前50个显著差异基因
select_genes <- rownames(sig_genes)[1:min(50, nrow(sig_genes))]
mat <- assay(vsd)[select_genes, ]

annotation_colors <- list(
  group = c(
    "Ct" = "#0072B2", 
    "C3" = "#D55E00" 
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

ggsave(Sys.getenv("MG_EXTERNAL_008"),p,width = 7,height = 6)
