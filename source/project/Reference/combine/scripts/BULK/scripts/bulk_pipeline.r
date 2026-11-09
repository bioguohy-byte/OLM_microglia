# MODULE : Bulk RNA-seq Differential Expression Analysis Pipeline -------------------
# Description: A standardized and reusable pipeline for bulk RNA-seq downstream 
#              analysis (Differential expression, PCA, Volcano Plot, Heatmap).
# Date: 2026-03-18

suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(ggrepel))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))

# 1 Parameters --------------------------------------------------------------
# 在此处修改参数，以适配不同的项目和数据
INPUT_FILE <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/processed_data/All_gene_counts.list")           # 输入的表达矩阵文件路径 (请修改为你的实际路径)
OUTPUT_DIR <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/results/DESeq2_Analysis_Results")  # 结果保存的文件夹名称

# 实验设计参数
CONTROL_KEYWORD <- "con"                 # 对照组样本列名的特征词 (例如: con1, con2)
TREAT_KEYWORD   <- "treat"               # 处理组样本列名的特征词 (例如: treat1, treat2)
CONTROL_NAME    <- "Control"             # 报告和图表中对照组的显示名称
TREAT_NAME      <- "Treatment"           # 报告和图表中处理组的显示名称

# 差异表达筛选阈值
PADJ_CUTOFF   <- 0.05                    # 显著性阈值 (FDR/p.adjust)
LOG2FC_CUTOFF <- 1.0                     # log2FoldChange的绝对值阈值 (1表示差异达到2倍)


# 2 Process Data ------------------------------------------------------------
cat("[INFO] Starting Process Data...\n")
# 创建输出目录
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# 2.1 读取表达矩阵
# 注意: 文件第一列列名为 '#ID'，R默认将 '#' 视为注释，因此必须设置 comment.char = ""
counts_data <- read.table(INPUT_FILE, header = TRUE, row.names = 1, 
                          sep = "\t", check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)

# 清理可能存在的小数点或缺失值 (确保counts是整数)
counts_data <- round(counts_data)
counts_data <- counts_data[rowSums(counts_data) > 0, ] # 过滤掉所有样本表达量均为0的基因

# 2.2 构建样本元信息 (colData)
sample_names <- colnames(counts_data)
group_list <- rep(NA, length(sample_names))
group_list[grepl(CONTROL_KEYWORD, sample_names, ignore.case = TRUE)] <- CONTROL_NAME
group_list[grepl(TREAT_KEYWORD, sample_names, ignore.case = TRUE)] <- TREAT_NAME

# 检查是否有未匹配到分组的样本
if(any(is.na(group_list))) {
  stop("Error: Some samples do not match the CONTROL_KEYWORD or TREAT_KEYWORD. Please check your parameters.")
}

col_data <- data.frame(
  row.names = sample_names,
  condition = factor(group_list, levels = c(CONTROL_NAME, TREAT_NAME)) # 确保Control为对照参考水平
)

# 2.3 构建DESeq2对象并运行分析
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = col_data,
                              design = ~ condition)

# 过滤低表达基因 (至少在一个组的样本数量以上中有>10的reads)
min_samples <- min(table(col_data$condition))
keep <- rowSums(counts(dds) >= 10) >= min_samples
dds <- dds[keep, ]

# 运行DESeq分析
dds <- DESeq(dds)

# 提取差异分析结果
res <- results(dds, contrast = c("condition", TREAT_NAME, CONTROL_NAME))
res_df <- as.data.frame(res) %>% rownames_to_column(var = "GeneID")

# 处理NA值
res_df$padj[is.na(res_df$padj)] <- 1
res_df$pvalue[is.na(res_df$pvalue)] <- 1

# 添加显著性标签 (Significant)
res_df <- res_df %>%
  mutate(
    Significance = case_when(
      padj < PADJ_CUTOFF & log2FoldChange > LOG2FC_CUTOFF ~ "Up",
      padj < PADJ_CUTOFF & log2FoldChange < -LOG2FC_CUTOFF ~ "Down",
      TRUE ~ "Not Sig"
    )
  )

cat(sprintf("[INFO] Differential expression summary: %d Up, %d Down\n", 
            sum(res_df$Significance == "Up"), 
            sum(res_df$Significance == "Down")))


# 定义比较方向的名称前缀，方便后续文件命名 (例如: Treatment_vs_Control)
COMP_PREFIX <- paste0(TREAT_NAME, "_vs_", CONTROL_NAME)

# 3 Generate Plots ----------------------------------------------------------
cat("[INFO] Generating Plots for", COMP_PREFIX, "...\n")
# 对数据进行标准化转换 (用于PCA和Heatmap)
vsd <- vst(dds, blind = FALSE)

# 3.1 PCA图 (Principal Component Analysis)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA Plot") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 3.2 火山图 (Volcano Plot)
# 选取差异最大的Top 10基因进行标记
top_genes <- res_df %>% filter(Significance != "Not Sig") %>% top_n(10, wt = abs(log2FoldChange))

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "black") +
  # 注意：如果你合并了上一步提供的 ID 转换代码，这里可以将 label = GeneID 改为 label = Symbol
  geom_text_repel(data = top_genes, aes(label = GeneID), color = "black", size = 3, max.overlaps = 20) +
  theme_bw() +
  ggtitle(paste0("Volcano Plot: ", COMP_PREFIX)) + # 标题动态显示是谁比谁
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

# 3.3 差异基因热图 (Heatmap)
sig_genes <- res_df %>% filter(Significance != "Not Sig") %>% pull(GeneID)

if (length(sig_genes) > 1) {
  # 选取显著差异基因的标准化表达量
  mat <- assay(vsd)[sig_genes, ]
  
  # 对数据进行Z-score标准化处理 (按行标准化)
  mat_scaled <- t(scale(t(mat)))
  
  # 热图注释
  anno_col <- data.frame(Group = col_data$condition)
  rownames(anno_col) <- rownames(col_data)
  
  # 如果差异基因超过50个，不显示基因名
  show_rownames_flag <- ifelse(length(sig_genes) <= 50, TRUE, FALSE)
  
  p_heatmap <- pheatmap(mat_scaled, 
                        annotation_col = anno_col, 
                        show_rownames = show_rownames_flag,
                        cluster_cols = TRUE,
                        scale = "none",
                        main = paste0("Significant DEGs: ", COMP_PREFIX), # 标题动态显示是谁比谁
                        silent = TRUE)
} else {
  cat("[WARNING] Not enough significant genes to generate heatmap.\n")
}


# 4 Save Results ------------------------------------------------------------
cat("[INFO] Saving Results to", OUTPUT_DIR, "...\n")

# 4.1 保存表格结果 (文件名全部加上比较方向前缀)
# 4.1.1 导出所有DESeq2结果
write.csv(res_df, file.path(OUTPUT_DIR, paste0("DESeq2_all_", COMP_PREFIX, ".csv")), row.names = FALSE)

# 4.1.2 导出显著差异基因列表
sig_res_df <- res_df %>% filter(Significance != "Not Sig") %>% arrange(padj)
write.csv(sig_res_df, file.path(OUTPUT_DIR, paste0("DESeq2_sig_DEGs_", COMP_PREFIX, ".csv")), row.names = FALSE)

# 4.1.3 导出Count标准化矩阵 (Normalized Counts) (归一化矩阵不涉及差异方向，故无需加前缀)
norm_counts <- counts(dds, normalized = TRUE) %>% as.data.frame() %>% rownames_to_column(var = "GeneID")
write.csv(norm_counts, file.path(OUTPUT_DIR, "DESeq2_normalized_counts.csv"), row.names = FALSE)

# 4.2 保存图表 (文件名全部加上比较方向前缀)
# 保存PCA (PCA是对整体数据的降维，不涉及方向，保持原名)
ggsave(filename = file.path(OUTPUT_DIR, "PCA_plot.pdf"), plot = p_pca, width = 6, height = 5)
ggsave(filename = file.path(OUTPUT_DIR, "PCA_plot.png"), plot = p_pca, width = 6, height = 5, dpi = 300)

# 保存火山图
ggsave(filename = file.path(OUTPUT_DIR, paste0("Volcano_plot_", COMP_PREFIX, ".pdf")), plot = p_volcano, width = 7, height = 6)
ggsave(filename = file.path(OUTPUT_DIR, paste0("Volcano_plot_", COMP_PREFIX, ".png")), plot = p_volcano, width = 7, height = 6, dpi = 300)

# 保存热图
if (length(sig_genes) > 1) {
  pdf(file.path(OUTPUT_DIR, paste0("Heatmap_", COMP_PREFIX, ".pdf")), width = 7, height = 8)
  grid::grid.newpage()
  grid::grid.draw(p_heatmap$gtable)
  dev.off()
  
  png(file.path(OUTPUT_DIR, paste0("Heatmap_", COMP_PREFIX, ".png")), width = 7, height = 8, units = "in", res = 300)
  grid::grid.newpage()
  grid::grid.draw(p_heatmap$gtable)
  dev.off()
}

cat("[INFO] Pipeline finished successfully! Check files with prefix:", COMP_PREFIX, "\n")