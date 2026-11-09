
# 步骤 0：清空环境与加载R包------------------

rm(list = ls()) # 清空环境变量，防止受之前的数据干扰

library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)


# 步骤 1：设置基础参数与创建文件夹------------------------

# 定义你的输入文件路径和输出文件夹名字
input_file <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/processed_data/All_gene_counts.list")
output_dir <- "bulk_analysis_output"

# 定义找差异基因的标准 (P值小于0.05，FoldChange大于2倍 即log2FC > 1)
p_adj_cutoff <- 0.05
log2fc_cutoff <- 1

# 在当前目录下新建几个文件夹存放结果，如果有了就不管它
dir.create(output_dir, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)  # 放图片
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE) # 放表格
dir.create(file.path(output_dir, "rdata"), showWarnings = FALSE)  # 放R数据


# 步骤 2：准备样本信息 (Metadata)---------------------

# 手动告诉软件，哪些样本是对照组(con)，哪些是处理组(treat)
sample_info <- data.frame(
  sample = c("con1", "con2", "con3", "con4", "treat1", "treat2", "treat3", "treat4"),
  condition = c("con", "con", "con", "con", "treat", "treat", "treat", "treat"),
  row.names = c("con1", "con2", "con3", "con4", "treat1", "treat2", "treat3", "treat4")
)

# 把condition变成因子，并把对照组(con)放在前面，这样算出来的Log2FC就是 treat 相比于 con 的变化
sample_info$condition <- factor(sample_info$condition, levels = c("con", "treat"))


# 步骤 3：读取并清理表达矩阵----------------------

# 读取基因表达量表格 (假设第一列是基因名，指定为行名 row.names = 1)
raw_df <- read.delim(input_file, sep = "\t", header = TRUE, row.names = 1)

# 把数据框变成矩阵格式，并且四舍五入变成整数（因为 DESeq2 只能吃整数的 count 数据）
count_mat <- round(as.matrix(raw_df))

# 确保表达矩阵的列名和我们的样本信息行名是对齐的
count_mat <- count_mat[, rownames(sample_info)]


# 步骤 4：过滤低表达基因--------------------------

# 规则：基因在至少2个样本中的表达量大于等于10，才保留它，否则删掉（节约计算时间）
keep <- rowSums(count_mat >= 10) >= 2
count_mat_filtered <- count_mat[keep, ]

# 打印一下过滤前后剩多少基因
print(paste("原始基因数:", nrow(count_mat), " | 过滤后保留基因数:", nrow(count_mat_filtered)))


# 步骤 5：使用 DESeq2 找差异基因 (核心步骤)

# 1. 把数据塞进 DESeq2 专用的对象里，告诉它我们要比较的是 condition 这一列
dds <- DESeqDataSetFromMatrix(countData = count_mat_filtered,
                              colData = sample_info,
                              design = ~ condition)

# 2. 一键执行 DESeq2 分析流程 (内部会自动进行标准化、离散度估计和统计检验)
dds <- DESeq(dds)

# 3. 提取对比结果：对比 "condition" 里 "treat" 相比于 "con" 的变化
res <- results(dds, contrast = c("condition", "treat", "con"), alpha = p_adj_cutoff)

# 把它转换成我们平时看的数据框(表格)格式，并把行名(基因名)单独拉出来做一列
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)


# 步骤 6：给基因打标签 (上调、下调、不显著)

# 先全部默认打上 "NotSig" (不显著) 的标签
res_df$regulation <- "NotSig"

# 找到 padj < 0.05 并且 log2FC >= 1 的，改成 "Up" (上调)
res_df$regulation[!is.na(res_df$padj) & res_df$padj < p_adj_cutoff & res_df$log2FoldChange >= log2fc_cutoff] <- "Up"

# 找到 padj < 0.05 并且 log2FC <= -1 的，改成 "Down" (下调)
res_df$regulation[!is.na(res_df$padj) & res_df$padj < p_adj_cutoff & res_df$log2FoldChange <= -log2fc_cutoff] <- "Down"

# 提取出只包含显著上调和下调基因的子表格，方便后续保存
deg_up <- res_df[res_df$regulation == "Up", ]
deg_down <- res_df[res_df$regulation == "Down", ]


# 步骤 7：数据转换 (为了后面画图好看)

# 提取 DESeq2 标准化后的表达矩阵
norm_counts <- counts(dds, normalized = TRUE)

# 对数据进行 VST 转换 (方差稳定变换)，用来画PCA和热图效果最好
vst_obj <- vst(dds, blind = FALSE)
vst_mat <- assay(vst_obj) 


# 步骤 8：画图并保存到 plots 文件夹

plot_dir <- file.path(output_dir, "plots")

# 画图1：PCA图 (看组内样本聚没聚在一起，组间分没分开)
pca_data <- plotPCA(vst_obj, intgroup = "condition", returnData = TRUE)
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4) + geom_text(vjust = -1) + theme_bw() + labs(title = "PCA Plot")
ggsave(file.path(plot_dir, "01_PCA_plot.pdf"), p_pca, width = 7, height = 6)

# 画图2：样本相关性热图
sample_dists <- dist(t(vst_mat))
pdf(file.path(plot_dir, "02_sample_distance.pdf"), width = 7, height = 6)
pheatmap(as.matrix(sample_dists), main = "Sample Distance Heatmap")
dev.off()

# 画图3：火山图 (一览所有基因的变化趋势)
p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "grey")) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(p_adj_cutoff), linetype = "dashed") +
  theme_bw() + labs(title = "Volcano Plot: treat vs con")
ggsave(file.path(plot_dir, "03_volcano_plot.pdf"), p_volcano, width = 7, height = 6)


# 步骤 9：保存表格和数据到 tables 文件夹

table_dir <- file.path(output_dir, "tables")

# 保存各种表格 (原始Count、标准化Count、所有结果、上调基因、下调基因)
write.csv(count_mat, file.path(table_dir, "raw_counts.csv"))
write.csv(norm_counts, file.path(table_dir, "normalized_counts.csv"))
write.csv(res_df, file.path(table_dir, "DEG_all_results.csv"), row.names = FALSE)
write.csv(deg_up, file.path(table_dir, "DEG_up.csv"), row.names = FALSE)
write.csv(deg_down, file.path(table_dir, "DEG_down.csv"), row.names = FALSE)

# 保存整个 R 工作空间 (下次直接 load 这个文件，所有变量都在，不用重新跑)
save.image(file = file.path(output_dir, "rdata", "bulk_analysis_workspace.RData"))

print("分析完成！去看看 bulk_analysis_output 文件夹吧！")