library(qs)
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(stringr)
library(ggplot2)
library(openxlsx)

# MODULE : GO Enrichment and Barplot for genes -------------------

## 1 Parameters --------------------------------------------------------------
# -- 数据输入参数 --
input_filepath  <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/results/bulk_analysis_output/tables/DEG_300_1000.xlsx")
target_cluster  <- "GO"       # 目标提取的 cluster 名称

# -- 过滤与分析阈值 --
p_val_adj_cut   <- 0.05                   # 显著性阈值
log2FC_cut      <- 0.5                    # Log2FC 阈值
db_org          <- "org.Hs.eg.db"         # 物种数据库: "org.Mm.eg.db" (鼠) / "org.Hs.eg.db" (人)
ontology_type   <- "BP"                   # GO 类型: "BP", "CC", "MF" 或 "ALL"
top_n_go        <- 10                     # 图中展示的 Top GO 通路数量

# -- 绘图与美化参数 --
cluster_color   <- "#33A02C"              # 柱状图的填充颜色
wrap_width      <- 60                     # 文本自动换行的字符数限制

# -- 输出参数 --
out_dir         <- "results/bulk_analysis_output/plots"  # 图片保存的文件夹路径
out_width       <- 8                      # 输出图片宽度
out_height      <- 5                      # 输出图片高度


## 2 Process Data ------------------------------------------------------------
# 动态加载物种数据库
require(db_org, character.only = TRUE)

# 1. 读取 Marker 数据框
marker <- read.xlsx(input_filepath)
marker_df <- marker$gene_id

# 3. ID 转换 (SYMBOL -> ENTREZID)
gene_conv <- bitr(marker$gene_id, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = db_org)

# 4. 运行 GO 富集分析
ego_cluster <- enrichGO(
  gene          = gene_conv$ENTREZID,
  OrgDb         = db_org,
  ont           = ontology_type,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.5,
  readable      = TRUE
)

# 5. 提取富集结果并准备作图数据
plot_data <- as.data.frame(ego_cluster)

if(nrow(plot_data) > 0) {
  plot_data <- plot_data %>%
    dplyr::arrange(p.adjust) %>%
    head(top_n_go) %>%
    dplyr::mutate(
      logP = -log10(p.adjust),
      # 对过长的 Description 自动换行
      Description = stringr::str_wrap(Description, width = wrap_width) 
    )
  
  # 锁定因子顺序，确保显著性最高的排在最上面
  plot_data$Description <- factor(plot_data$Description, levels = rev(plot_data$Description))
} else {
  warning(paste("Warning: No significant GO terms found for", target_cluster))
}


## 3 Generate Plots ----------------------------------------------------------
# 初始化画图对象为空
p_barplot <- NULL

if(nrow(plot_data) > 0) {
  
  # 计算文字的 X 轴起始坐标 (放在坐标轴起点偏右一点点，确保整齐对齐)
  text_x_pos <- max(plot_data$logP) * 0.02 
  
  p_barplot <- ggplot(plot_data, aes(x = logP, y = Description)) +
    geom_col(width = 0.8, alpha = 0.8, fill = cluster_color, color = "black", linewidth = 0.3) +
    # 文字标签：固定在图表左侧对齐
    geom_text(aes(x = text_x_pos, label = Description), 
              size = 4.5, hjust = 0, color = "black") + 
    theme_minimal() +
    labs(
      title = paste0("Cluster ", target_cluster, " Top ", top_n_go, " GO Terms (", ontology_type, ")"), 
      x = "-log10(p.adjust)", 
      y = NULL
    ) +
    theme(
      axis.text.y  = element_blank(), # 隐藏原生 Y 轴文字
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(), # 隐藏网格线，画面更干净
      axis.line    = element_line(color = "black"),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    # 向右延展 X 轴空间，防止特别长的通路名字越界被裁切
    scale_x_continuous(expand = expansion(mult = c(0, 0.4))) 
}


## 4 Save Results ------------------------------------------------------------
if(!is.null(p_barplot)) {
  # 如果文件夹不存在，自动创建
  if(!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # 替换 cluster 名字里的非法字符（如果有），防止保存报错
  safe_name <- gsub("[^a-zA-Z0-9_]", "_", target_cluster) 
  file_base <- file.path(out_dir, paste0("Cluster_", safe_name, "_GO_", ontology_type))
  
  # 保存 PDF
  ggsave(filename = paste0(file_base, ".pdf"), 
         plot     = p_barplot, 
         width    = out_width, 
         height   = out_height)
  
  # 保存 PNG
  ggsave(filename = paste0(file_base, ".png"), 
         plot     = p_barplot, 
         width    = out_width, 
         height   = out_height, 
         dpi      = 300)
  
  print(paste("✅ GO barplot successfully generated and saved for:", target_cluster))
}