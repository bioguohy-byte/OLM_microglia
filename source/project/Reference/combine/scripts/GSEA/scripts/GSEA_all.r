setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/"))
# MODULE 10: GSEA Analysis (All GO Terms) ------------------------------------

## 10.1 Parameters -------------------------------------------------------------
input_file      <- 'data/processed_data/mg_cca_meta_v2.qs'
dir_data        <- "scripts/GSEA/data/"
dir_res         <- "scripts/GSEA/results/"
custom_set_file <- paste0(dir_data, 'my_gene_set.xlsx')
go_gmt_file     <- paste0(dir_data, "c5.go.v2023.2.Hs.symbols.gmt")

target_group_by <- 'CellType1'
target_cluster  <- 'MG0'  # Modify this to change the cluster

output_qs       <- paste0(dir_data, 'GSEA_kk_ALL_', target_cluster, '.qs')
# 导出的完整表格路径
output_txt      <- paste0(dir_res, target_cluster, "_ALL_GO_GSEA_result.txt") 
output_csv      <- paste0(dir_res, target_cluster, "_ALL_GO_GSEA_result.csv") # 额外存一个CSV方便Excel查看

dir.create(dir_res, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_data, recursive = TRUE, showWarnings = FALSE)

## 10.2 Process Data -----------------------------------------------------------
library(Seurat)
library(clusterProfiler)
library(openxlsx)
library(dplyr)
library(qs)

if(!file.exists(input_file)) stop("Input file not found.")
tmp_obj <- qread(input_file)

# 1. 准备自定义 Gene Sets
gene_sets_raw <- read.xlsx(custom_set_file, colNames = FALSE, na.strings = "")
temp_gmt_path <- paste0(dir_data, 'temp_my_genesets.gmt')
write.table(gene_sets_raw, temp_gmt_path, sep = "\t", row.names = F, col.names = F, quote = F)

custom_gmt <- read.gmt(temp_gmt_path)
custom_gmt$gene <- toupper(custom_gmt$gene)

# 2. 准备完整的 GO Gene Sets (不再进行 target_go_terms 过滤)
gmt_go <- read.gmt(go_gmt_file)
processed_go <- gmt_go %>%
  # 使用正则匹配去掉前缀 (例如 GOBP_, GOCC_, GOMF_)
  mutate(term = sub("^GO[A-Z]*_", "", term)) %>% 
  # 将下划线替换为空格，并转为小写，方便阅读
  mutate(term = tolower(gsub("_", " ", term)))

# 将自定义基因集和全部GO基因集合并
combined_gmt <- rbind(custom_gmt, processed_go)

# 3. 准备 Ranked Gene List
Idents(tmp_obj) <- target_group_by

deg_genes <- FindMarkers(
  tmp_obj, 
  ident.1 = target_cluster, 
  min.pct = 0.1, 
  logfc.threshold = 0,
  only.pos = FALSE
)

deg_genes$gene <- rownames(deg_genes)
rt_sorted <- deg_genes %>%
  arrange(desc(avg_log2FC)) %>%
  distinct(gene, .keep_all = TRUE)

logFC_vector <- rt_sorted$avg_log2FC
names(logFC_vector) <- toupper(rt_sorted$gene)

# 4. 运行全集 GSEA
kk <- GSEA(
  geneList      = logFC_vector, 
  TERM2GENE     = combined_gmt, 
  pvalueCutoff  = 1,     # 保持为1，这样会输出所有通路，方便你自己在表格里用Excel筛选
  pAdjustMethod = 'BH',
  minGSSize     = 10,    # 建议设为10或15，过滤掉包含基因太少没有统计意义的通路
  maxGSSize     = 500,
  scoreType     = "pos", # 注意：如果你的排序包含负值，且你想看下调通路，这里可能需要去掉或改为"std"
  eps           = 0,
  verbose       = TRUE
)

## 10.3 Save Results -----------------------------------------------------------
# 保存对象以便后续画图直接读取，不用重新跑GSEA
qsave(kk, file = output_qs)

# 导出为 txt 和 csv 格式，方便用 Excel 打开慢慢筛选
gsea_result_df <- as.data.frame(kk)
write.table(gsea_result_df, file = output_txt, sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(gsea_result_df, file = output_csv, row.names = FALSE)

if(file.exists(temp_gmt_path)) unlink(temp_gmt_path)

message(paste0("Module 10 (Full GSEA) Done. Please check the result table for cluster: ", target_cluster))