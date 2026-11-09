# ==============================================================================
# 0. 加载必要的包
# ==============================================================================
library(Seurat)
library(data.table)
library(qs)

# 定义文件路径 (方便后续修改)
path_expr <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Fullard_et_al./data/raw_data/GSE164485_gene_expression_20210104.csv.gz")
path_meta <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Fullard_et_al./data/raw_data/GSE164485_meta_data_20210104.csv.gz")

# ==============================================================================
# 1. 读取并处理表达矩阵 (Expression Matrix)
# ==============================================================================
message("正在读取表达矩阵...")
# 使用 fread 快速读取
expr_data <- fread(path_expr)
expr_data <- as.data.frame(expr_data)

# 设置行名 (假设第一列是基因名/或细胞名，取决于原始文件)
rownames(expr_data) <- expr_data[, 1]
expr_data <- expr_data[, -1] 

# 【关键步骤：转置矩阵】
# 原始数据是 "行=细胞，列=基因"，Seurat 需要 "行=基因，列=细胞"
message("正在转置矩阵 (这可能需要一点时间)...")
expr_data <- t(expr_data)

# 检查维度 (行数应为基因数，列数应为细胞数 ~68k)
print(dim(expr_data))

# ==============================================================================
# 2. 读取元数据并创建 Seurat 对象
# ==============================================================================
message("正在读取元数据...")
meta_data <- read.csv(path_meta, row.names = 1)

message("正在创建 Seurat 对象...")
seurat_obj <- CreateSeuratObject(counts = expr_data, 
                                 meta.data = meta_data, 
                                 project = "GSE164485")

# 【内存优化】创建对象后，原始的巨大矩阵就不需要了，删除以释放内存
rm(expr_data)
gc() # 强制进行垃圾回收

# ==============================================================================
# 3. 导入作者原本的 UMAP 坐标
# ==============================================================================
message("正在导入 UMAP 坐标...")

# 从 meta.data 提取坐标 (确保列名 UMAP1/UMAP2 与 CSV 中完全一致)
# 如果报错，请检查 colnames(seurat_obj@meta.data) 确认列名
umap_coords <- seurat_obj@meta.data[, c("UMAP1", "UMAP2")]
umap_matrix <- as.matrix(umap_coords)

# 创建 DimReduc 对象并添加到 Seurat 对象中
seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = umap_matrix, 
                                             key = "UMAP_", 
                                             assay = DefaultAssay(seurat_obj))

# 验证对象结构
print(seurat_obj)

# ==============================================================================
# 4. 数据标准化 (Normalization)
# ==============================================================================
# FeaturePlot 需要标准化后的数据 (data slot)，否则表达量展示不准确
message("正在进行数据标准化...")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
qsave(seurat_obj,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Fullard_et_al./data/processed_data/Fullard.qs"))
# tmp_obj <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Fullard_et_al./data/processed_data/Fullard.qs"))

# ==============================================================================
# 5. 可视化作图
# ==============================================================================
message("开始绘图...")

# 5.1 聚类图 (DimPlot)
# group.by 参数可以是 metadata 中的任意一列，例如 "cell_type", "Cluster", "orig.ident"
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") + ggtitle("Sample Identity")

# 5.2 基因表达图 (FeaturePlot)
# 这里的 TRIM2 是你想查看的基因
p2 <- FeaturePlot(seurat_obj, features = "TRIM2", reduction = "umap")

# 显示图片
print(p1)
print(p2)

message("分析流程完成！")