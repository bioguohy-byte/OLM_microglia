library(Seurat)
library(ggplot2)
library(dplyr)

obj_with_M0_threshold_labels_v5 <- readRDS(Sys.getenv("MG_EXTERNAL_029"))

obj_with_M0_threshold_labels_v5$Celltype <- obj_with_M0_threshold_labels_v5$SCT_snn_res.0.8
obj_with_M0_threshold_labels_v5@meta.data <- obj_with_M0_threshold_labels_v5@meta.data %>%
  mutate(Celltype = case_when(
    Celltype %in% c("13") ~ "MG0-like",
    Celltype %in% c("1", "2","3","4","5","6","7","8","9","10","11","12") ~ "Other MG",
    TRUE ~ Celltype
  ))

pdf(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/results/figure/spinal_umap_v5.pdf"),width = 5,height = 4)
DimPlot(obj_with_M0_threshold_labels_v5,group = 'Celltype')
dev.off()





#!/usr/bin/env Rscript

# MODULE : FeatureStatPlot Generation ----------------------------------------
# Description: 基于 Seurat 对象，动态生成带有显著性标记的基因表达小提琴图，
#              支持自定义分组顺序、指定对照组（Reference）并自动配对比较。

rm(list = ls())
gc()

# 载入依赖包
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(SCP)
})

## 1 Parameters --------------------------------------------------------------

# --- I/O Paths (输入与输出路径) ---
IN_RDS_PATH  <- Sys.getenv("MG_EXTERNAL_029")

# 输出路径
OUT_PDF_PATH <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/results/figure/Violin_Stat_Markers_v5.pdf")

# --- Grouping & Comparisons (细胞分组与统计比较策略) ---
GROUP_COL    <- "Celltype"
TARGET_ORDER <- c("MG0-like", "Other MG")
REF_GROUP    <- "MG0-like"

# --- Plotting Variables (绘图特征与视觉参数) ---
TARGET_GENES <- c("TRIM2", "RNF220", "MBP", "MOBP", "MOG", "OLIG1", "BCAS1", "ST18")

custom_colors <- c(
  "MG0-like" = "#F8766D",
  "Other MG" = "#00BFC4"
)
PALETTE_COLORS <- custom_colors

# 图片保存尺寸
PLOT_WIDTH  <- 12
PLOT_HEIGHT <- 6


## Helper -------------------------------------------------------------------

# 为 SCP 构造一个兼容 Seurat v4 Assay 的绘图 assay，避免 Assay5 -> Assay 直接转换报错
make_scp_compatible_assay <- function(seu, src_assay = "RNA", dst_assay = "RNA_SCP") {
  assay_obj <- seu[[src_assay]]
  
  if (!inherits(assay_obj, "Assay5")) {
    message("Assay '", src_assay, "' is already compatible.")
    DefaultAssay(seu) <- src_assay
    return(list(seu = seu, assay = src_assay))
  }
  
  message("Detected Assay5 in '", src_assay, "'. Building SCP-compatible assay...")
  
  layer_names <- Layers(assay_obj)
  message("Layers before processing: ", paste(layer_names, collapse = ", "))
  
  # 如果是 split layers（如 counts.xxx / data.xxx），先合并
  if (any(grepl("^(counts|data)\\.", layer_names))) {
    seu <- suppressWarnings(JoinLayers(seu, assay = src_assay))
    layer_names <- Layers(seu[[src_assay]])
    message("Layers after JoinLayers: ", paste(layer_names, collapse = ", "))
  }
  
  pick_layer <- function(prefix) {
    current_layers <- Layers(seu[[src_assay]])
    candidates <- unique(c(
      prefix,
      grep(paste0("^", gsub("\\.", "\\\\.", prefix)), current_layers, value = TRUE)
    ))
    
    for (lay in candidates) {
      mat <- suppressWarnings(
        tryCatch(
          GetAssayData(seu, assay = src_assay, layer = lay),
          error = function(e) NULL
        )
      )
      if (!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
        return(mat)
      }
    }
    return(NULL)
  }
  
  counts_mat <- pick_layer("counts")
  data_mat   <- pick_layer("data")
  scale_mat  <- pick_layer("scale.data")
  
  if (is.null(counts_mat) && is.null(data_mat)) {
    stop("RNA assay 中没有可用的 counts/data layer，无法生成兼容 SCP 的 assay。")
  }
  
  # CreateAssayObject 必须先给 counts；如果没有 counts，就先用 data 顶上用于作图兼容
  seed_mat <- if (!is.null(counts_mat)) counts_mat else data_mat
  compat_assay <- CreateAssayObject(counts = seed_mat)
  
  if (!is.null(data_mat)) {
    compat_assay <- SetAssayData(
      object = compat_assay,
      slot = "data",
      new.data = data_mat
    )
  } else {
    compat_assay <- SetAssayData(
      object = compat_assay,
      slot = "data",
      new.data = seed_mat
    )
    message("No 'data' layer found. Using counts as data for plotting.")
  }
  
  if (!is.null(scale_mat) && nrow(scale_mat) > 0 && ncol(scale_mat) > 0) {
    compat_assay <- SetAssayData(
      object = compat_assay,
      slot = "scale.data",
      new.data = scale_mat
    )
  }
  
  seu[[dst_assay]] <- compat_assay
  DefaultAssay(seu) <- dst_assay
  
  message("Created compatible assay: ", dst_assay)
  return(list(seu = seu, assay = dst_assay))
}


## 2 Process Data ------------------------------------------------------------
message("--- Step 2: Loading and Processing Data ---")

# 1. 读取 Seurat 对象
seu <- readRDS(IN_RDS_PATH)

# 2. 为 SCP 构建兼容 assay
compat_res <- make_scp_compatible_assay(
  seu = seu,
  src_assay = "RNA",
  dst_assay = "RNA_SCP"
)
seu <- compat_res$seu
DefaultAssay(seu) <- compat_res$assay

# 3. 整合新的 Celltype 注释逻辑
seu$Celltype <- seu$SCT_snn_res.0.8
seu@meta.data <- seu@meta.data %>%
  mutate(Celltype = case_when(
    Celltype %in% c("13") ~ "MG0-like",
    Celltype %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12") ~ "Other MG",
    TRUE ~ as.character(Celltype)
  ))

# 4. 安全设置分组顺序，并去掉空 level
group_vec <- as.character(seu[[GROUP_COL]][, 1])
present_groups <- TARGET_ORDER[TARGET_ORDER %in% unique(group_vec)]
seu[[GROUP_COL]] <- factor(group_vec, levels = present_groups)

message("Current cluster levels: ", paste(levels(seu[[GROUP_COL]][, 1]), collapse = ", "))

# 5. 动态生成比较组合
all_clusters <- levels(seu[[GROUP_COL]][, 1])
compare_groups <- setdiff(all_clusters, REF_GROUP)

if (!(REF_GROUP %in% all_clusters)) {
  stop("Reference group '", REF_GROUP, "' 不存在于当前分组中。")
}
if (length(compare_groups) == 0) {
  stop("没有可用于与 Reference group 比较的其它分组。")
}

my_comparisons <- lapply(compare_groups, function(x) c(x, REF_GROUP))
message("Comparisons generated against Reference: ", REF_GROUP)


## 3 Generate Plots ----------------------------------------------------------
message("--- Step 3: Generating Violin Plots ---")

p_violin <- FeatureStatPlot(
  srt         = seu,
  group.by    = GROUP_COL,
  stat.by     = TARGET_GENES,
  add_box     = FALSE,
  comparisons = my_comparisons,
  bg.by       = GROUP_COL,
  palcolor    = PALETTE_COLORS,
  bg_palcolor = PALETTE_COLORS,
  ncol        = 4
)


## 4 Save Results ------------------------------------------------------------
message("--- Step 4: Saving Plot Results ---")

out_dir <- dirname(OUT_PDF_PATH)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(
  filename = OUT_PDF_PATH,
  plot     = p_violin,
  width    = PLOT_WIDTH,
  height   = PLOT_HEIGHT
)

message("=====================================================================")
message("SUCCESS: Plot saved to -> ", OUT_PDF_PATH)
message("=====================================================================")
