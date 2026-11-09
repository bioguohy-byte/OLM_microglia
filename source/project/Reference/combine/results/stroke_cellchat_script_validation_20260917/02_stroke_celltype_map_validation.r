# Validation-only copy of:
#   combine/scripts/CellChat/scripts/stroke_celltype_map.r
# The mapping logic is unchanged; only input/output paths point to validation files.

library(Seurat)
library(qs)

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
validation_root <- file.path(project_root, "combine/results/stroke_cellchat_script_validation_20260917")

seurat_obj <- qread(file.path(validation_root, "Yao_stroke_reproduced.qs"))
mg_meta <- qread(file.path(project_root, "combine/data/processed_data/mg_cca_meta_v2.qs"))

transfer_sub_labels <- function(seu_main,
                                seu_sub,
                                main_col = "CellType1",
                                sub_col = "CellType2",
                                new_col = "cellchat_type") {
  if (!main_col %in% colnames(seu_main@meta.data)) {
    stop(paste("主对象中找不到列名:", main_col))
  }
  if (!sub_col %in% colnames(seu_sub@meta.data)) {
    stop(paste("子对象中找不到列名:", sub_col))
  }
  seu_main[[new_col]] <- as.character(seu_main[[main_col, drop = TRUE]])
  mapping_vec <- setNames(as.character(seu_sub[[sub_col, drop = TRUE]]), colnames(seu_sub))
  common_cells <- intersect(names(mapping_vec), colnames(seu_main))
  if (length(common_cells) > 0) {
    seu_main@meta.data[common_cells, new_col] <- mapping_vec[common_cells]
    message(sprintf("已更新 %d 个共有细胞的标签。", length(common_cells)))
  } else {
    warning("两个对象之间没有检测到共有细胞 ID。")
  }
  seu_main
}

# Same call and label order as the original script.
seurat_obj <- transfer_sub_labels(
  seurat_obj,
  mg_meta,
  main_col = "combine_celltype",
  sub_col = "CellType1",
  new_col = "cellchat_type"
)

Idents(seurat_obj) <- "cellchat_type"
desired_order <- c("MG0", "MG1", "MG2", "MG3", "MG4", "MG5", "Neuron",
                   "OLIG", "OPC", "Astrocyte", "Vascular", "T cell", "Ependymal")
seurat_obj$cellchat_type <- factor(seurat_obj$cellchat_type, levels = desired_order)

output_path <- file.path(validation_root, "stroke_cellchat_reproduced.qs")
qsave(seurat_obj, output_path)
cat("Reproduced CellChat input object:\n", output_path, "\n", sep = "")
print(table(seurat_obj$cellchat_type, useNA = "ifany"))
