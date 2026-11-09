#!/usr/bin/env Rscript

# Audit-only CellType evidence package for all non-Yao human datasets.
# It never overwrites source objects and never assigns a CellType automatically.

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(harmony)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: run_yao_marker_celltype_audit.R <manifest.tsv> <output_dir>")
manifest <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
out_dir <- args[[2]]
dir.create(file.path(out_dir, "evidence"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "checkpoints"), recursive = TRUE, showWarnings = FALSE)
set.seed(42)

# Exact marker panel in Yao_et_al./scripts/00.yao_alldata.r, section 10.
yao_markers <- data.frame(
  Gname = c("AQP4","ITGB4","ETNPPL","SLC1A2","RFX4","GLI3","GRAMD2B","GFAP",
            "VCAN","C1QL1","SUSD5","COL9A1",
            "MOG","MOBP","SLC5A11","PICK1","HAPLN2","CNDP1","MYRF","CPM","TF","TMC6","CDK18",
            "PTPRC",
            "CSF1R","MRC1","IL18","C3","SLC2A5","BIN2","DOCK2","ENTPD1",
            "CFAP54","DNAH11",
            "DCN","LAMA2","RERGL","ABCG2","EMCN","LEF1","FLT1","SLC6A12","CFH","SLCO1B3",
            "RBFOX3","SNAP25","SYT1","GRIN1","CD247"),
  Type = c(rep("Astr",8),rep("OPC",4),rep("Oli",11),rep("I",1),rep("Micro",8),
           rep("Ependymal",2),rep("Vascular",10),rep("Neurons",4),rep("T cells",1)),
  stringsAsFactors = FALSE
)

numeric_cluster_order <- function(x) {
  x <- as.character(x); as_num <- suppressWarnings(as.numeric(x))
  if (all(!is.na(as_num))) as.character(sort(unique(as_num))) else sort(unique(x))
}

ensure_normalized_rna <- function(obj) {
  if (!"RNA" %in% Assays(obj)) stop("Object has no RNA assay for Yao marker DotPlot.")
  DefaultAssay(obj) <- "RNA"
  data_ok <- tryCatch({
    d <- GetAssayData(obj, assay = "RNA", layer = "data")
    nrow(d) > 0L && ncol(d) > 0L
  }, error = function(e) FALSE)
  if (!data_ok) obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
  obj
}

match_markers <- function(obj, table) {
  features <- rownames(obj[[DefaultAssay(obj)]])
  lookup <- setNames(features, toupper(features))
  table$feature <- unname(lookup[toupper(table$Gname)])
  table$present <- !is.na(table$feature)
  table
}

# Faithful to the original Yao marker-DotPlot logic/palette. The only added
# guard is to drop marker genes absent from a particular study before plotting.
make_yao_dotplot <- function(obj, cluster_col, markers, label_size = 10) {
  markers <- markers[markers$present, , drop = FALSE]
  if (!nrow(markers)) stop("None of the Yao marker genes is present.")
  all_dot <- DotPlot(obj, features = unique(markers$feature), group.by = cluster_col, scale = TRUE)
  lim <- range(all_dot$data$avg.exp.scaled, na.rm = TRUE)
  if (!all(is.finite(lim))) lim <- c(-1, 1)
  if (diff(lim) == 0) lim <- lim + c(-1, 1)
  plots <- list(); widths <- numeric()
  for (type_name in unique(markers$Type)) {
    cur <- markers[markers$Type == type_name, , drop = FALSE]
    dot <- DotPlot(obj, features = unique(cur$feature), group.by = cluster_col, scale = TRUE)$data
    dot <- merge(dot, cur[, c("feature", "Gname", "Type")], by.x = "features.plot", by.y = "feature", all.x = TRUE, sort = FALSE)
    dot$Gname <- factor(dot$Gname, levels = cur$Gname)
    p <- ggplot(dot, aes(x = Gname, y = id, size = pct.exp, colour = avg.exp.scaled)) +
      geom_point() + scale_size_continuous(range = c(1, 4), limits = c(0, 100)) +
      scale_color_gradientn(colors = rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", "#ffe6b7", "grey90", "white")), limits = lim) +
      facet_wrap(~Type) +
      guides(colour = guide_colourbar(title = "Average Expression"), size = guide_legend(title = "Percent Expressed")) +
      theme_classic() + theme(panel.grid = element_blank(), panel.background = element_blank(), axis.title = element_blank(),
                              axis.text = element_text(size = label_size, colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
                              legend.text = element_text(size = label_size), legend.title = element_text(size = label_size),
                              strip.background = element_rect(color = "grey90", fill = "grey90"), strip.text = element_text(size = label_size, colour = "black"))
    if (length(plots)) p <- p + theme(axis.title.y = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    plots[[length(plots) + 1L]] <- p; widths <- c(widths, length(unique(cur$feature)))
  }
  Reduce(`+`, plots) + plot_layout(ncol = length(plots), widths = widths, guides = "collect")
}

# Su has no original stored UMAP/cluster. The original Yao script contains an
# explicit Harmony integration branch (section 9.3). It is preferable here:
# 16 batches would require 120 all-pair RPCA comparisons.
run_su_harmony_audit <- function(source_obj, batch_col, checkpoint_path) {
  if (!batch_col %in% colnames(source_obj@meta.data)) stop("Su batch column absent: ", batch_col)
  message("Su: new audit-only Harmony clustering following Yao Harmony parameters.")
  source_obj <- ensure_normalized_rna(source_obj)
  source_obj <- FindVariableFeatures(source_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  # PCA/Harmony only consume the selected 2,000 variable features. Limiting
  # scaling to that same set avoids a 22,000+ gene dense matrix without
  # changing the features used for Su's audit clustering.
  source_obj <- ScaleData(source_obj, features = VariableFeatures(source_obj), verbose = FALSE)
  source_obj <- RunPCA(source_obj, npcs = 30, verbose = FALSE)
  source_obj <- RunHarmony(source_obj, group.by.vars = batch_col, verbose = FALSE)
  source_obj <- RunUMAP(source_obj, seed.use = 42, reduction = "harmony", dims = 1:30, verbose = FALSE)
  source_obj <- FindNeighbors(source_obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
  # Yao's final CellType mapping reads the resolution-0.5 partition. Compute
  # that exact partition only; the other 19 exploratory resolutions in the
  # original script are not used by this audit.
  source_obj <- FindClusters(source_obj, resolution = 0.5, verbose = FALSE)
  res_col <- grep("_snn_res\\.0\\.5$", colnames(source_obj@meta.data), value = TRUE)
  if (length(res_col) != 1L) stop("Could not uniquely identify the Harmony resolution-0.5 column.")
  source_obj$audit_harmony_snn_res.0.5 <- source_obj@meta.data[[res_col]]
  qsave(source_obj, checkpoint_path)
  source_obj
}

cluster_rows <- list(); marker_rows <- list(); provenance_rows <- list()
for (i in seq_len(nrow(manifest))) {
  r <- manifest[i, , drop = FALSE]; study <- r$study[[1]]
  message("\n===== ", study, " =====")
  obj <- qread(r$object_path[[1]]); cluster_col <- r$cluster_col[[1]]; reduction <- r$reduction[[1]]
  if (study == "Su") {
    checkpoint <- file.path(out_dir, "checkpoints", "Su_audit_harmony_clustered.qs")
    if (file.exists(checkpoint)) { rm(obj); gc(); obj <- qread(checkpoint) } else obj <- run_su_harmony_audit(obj, r$batch_col[[1]], checkpoint)
  }
  if (!cluster_col %in% colnames(obj@meta.data)) stop(study, ": missing cluster column ", cluster_col)
  if (!reduction %in% Reductions(obj)) stop(study, ": missing reduction ", reduction)
  levels_use <- numeric_cluster_order(obj@meta.data[[cluster_col]])
  obj@meta.data[[cluster_col]] <- factor(as.character(obj@meta.data[[cluster_col]]), levels = levels_use)
  p_umap <- DimPlot(obj, reduction = reduction, group.by = cluster_col, label = TRUE, repel = TRUE, pt.size = 0.1, raster = TRUE) +
    ggtitle(paste0(study, ": audit clusters (no CellType assignment)")) + theme_classic() + theme(legend.position = "none", plot.title = element_text(face = "bold"))
  obj <- ensure_normalized_rna(obj); markers <- match_markers(obj, yao_markers)
  p_dot <- make_yao_dotplot(obj, cluster_col, markers, label_size = 10)
  # The original Yao panel used a 20 x 7 in canvas. Here the same 49-gene
  # panel is placed beside an UMAP, so a wider canvas prevents truncation of
  # marker-family facet titles such as Ependymal and T cells.
  p_all <- p_umap + p_dot + plot_layout(widths = c(1, 3.5)) + plot_annotation(title = paste0(study, ": Yao-marker CellType audit — reviewer assignment required"))
  ggsave(file.path(out_dir, "evidence", paste0(study, "_Yao_marker_cluster_audit.pdf")), p_all, width = 28, height = 8, limitsize = FALSE)
  ggsave(file.path(out_dir, "evidence", paste0(study, "_Yao_marker_cluster_audit.png")), p_all, width = 28, height = 8, dpi = 220, limitsize = FALSE)
  label_col <- r$existing_label_col[[1]]
  if (!is.na(label_col) && label_col %in% colnames(obj@meta.data)) {
    p_context <- DimPlot(obj, reduction = reduction, group.by = label_col, label = TRUE, repel = TRUE, pt.size = 0.1, raster = TRUE) +
      ggtitle(paste0(study, ": existing metadata label (context only; not adopted)")) + theme_classic()
    ggsave(file.path(out_dir, "evidence", paste0(study, "_existing_metadata_label_context.pdf")), p_context, width = 11, height = 8.5, limitsize = FALSE)
  }
  counts <- as.data.frame(table(obj@meta.data[[cluster_col]]), stringsAsFactors = FALSE); colnames(counts) <- c("audit_cluster", "n_cells")
  counts$study <- study; counts$cluster_source <- r$cluster_source[[1]]; counts$reviewer_final_CellType <- ""; counts$reviewer_notes <- ""; counts$review_status <- "pending"
  cluster_rows[[length(cluster_rows) + 1L]] <- counts[, c("study","audit_cluster","n_cells","cluster_source","reviewer_final_CellType","reviewer_notes","review_status")]
  markers$study <- study; marker_rows[[length(marker_rows) + 1L]] <- markers[, c("study","Type","Gname","feature","present")]
  provenance_rows[[length(provenance_rows) + 1L]] <- data.frame(study = study, object_path = r$object_path[[1]], cluster_column = cluster_col, reduction = reduction, batch_column = r$batch_col[[1]], existing_label_column = label_col, cluster_source = r$cluster_source[[1]], stringsAsFactors = FALSE)
  rm(obj); gc()
}
write.csv(do.call(rbind, cluster_rows), file.path(out_dir, "tables", "reviewer_cluster_assignment_template.csv"), row.names = FALSE)
write.csv(do.call(rbind, marker_rows), file.path(out_dir, "tables", "yao_marker_availability_by_study.csv"), row.names = FALSE)
write.csv(do.call(rbind, provenance_rows), file.path(out_dir, "tables", "cluster_and_object_provenance.csv"), row.names = FALSE)
writeLines(c("# Yao-marker CellType audit package", "", "Original .qs objects were opened read-only and not modified.", "Every *_Yao_marker_cluster_audit.pdf pairs a cluster UMAP with the exact marker panel in Yao_et_al./scripts/00.yao_alldata.r, section 10.", "No CellType labels have been generated. reviewer_cluster_assignment_template.csv has intentionally blank reviewer fields.", "Existing metadata labels are shown only in separately named context PDFs and never populate a reviewer decision.", "Su had no stored UMAP/cluster. Its separate audit-only Harmony checkpoint uses Batch, 2,000 variable genes, PCA 1:30, Harmony batch correction, seed 42 and resolution 0.5 following the Yao Harmony branch."), file.path(out_dir, "README.md"))
