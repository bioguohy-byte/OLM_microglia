#!/usr/bin/env Rscript

# Corrected, read-only Chia CellType evidence package.
# The source object is never saved or modified; reviewer labels are deliberately blank.

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: build_chia_corrected_celltype_evidence.R <Chia.qs> <out_dir>")
source_path <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
set.seed(42)

# Exact Yao marker list, Yao_et_al./scripts/00.yao_alldata.r section 10.
marker_tbl <- data.frame(
  Gname = c("AQP4","ITGB4","ETNPPL","SLC1A2","RFX4","GLI3","GRAMD2B","GFAP",
            "VCAN","C1QL1","SUSD5","COL9A1",
            "MOG","MOBP","SLC5A11","PICK1","HAPLN2","CNDP1","MYRF","CPM","TF","TMC6","CDK18",
            "PTPRC","CSF1R","MRC1","IL18","C3","SLC2A5","BIN2","DOCK2","ENTPD1",
            "CFAP54","DNAH11","DCN","LAMA2","RERGL","ABCG2","EMCN","LEF1","FLT1","SLC6A12","CFH","SLCO1B3",
            "RBFOX3","SNAP25","SYT1","GRIN1","CD247"),
  Type = c(rep("Astrocyte",8), rep("OPC",4), rep("Oligodendrocyte",11), rep("Immune",1),
           rep("Microglia",8), rep("Ependymal",2), rep("Vascular",10), rep("Neurons",4), rep("T cells",1)),
  stringsAsFactors = FALSE
)

obj <- qread(source_path)
cluster_col <- "seurat_clusters"
label_col <- "cell.type"

if (!cluster_col %in% colnames(obj@meta.data)) stop("Missing final cluster column: ", cluster_col)
if (!label_col %in% colnames(obj@meta.data)) stop("Missing historical label column: ", label_col)
if (!"umap" %in% Reductions(obj)) stop("Missing stored UMAP.")
if (!identical(rownames(obj@meta.data), colnames(obj))) stop("Metadata barcodes do not align with the object.")
if (!identical(rownames(Embeddings(obj, "umap")), colnames(obj))) stop("UMAP barcodes do not align with the object.")
if (!identical(as.character(Idents(obj)), as.character(obj@meta.data[[cluster_col]]))) {
  stop("Active identities are not the final seurat_clusters field.")
}

cluster_levels <- as.character(sort(as.integer(unique(as.character(obj@meta.data[[cluster_col]])))))
obj@meta.data[[cluster_col]] <- factor(as.character(obj@meta.data[[cluster_col]]), levels = cluster_levels)

DefaultAssay(obj) <- "RNA"
rna_data <- GetAssayData(obj, assay = "RNA", layer = "data")
feature_lookup <- setNames(rownames(rna_data), toupper(rownames(rna_data)))
marker_tbl$feature <- unname(feature_lookup[toupper(marker_tbl$Gname)])
marker_tbl$present <- !is.na(marker_tbl$feature)
if (!all(marker_tbl$present)) stop("Not all Yao marker genes are present in the Chia RNA assay.")

make_dotplot <- function(object, markers, cluster_column) {
  global_dot <- DotPlot(object, features = unique(markers$feature), group.by = cluster_column, scale = TRUE)$data
  lim <- range(global_dot$avg.exp.scaled, na.rm = TRUE)
  if (!all(is.finite(lim)) || diff(lim) == 0) lim <- c(-1, 1)
  panels <- list(); widths <- numeric()
  for (type_name in unique(markers$Type)) {
    cur <- markers[markers$Type == type_name, , drop = FALSE]
    dot <- DotPlot(object, features = unique(cur$feature), group.by = cluster_column, scale = TRUE)$data
    dot <- merge(dot, cur[, c("feature", "Gname")], by.x = "features.plot", by.y = "feature", all.x = TRUE, sort = FALSE)
    dot$Gname <- factor(dot$Gname, levels = cur$Gname)
    dot$id <- factor(as.character(dot$id), levels = cluster_levels)
    dot$Type <- type_name
    p <- ggplot(dot, aes(x = Gname, y = id, size = pct.exp, colour = avg.exp.scaled)) +
      geom_point() +
      scale_size_continuous(range = c(0.7, 3.6), limits = c(0, 100)) +
      scale_colour_gradientn(colors = rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", "#ffe6b7", "grey90", "white")), limits = lim) +
      facet_wrap(~Type) +
      guides(colour = guide_colourbar(title = "Average expression"), size = guide_legend(title = "Percent expressed")) +
      theme_classic(base_size = 11) +
      theme(panel.grid = element_blank(), axis.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(size = 7),
            strip.background = element_rect(color = "grey88", fill = "grey88"),
            strip.text = element_text(size = 10, colour = "black"),
            legend.title = element_text(size = 9), legend.text = element_text(size = 8))
    if (length(panels)) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
    panels[[length(panels) + 1L]] <- p
    widths <- c(widths, length(unique(cur$feature)))
  }
  Reduce(`+`, panels) + plot_layout(ncol = length(panels), widths = widths, guides = "collect")
}

cluster_counts <- as.data.frame(table(obj@meta.data[[cluster_col]]), stringsAsFactors = FALSE)
colnames(cluster_counts) <- c("audit_cluster", "n_cells")
cluster_counts$study <- "Chia"
cluster_counts$cluster_source <- "Original project final RNA-graph clustering (seurat_clusters; 47 clusters; active identities)"
cluster_counts$reviewer_final_CellType <- ""
cluster_counts$reviewer_notes <- ""
cluster_counts$review_status <- "pending"
cluster_counts <- cluster_counts[, c("study", "audit_cluster", "n_cells", "cluster_source", "reviewer_final_CellType", "reviewer_notes", "review_status")]

umap_plot <- DimPlot(
  obj, reduction = "umap", group.by = cluster_col, label = TRUE, repel = TRUE,
  pt.size = 0.32, raster = FALSE, label.size = 4
) +
  coord_fixed() +
  ggtitle("Chia: final RNA-graph clusters (47; labels pending reviewer assignment)") +
  theme_classic(base_size = 15) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 18))

note_plot <- ggplot() +
  annotate("text", x = 0, y = 0.95, hjust = 0, vjust = 1, size = 5.2,
           label = paste(
             "Review scope", "",
             "• Source: historical Chia.qs, read-only", "• UMAP: original stored coordinates", "• Audit units: final RNA-graph seurat_clusters (47)",
             "• Marker evidence: exact 49-gene Yao panel", "• Existing cell.type is shown on page 3 only as context", "• No CellType assignment is made by this script", "",
             "Use the accompanying CSV to record cluster → CellType decisions.", sep = "\n")) +
  xlim(0, 1) + ylim(0, 1) + theme_void()

dot_plot <- make_dotplot(obj, marker_tbl, cluster_col) +
  plot_annotation(title = "Chia: Yao marker evidence by final RNA-graph cluster", subtitle = "Color = scaled average RNA expression; dot size = percentage of cells expressing the gene")

author_plot <- DimPlot(
  obj, reduction = "umap", group.by = label_col, label = TRUE, repel = TRUE,
  pt.size = 0.32, raster = FALSE, label.size = 4
) +
  coord_fixed() +
  ggtitle("Chia: existing cell.type metadata (context only; not adopted as reviewer decision)") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 18))

pdf_file <- file.path(out_dir, "Chia_final_RNAcluster_Yao_CellType_reviewer_evidence.pdf")
grDevices::pdf(pdf_file, width = 32, height = 14, onefile = TRUE, useDingbats = FALSE)
print(umap_plot + note_plot + plot_layout(widths = c(3, 1)))
print(dot_plot)
print(author_plot)
grDevices::dev.off()

ggsave(file.path(out_dir, "Chia_final_RNAcluster_UMAP.png"), umap_plot, width = 15, height = 13, dpi = 240, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "Chia_final_RNAcluster_Yao_DotPlot.png"), dot_plot, width = 32, height = 14, dpi = 200, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "Chia_existing_celltype_context.png"), author_plot, width = 15, height = 13, dpi = 240, bg = "white", limitsize = FALSE)
write.csv(cluster_counts, file.path(out_dir, "tables", "Chia_reviewer_cluster_assignment_template.csv"), row.names = FALSE)
write.csv(marker_tbl[, c("Type", "Gname", "feature", "present")], file.path(out_dir, "tables", "Chia_Yao_marker_panel.csv"), row.names = FALSE)
write.csv(data.frame(
  source_object = source_path,
  audit_cluster_column = cluster_col,
  active_identity_check = TRUE,
  reduction = "umap",
  marker_genes_present = paste0(sum(marker_tbl$present), "/", nrow(marker_tbl)),
  source_object_modified = FALSE,
  stringsAsFactors = FALSE
), file.path(out_dir, "tables", "Chia_corrected_audit_provenance.csv"), row.names = FALSE)
writeLines(c(
  "# Corrected Chia CellType reviewer evidence", "",
  "This supersedes the earlier Chia SCT-25 audit only.",
  "Reason: SCT_snn_res.0.8 contains 25 broad mixed clusters. The historical final active partition is seurat_clusters / RNA_snn_res.0.8 (47 clusters).",
  "No original object, metadata, UMAP, or cell label was saved or altered.",
  "The reviewer assignment template is intentionally blank. Existing cell.type labels are displayed in the PDF only as metadata context and are not propagated into reviewer assignments."
), file.path(out_dir, "README.md"))
