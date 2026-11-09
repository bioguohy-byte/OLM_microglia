#!/usr/bin/env Rscript

# Provisional marker-based CellType suggestions for reviewer audit.
# Source objects are read-only. This script never writes CellType labels to a Seurat object.

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: build_provisional_Yao_CellType_suggestions.R <manifest.tsv> <out_dir>")
manifest <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "png"), recursive = TRUE, showWarnings = FALSE)
set.seed(42)

# Exact marker panel from Yao_et_al./scripts/00.yao_alldata.r, section 10.
markers <- data.frame(
  Gname = c("AQP4","ITGB4","ETNPPL","SLC1A2","RFX4","GLI3","GRAMD2B","GFAP",
            "VCAN","C1QL1","SUSD5","COL9A1",
            "MOG","MOBP","SLC5A11","PICK1","HAPLN2","CNDP1","MYRF","CPM","TF","TMC6","CDK18",
            "PTPRC","CSF1R","MRC1","IL18","C3","SLC2A5","BIN2","DOCK2","ENTPD1",
            "CFAP54","DNAH11","DCN","LAMA2","RERGL","ABCG2","EMCN","LEF1","FLT1","SLC6A12","CFH","SLCO1B3",
            "RBFOX3","SNAP25","SYT1","GRIN1","CD247"),
  CellType = c(rep("Astr", 8), rep("OPC", 4), rep("Oli", 11), rep("I", 1), rep("Micro", 8),
               rep("Ependymal", 2), rep("Vascular", 10), rep("Neurons", 4), rep("T cells", 1)),
  stringsAsFactors = FALSE
)

# Single-marker groups are reported but never proposed automatically.
multi_marker_types <- c("Astr", "OPC", "Oli", "Micro", "Ependymal", "Vascular", "Neurons")

numeric_cluster_order <- function(x) {
  x <- as.character(x)
  n <- suppressWarnings(as.numeric(x))
  if (all(!is.na(n))) as.character(sort(unique(n))) else sort(unique(x))
}

match_markers <- function(obj, marker_table) {
  features <- rownames(obj[["RNA"]])
  lookup <- setNames(features, toupper(features))
  marker_table$feature <- unname(lookup[toupper(marker_table$Gname)])
  marker_table$present <- !is.na(marker_table$feature)
  marker_table
}

score_markers <- function(dot, marker_table, study, cluster_levels) {
  joined <- merge(dot, marker_table[, c("feature", "Gname", "CellType")],
                  by.x = "features.plot", by.y = "feature", all.x = TRUE, sort = FALSE)
  joined$cluster <- factor(as.character(joined$id), levels = cluster_levels)
  joined <- joined[!is.na(joined$CellType), , drop = FALSE]
  keys <- split(seq_len(nrow(joined)), interaction(joined$cluster, joined$CellType, drop = TRUE, lex.order = TRUE))
  out <- lapply(keys, function(ii) {
    d <- joined[ii, , drop = FALSE]
    data.frame(
      study = study,
      cluster = as.character(d$cluster[[1]]),
      CellType = d$CellType[[1]],
      n_markers = nrow(d),
      mean_scaled_expression = mean(d$avg.exp.scaled, na.rm = TRUE),
      median_scaled_expression = median(d$avg.exp.scaled, na.rm = TRUE),
      mean_percent_expressed = mean(d$pct.exp, na.rm = TRUE),
      median_percent_expressed = median(d$pct.exp, na.rm = TRUE),
      supported_markers = sum(d$avg.exp.scaled >= 0.5 & d$pct.exp >= 10, na.rm = TRUE),
      support_fraction = mean(d$avg.exp.scaled >= 0.5 & d$pct.exp >= 10, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

propose_clusters <- function(stats) {
  # A group that becomes single-marker after study-specific feature matching is
  # intentionally treated as manual-review-only, just like I/PTPRC and T/CD247.
  eligible <- stats[stats$CellType %in% multi_marker_types & stats$n_markers >= 2L, , drop = FALSE]
  out <- lapply(split(eligible, eligible$cluster), function(d) {
    d <- d[order(-d$mean_scaled_expression, -d$median_scaled_expression, -d$support_fraction), , drop = FALSE]
    top <- d[1, , drop = FALSE]
    second_score <- if (nrow(d) >= 2L) d$mean_scaled_expression[[2]] else NA_real_
    min_supported <- max(2L, ceiling(top$n_markers[[1]] * 0.25))
    margin <- top$mean_scaled_expression[[1]] - second_score
    strong <- top$mean_scaled_expression[[1]] >= 0.35 &&
      top$median_scaled_expression[[1]] >= 0 &&
      top$supported_markers[[1]] >= min_supported &&
      top$mean_percent_expressed[[1]] >= 10 &&
      is.finite(margin) && margin >= 0.25
    data.frame(
      study = top$study[[1]],
      audit_cluster = top$cluster[[1]],
      proposed_CellType = if (strong) top$CellType[[1]] else "",
      top_marker_group = top$CellType[[1]],
      top_mean_scaled_expression = top$mean_scaled_expression[[1]],
      top_median_scaled_expression = top$median_scaled_expression[[1]],
      top_mean_percent_expressed = top$mean_percent_expressed[[1]],
      supported_markers = top$supported_markers[[1]],
      required_supported_markers = min_supported,
      second_best_marker_group = if (nrow(d) >= 2L) d$CellType[[2]] else "",
      second_best_mean_scaled_expression = second_score,
      score_margin = margin,
      proposal_status = if (strong) "suggested - reviewer confirmation required" else "ambiguous - manual review required",
      reviewer_final_CellType = "",
      reviewer_notes = "",
      review_status = "pending",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

make_annotated_dotplot <- function(dot, marker_table, proposal, study, cluster_levels) {
  joined <- merge(dot, marker_table[, c("feature", "Gname", "CellType")],
                  by.x = "features.plot", by.y = "feature", all.x = TRUE, sort = FALSE)
  joined$Gname <- factor(joined$Gname, levels = marker_table$Gname)
  joined$id <- factor(as.character(joined$id), levels = cluster_levels)
  limits <- range(joined$avg.exp.scaled, na.rm = TRUE)
  if (!all(is.finite(limits)) || diff(limits) == 0) limits <- c(-1, 1)
  plots <- list(); widths <- numeric()
  for (type_name in unique(marker_table$CellType)) {
    cur <- marker_table[marker_table$CellType == type_name, , drop = FALSE]
    d <- joined[joined$CellType == type_name, , drop = FALSE]
    d$Gname <- factor(d$Gname, levels = cur$Gname)
    rect_clusters <- proposal$audit_cluster[proposal$proposed_CellType == type_name & grepl("^suggested", proposal$proposal_status)]
    rect_data <- if (length(rect_clusters)) {
      data.frame(
        xmin = rep(0.5, length(rect_clusters)),
        xmax = rep(length(unique(cur$Gname)) + 0.5, length(rect_clusters)),
        ymin = match(rect_clusters, cluster_levels) - 0.46,
        ymax = match(rect_clusters, cluster_levels) + 0.46,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(xmin = numeric(), xmax = numeric(), ymin = numeric(), ymax = numeric())
    }
    rect_data <- rect_data[is.finite(rect_data$ymin), , drop = FALSE]
    p <- ggplot(d, aes(x = Gname, y = id, size = pct.exp, colour = avg.exp.scaled)) +
      geom_point() +
      scale_size_continuous(range = c(0.8, 3.8), limits = c(0, 100)) +
      scale_colour_gradientn(colors = rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", "#ffe6b7", "grey90", "white")), limits = limits) +
      facet_wrap(~CellType) +
      guides(colour = guide_colourbar(title = "Average expression"), size = guide_legend(title = "Percent expressed")) +
      theme_classic(base_size = 11) +
      theme(panel.grid = element_blank(), axis.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(size = if (length(cluster_levels) > 30) 6 else 8),
            strip.background = element_rect(color = "grey88", fill = "grey88"),
            strip.text = element_text(size = 10, colour = "black"),
            legend.title = element_text(size = 9), legend.text = element_text(size = 8))
    if (nrow(rect_data)) p <- p + geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                             inherit.aes = FALSE, fill = NA, colour = "#E31A1C", linewidth = 0.8)
    if (length(plots)) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
    plots[[length(plots) + 1L]] <- p
    widths <- c(widths, length(unique(cur$Gname)))
  }
  Reduce(`+`, plots) + plot_layout(ncol = length(plots), widths = widths, guides = "collect") +
    plot_annotation(
      title = paste0(study, ": Yao-marker DotPlot with provisional suggestions"),
      subtitle = "Red outline: a multi-marker statistical suggestion only; reviewer confirmation required. Unboxed rows are ambiguous or require manual review."
    )
}

all_stats <- list(); all_proposals <- list(); all_counts <- list(); all_marker_availability <- list(); page_plots <- list()
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, , drop = FALSE]
  study <- row$study[[1]]
  message("Processing ", study)
  obj <- qread(row$object_path[[1]])
  cluster_col <- row$cluster_col[[1]]
  if (!cluster_col %in% colnames(obj@meta.data)) stop(study, ": missing cluster column ", cluster_col)
  if (!row$reduction[[1]] %in% Reductions(obj)) stop(study, ": missing stored UMAP")
  if (!"RNA" %in% Assays(obj)) stop(study, ": RNA assay absent")
  DefaultAssay(obj) <- "RNA"
  marker_table <- match_markers(obj, markers)
  marker_table$study <- study
  all_marker_availability[[study]] <- marker_table[, c("study", "CellType", "Gname", "feature", "present")]
  marker_table <- marker_table[marker_table$present, , drop = FALSE]
  if (!nrow(marker_table)) stop(study, ": no Yao marker is present in the RNA assay.")
  levels_use <- numeric_cluster_order(obj@meta.data[[cluster_col]])
  obj@meta.data[[cluster_col]] <- factor(as.character(obj@meta.data[[cluster_col]]), levels = levels_use)
  dot <- DotPlot(obj, features = marker_table$feature, group.by = cluster_col, scale = TRUE)$data
  stats <- score_markers(dot, marker_table, study, levels_use)
  proposal <- propose_clusters(stats)
  proposal <- proposal[match(levels_use, proposal$audit_cluster), , drop = FALSE]
  proposal$n_cells <- as.integer(table(obj@meta.data[[cluster_col]])[proposal$audit_cluster])
  proposal$cluster_source <- row$cluster_source[[1]]
  proposal <- proposal[, c("study", "audit_cluster", "n_cells", "cluster_source", "proposed_CellType", "top_marker_group",
                           "top_mean_scaled_expression", "top_median_scaled_expression", "top_mean_percent_expressed",
                           "supported_markers", "required_supported_markers", "second_best_marker_group",
                           "second_best_mean_scaled_expression", "score_margin", "proposal_status",
                           "reviewer_final_CellType", "reviewer_notes", "review_status")]
  p_umap <- DimPlot(obj, reduction = row$reduction[[1]], group.by = cluster_col, label = TRUE, repel = TRUE,
                    pt.size = 0.28, raster = FALSE, label.size = if (length(levels_use) > 30) 3.5 else 4) +
    coord_fixed() + ggtitle(paste0(study, ": audit clusters")) +
    theme_classic(base_size = 14) + theme(legend.position = "none", plot.title = element_text(face = "bold", size = 16))
  p_dot <- make_annotated_dotplot(dot, marker_table, proposal, study, levels_use)
  note_text <- paste0(
    "Proposal summary\n\nSuggested: ", sum(grepl("^suggested", proposal$proposal_status)), "/", nrow(proposal), " clusters\n",
    "Ambiguous/manual: ", sum(grepl("^ambiguous", proposal$proposal_status)), " clusters\n\n",
    "Statistics:\n• mean / median scaled expression across each marker group\n• percentage of cells expressing markers\n• number of markers with scaled expression ≥0.5 and pct.exp ≥10\n• margin from the second-best multi-marker group\n\n",
    "Marker availability in this study: ", nrow(marker_table), "/", nrow(markers), ".\n",
    "Automatic suggestions exclude I/PTPRC, T-cell/CD247 and any marker group reduced to one available gene.\n",
    "Every proposal requires reviewer confirmation; no source-object label is changed."
  )
  p_note <- ggplot() + annotate("text", x = 0, y = 1, hjust = 0, vjust = 1, label = note_text, size = 4.2) +
    xlim(0, 1) + ylim(0, 1) + theme_void()
  page_plots[[study]] <- list(umap = p_umap, dot = p_dot, note = p_note)
  ggsave(file.path(out_dir, "png", paste0(study, "_provisional_suggestions.png")), p_umap + p_dot + p_note + plot_layout(widths = c(1.2, 4.2, 1.1)),
         width = 34, height = if (length(levels_use) > 30) 14 else 10, dpi = 180, bg = "white", limitsize = FALSE)
  all_stats[[study]] <- stats
  all_proposals[[study]] <- proposal
  all_counts[[study]] <- data.frame(study = study, audit_cluster = levels_use, n_cells = as.integer(table(obj@meta.data[[cluster_col]])[levels_use]), stringsAsFactors = FALSE)
  rm(obj); gc()
}

pdf_path <- file.path(out_dir, "All_studies_provisional_Yao_CellType_suggestions.pdf")
grDevices::pdf(pdf_path, width = 34, height = 14, onefile = TRUE, useDingbats = FALSE)
for (study in names(page_plots)) {
  p <- page_plots[[study]]$umap + page_plots[[study]]$dot + page_plots[[study]]$note + plot_layout(widths = c(1.2, 4.2, 1.1))
  print(p)
}
grDevices::dev.off()

stats_out <- do.call(rbind, all_stats)
proposal_out <- do.call(rbind, all_proposals)
write.csv(stats_out, file.path(out_dir, "tables", "cluster_marker_group_statistics.csv"), row.names = FALSE)
write.csv(proposal_out, file.path(out_dir, "tables", "provisional_cluster_CellType_suggestions.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_marker_availability), file.path(out_dir, "tables", "marker_availability_by_study.csv"), row.names = FALSE)
write.csv(data.frame(
  marker_group = unique(markers$CellType),
  n_markers = as.integer(table(markers$CellType)[unique(markers$CellType)]),
  automatic_proposal_eligible = unique(markers$CellType) %in% multi_marker_types,
  stringsAsFactors = FALSE
), file.path(out_dir, "tables", "proposal_algorithm_marker_groups.csv"), row.names = FALSE)
writeLines(c(
  "# Provisional Yao-marker CellType suggestions", "",
  "Source Seurat objects were opened read-only. No labels were written to source metadata or persisted objects.",
  "For each study and cluster, the statistics table reports the mean/median per-gene scaled DotPlot expression, expression prevalence, multi-marker support count and score margin.",
  "A suggested label requires all of: mean scaled expression >= 0.35; median scaled expression >= 0; mean prevalence >= 10%; at least max(2, ceil(25% of group markers)) supported markers (scaled expression >=0.5 and pct.exp >=10); and a >=0.25 mean-score margin over the second-best multi-marker group.",
  "The one-marker I/PTPRC and T-cell/CD247 groups, and any group reduced to one present feature in a study, are never auto-suggested and require manual review.",
  "Red outlines on DotPlots identify suggested clusters. Every result is provisional and requires reviewer confirmation."
), file.path(out_dir, "README.md"))
