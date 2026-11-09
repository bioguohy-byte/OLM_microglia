#!/usr/bin/env Rscript

# GSE148676 microglia-candidate object:
# normalized Trim2/Rnf220 expression plus oligodendroglial markers.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(Seurat)
})

project <- Sys.getenv("MG_LPC_ROOT")
object_path <- file.path(
  project,
  "objects",
  "simple_microglia",
  "mouse_GSE148676_simple_microglia_seurat_v5.rds"
)
figure_dir <- file.path(project, "figures", "simple_microglia")
individual_figure_dir <- file.path(
  figure_dir,
  "mouse_GSE148676_individual_featureplots"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_figure_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading ", object_path)
microglia <- readRDS(object_path)
DefaultAssay(microglia) <- "RNA"

requested_genes <- c("Trim2", "Rnf220", "St18", "Olig2", "Mbp", "Mog")
available_genes <- rownames(microglia)
matched_genes <- vapply(
  requested_genes,
  function(gene) {
    hits <- available_genes[toupper(available_genes) == toupper(gene)]
    if (length(hits) == 0) {
      stop("Gene is absent from the RNA assay: ", gene)
    }
    hits[[1]]
  },
  character(1)
)
names(matched_genes) <- requested_genes
message(
  "Matched features: ",
  paste(names(matched_genes), matched_genes, sep = "=", collapse = ", ")
)

umap_frame <- as.data.frame(Embeddings(microglia, "umap"))
colnames(umap_frame)[1:2] <- c("UMAP_1", "UMAP_2")
cell_order <- rownames(umap_frame)

normalized_matrix <- LayerData(
  microglia[["RNA"]],
  layer = "data"
)[matched_genes, cell_order, drop = FALSE]

expression_values <- list()
gene_plots <- list()
point_size <- 0.20

for (gene_label in requested_genes) {
  feature_name <- matched_genes[[gene_label]]
  expression <- as.numeric(normalized_matrix[feature_name, ])
  expression_values[[gene_label]] <- expression

  plot_frame <- umap_frame
  plot_frame$expression <- expression
  positive_frame <- plot_frame[plot_frame$expression > 0, , drop = FALSE]
  positive_frame <- positive_frame[
    order(positive_frame$expression),
    ,
    drop = FALSE
  ]

  if (nrow(positive_frame) == 0) {
    upper_limit <- 1
  } else {
    upper_limit <- as.numeric(
      quantile(
        positive_frame$expression,
        0.99,
        names = FALSE
      )
    )
  }

  gene_plots[[gene_label]] <- ggplot(
    plot_frame,
    aes(x = UMAP_1, y = UMAP_2)
  ) +
    geom_point(
      color = "#B8B8B8",
      size = point_size,
      alpha = 0.60,
      shape = 16,
      stroke = 0
    ) +
    geom_point(
      data = positive_frame,
      aes(color = expression),
      size = point_size * 1.55,
      alpha = 1,
      shape = 16,
      stroke = 0
    ) +
    scale_color_gradient(
      low = "#FFD166",
      high = "#B2182B",
      limits = c(0, upper_limit),
      oob = scales::squish
    ) +
    labs(
      title = gene_label,
      color = "Log-normalized\nexpression",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    coord_equal() +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
}

joint_expression <- sqrt(
  expression_values[["Trim2"]] * expression_values[["Rnf220"]]
)
joint_frame <- umap_frame
joint_frame$expression <- joint_expression
joint_positive <- joint_frame[
  joint_frame$expression > 0,
  ,
  drop = FALSE
]
joint_positive <- joint_positive[
  order(joint_positive$expression),
  ,
  drop = FALSE
]
joint_upper <- if (nrow(joint_positive) == 0) {
  1
} else {
  as.numeric(
    quantile(joint_positive$expression, 0.99, names = FALSE)
  )
}

joint_plot <- ggplot(
  joint_frame,
  aes(x = UMAP_1, y = UMAP_2)
) +
  geom_point(
    color = "#B8B8B8",
    size = point_size,
    alpha = 0.60,
    shape = 16,
    stroke = 0
  ) +
  geom_point(
    data = joint_positive,
    aes(color = expression),
    size = point_size * 1.75,
    alpha = 1,
    shape = 16,
    stroke = 0
  ) +
  scale_color_gradient(
    low = "#FFD166",
    high = "#B2182B",
    limits = c(0, joint_upper),
    oob = scales::squish
  ) +
  labs(
    title = "Trim2 + Rnf220",
    color = "Joint log-normalized\nexpression",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  coord_equal() +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

# Export every marker panel as a standalone figure. These files retain the
# same UMAP coordinates and normalized RNA data used in the combined panel.
for (gene_label in requested_genes) {
  individual_path <- file.path(
    individual_figure_dir,
    paste0(
      "mouse_GSE148676_",
      gene_label,
      "_log_normalized_FeaturePlot.png"
    )
  )
  ggsave(
    individual_path,
    gene_plots[[gene_label]],
    width = 7.2,
    height = 6.4,
    dpi = 300,
    bg = "white"
  )
  message("Wrote ", individual_path)

  individual_pdf_path <- sub("\\.png$", ".pdf", individual_path)
  ggsave(
    individual_pdf_path,
    gene_plots[[gene_label]],
    width = 7.2,
    height = 6.4,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  message("Wrote ", individual_pdf_path)
}

joint_output_path <- file.path(
  individual_figure_dir,
  paste0(
    "mouse_GSE148676_Trim2_Rnf220_joint_",
    "log_normalized_FeaturePlot.png"
  )
)
ggsave(
  joint_output_path,
  joint_plot,
  width = 7.2,
  height = 6.4,
  dpi = 300,
  bg = "white"
)
message("Wrote ", joint_output_path)

joint_pdf_path <- sub("\\.png$", ".pdf", joint_output_path)
ggsave(
  joint_pdf_path,
  joint_plot,
  width = 7.2,
  height = 6.4,
  device = grDevices::cairo_pdf,
  bg = "white"
)
message("Wrote ", joint_pdf_path)

target_row <- (
  gene_plots[["Trim2"]] |
    gene_plots[["Rnf220"]] |
    joint_plot
)
lineage_row <- (
  gene_plots[["St18"]] |
    gene_plots[["Olig2"]] |
    gene_plots[["Mbp"]] |
    gene_plots[["Mog"]]
)
expanded_plot <- (
  target_row / lineage_row
) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = paste(
      "mouse_GSE148676 log-normalized expression:",
      "Trim2/Rnf220 and oligodendroglial markers"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 18)
    )
  )

output_path <- file.path(
  figure_dir,
  paste0(
    "mouse_GSE148676_02B_normalized_",
    "Trim2_Rnf220_St18_Olig2_Mbp_Mog_FeaturePlot.png"
  )
)
ggsave(
  output_path,
  expanded_plot,
  width = 18,
  height = 10.125,
  dpi = 300,
  bg = "white"
)

message("Wrote ", output_path)
