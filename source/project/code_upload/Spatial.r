# Spatial ----------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(spatstat)
library(qs)
library(SCP)

file_master_rds <- 'data/raw_data/merfish_integrated_v2.rds'
file_in_merfish_immune <- 'data/raw_data/merfish_immune.combined_with3.rds'

col_master_celltype <- 'cell.type_manual2'
target_colors <- c("ST_MG0" = "red", "OLG" = "blue", "OPC" = "yellow")
color_other <- "gray80"

sample_configs <- list(
  "s2r1" = list(
    in_rds  = "data/raw_data/merfish_s2r1.rds",
    out_rds = "data/processed_data/merfish_s2r1_v2.rds",
    coord   = c(x_min = -2000, x_max = -1500, y_min = -8000, y_max = -7500)
  ),
  "s2r2" = list(
    in_rds  = "data/raw_data/merfish_s2r2.rds",
    out_rds = "data/processed_data/merfish_s2r2_v2.rds",
    coord   = c(x_min = 1000, x_max = 1500, y_min = -9000, y_max = -8500)
  ),
  "s2r3" = list(
    in_rds  = "data/raw_data/merfish_s2r3.rds",
    out_rds = "data/processed_data/merfish_s2r3_v2.rds",
    coord   = c(x_min = 1000, x_max = 1500, y_min = -5500, y_max = -4500)
  )
)

master_merfish <- readRDS(file_master_rds)
master_celltype_df <- as.data.frame(master_merfish[[col_master_celltype]])
colnames(master_celltype_df) <- 'CellType'
rm(master_merfish)

if(!dir.exists('results/figure/')) dir.create('results/figure/', recursive = TRUE)
if(!dir.exists('data/processed_data/')) dir.create('data/processed_data/', recursive = TRUE)

# Fig.01 ---------------
for (smpl in names(sample_configs)) {
  cfg <- sample_configs[[smpl]]
  
  obj <- readRDS(cfg$in_rds)
  
  smpl_celltype <- master_celltype_df[rownames(master_celltype_df) %in% colnames(obj), , drop = FALSE]
  obj$CellType <- smpl_celltype
  
  df_meta <- obj@meta.data %>% select(rotated_x, rotated_y, CellType)
  cells_keep <- df_meta %>% 
    filter(rotated_x >= cfg$coord["x_min"] & rotated_x <= cfg$coord["x_max"] & 
             rotated_y >= cfg$coord["y_min"] & rotated_y <= cfg$coord["y_max"]) %>% 
    rownames()
  
  obj_sub <- subset(obj, cells = cells_keep)
  
  current_celltypes <- unique(obj_sub$CellType)
  colors_mapped <- sapply(current_celltypes, function(x) {
    if (x %in% names(target_colors)) return(target_colors[[x]])
    else return(color_other)
  })
  
  p <- DimPlot(obj_sub, reduction = 'spatial', group.by = 'CellType', pt.size = 0.5) + 
    scale_color_manual(values = colors_mapped) +
    ggtitle(paste("Spatial Plot -", smpl))
  
  saveRDS(obj, file = cfg$out_rds)
  
  file_base <- paste0('results/figure/', smpl, "_plot_zoom")
  ggsave(paste0(file_base, ".pdf"), plot = p, width = 4, height = 4)
  ggsave(paste0(file_base, ".png"), plot = p, width = 4, height = 4)
}

merfish <- readRDS(file_in_merfish_immune)
mouse_mg_markers <- qread('../MultiSpecies/data/saved_data/mouse_mg_markers.qs')

celltype_dict <- c(
  'Myelin_MG' = 'ST_MG0',
  'MG_0'      = 'ST_MG1',
  'MG_1'      = 'ST_MG2',
  'MG_2'      = 'ST_MG3',
  'MG_3'      = 'ST_MG4',
  'MG_5'      = 'ST_MG5'
)

merfish@meta.data[['plot_celltype']] <- as.character(merfish@meta.data[['cell.type_manual']])
for (old_name in names(celltype_dict)) {
  new_name <- celltype_dict[[old_name]]
  merfish@meta.data[['plot_celltype']][merfish@meta.data[['plot_celltype']] %in% old_name] <- new_name
}

filtered_data <- mouse_mg_markers %>% filter(.data[['cluster']] == 'cluster2')
int_gene <- intersect(filtered_data[['gene']], rownames(merfish))

merfish <- AddModuleScore(object = merfish,
                          features = list(int_gene),
                          ctrl = 10,
                          name = "myelin")

saveRDS(merfish, file = 'data/processed_data/merfish_immune.combined_with4.rds')

st_mg_colors_npg <- c(
  "ST_MG0" = "#E64B35", 
  "ST_MG1" = "#4DBBD5", 
  "ST_MG2" = "#00A087", 
  "ST_MG3" = "#3C5488", 
  "ST_MG4" = "#F39B7F", 
  "ST_MG5" = "#7E6148"
)

# Fig.02 ---------------
p1 <- DimPlot(merfish, group.by = 'plot_celltype', cols = st_mg_colors_npg)

target_order <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")
merfish$plot_celltype <- factor(merfish$plot_celltype, levels = target_order)

celltype <- levels(merfish$plot_celltype)
my_comparisons <- lapply(celltype[-1], function(x) c(x, celltype[1]))

p2 <- FeatureStatPlot(
  srt = merfish, group.by = "plot_celltype",
  stat.by = c("myelin1"), add_box = FALSE,
  comparisons = my_comparisons,
  bg.by = 'plot_celltype',
  palcolor = st_mg_colors_npg,
  bg_palcolor = st_mg_colors_npg
)

combined_plot <- p1 | p2
ggsave(filename = 'results/figure/module_score_myelin.pdf', plot = combined_plot, width = 10, height = 5)

obj_l <- readRDS("data/processed_data/merfish_s2r1_v2.rds")
coord_filter <- c(x_min = -3000, x_max = -1500, y_min = -8000, y_max = -5300)

df_l <- obj_l@meta.data %>% 
  select(rotated_x, rotated_y, CellType) %>%
  filter(rotated_x >= coord_filter["x_min"] & rotated_x <= coord_filter["x_max"] & 
           rotated_y >= coord_filter["y_min"] & rotated_y <= coord_filter["y_max"])

df_l$rotated_x <- as.numeric(df_l$rotated_x)
df_l$rotated_y <- as.numeric(df_l$rotated_y)

pp <- ppp(df_l$rotated_x, df_l$rotated_y, 
          window = owin(range(df_l$rotated_x), range(df_l$rotated_y)), 
          marks = factor(df_l$CellType))

centered_L <- function(pp, i, j, r_max, step, nsim) {
  r_seq <- seq(0, r_max, by = step)
  L_sim <- envelope(pp, Lcross, i = i, j = j, nsim = nsim, r = r_seq, 
                    correction = "border", simulate = expression(rlabel(pp)),
                    verbose = FALSE) 
  
  if (is.null(L_sim$lo) || is.null(L_sim$hi) || is.null(L_sim$mmean)) return(NULL)
  
  data.frame(
    r  = r_seq, 
    L  = L_sim$obs - L_sim$mmean, 
    lo = L_sim$lo - L_sim$mmean, 
    hi = L_sim$hi - L_sim$mmean
  )
}

center_cell <- "OLG"
target_cells <- c("ST_MG0", "ST_MG1", "ST_MG2", "ST_MG3", "ST_MG4", "ST_MG5")
results <- list()

for (target in target_cells) {
  if (sum(df_l$CellType == center_cell) >= 10 && sum(df_l$CellType == target) >= 10) {
    res <- centered_L(pp, center_cell, target, r_max = 250, step = 1, nsim = 200)
    if (!is.null(res)) {
      res$cell_type1 <- center_cell
      res$cell_type2 <- target
      results[[paste(center_cell, target, sep = "_vs_")]] <- res
    }
  }
}

df_combined <- bind_rows(results)
df_filtered <- df_combined %>% filter(r >= 5)

# Fig.03 ---------------
p_l_func <- ggplot(df_filtered, aes(x = r)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = cell_type2), alpha = 0.3) +
  geom_line(aes(y = L, color = cell_type2), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ cell_type2, scales = "fixed", ncol = length(target_cells)) + 
  labs(x = "Distance (μm)", y = "Centered L(r)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1) 
  ) +
  scale_color_manual(values = st_mg_colors_npg) +
  scale_fill_manual(values = st_mg_colors_npg)

ggsave(filename = "results/figure/s2r1_OLG_vs_other_MG_subclusters_L_function_plot_single_row.pdf", 
       plot = p_l_func, width = 12, height = 3.5, units = "in", device = cairo_pdf, dpi = 300)