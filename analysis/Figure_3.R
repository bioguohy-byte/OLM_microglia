#!/usr/bin/env Rscript
# Figure-ordered code from the reviewed manuscript analysis package.
# Each panel contains its original quoted R code.
# Source this chapter, then use mg_list_panels() or mg_show_panel("PANEL_ID").
# Source paths in panel records identify the original analysis scripts.
.mg_panels <- new.env(parent=emptyenv())
.mg_order <- character()
mg_register_panel <- function(id, code) {
  if (!id %in% .mg_order) .mg_order <<- c(.mg_order,id)
  assign(id, code, envir=.mg_panels)
}
mg_list_panels <- function() .mg_order
mg_panel_code <- function(id) {
  if (!exists(id,envir=.mg_panels,inherits=FALSE)) stop("Unknown panel: ",id)
  get(id,envir=.mg_panels,inherits=FALSE)
}
mg_show_panel <- function(id) {
  blocks <- mg_panel_code(id)
  if (is.null(blocks)) {cat(id, ": No R block is included for this panel.\n");return(invisible(NULL))}
  for (block in blocks) {
    cat("# Source: ",block$source," [actual ",block$ranges,"; requested ",block$requested_ranges,"]\n# Role: ",block$role,"\n",sep="")
    cat(paste(deparse(block$code,width.cutoff=120L),collapse="\n"),"\n\n")
  }
  invisible(blocks)
}


# Figure 3

# F3A  Human/monkey/mouse integrated microglial UMAP, species-colored UMAP and OLM-like highlights; decorative species icons and manually drawn red circles
# Kind: computational
mg_register_panel("F3A", list(
  list(source="source/project/code_upload/MultiSpecies.r", ranges="7-60", requested_ranges="7-60", role="", code=quote({
manual_colors <- c("MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A", 
                   "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628")

custom_colors <- c("MG0" = "#4DAF4A", "Mouse_MG0_like" = "#33A02C",  
                   "Monkey_MG0_like" = "#B2DF8A", "Other_MG" = "#D9D9D9")

species_colors <- c("Human" = "#E69F00", "Mouse" = "#56B4E9", "Monkey" = "#009E73")

integrated_mg <- qread('data/processed_data/02.cross_species_integrated_noref_Human_Mouse_Monkey.qs')

integrated_mg$new_celltype <- as.character(integrated_mg$CellType1)

target_cells_mouse <- integrated_mg$integrated_snn_res.0.1 == '1' & integrated_mg$Species == 'Mouse'
integrated_mg$new_celltype[target_cells_mouse] <- "Mouse_MG0_like"

target_cells_monkey <- integrated_mg$integrated_snn_res.0.1 == '1' & integrated_mg$Species == 'Monkey'
integrated_mg$new_celltype[target_cells_monkey] <- "Monkey_MG0_like"

keep_types <- c("MG0", "Mouse_MG0_like", "Monkey_MG0_like")
integrated_mg$new_celltype <- ifelse(
  integrated_mg$new_celltype %in% keep_types,
  integrated_mg$new_celltype, 
  "Other_MG"
)
Idents(integrated_mg) <- "new_celltype"

# Fig.01 ---------------
p1 <- DimPlot(integrated_mg, 
              reduction = "umap", 
              group.by = "new_celltype", 
              cols = custom_colors,
              split.by = 'Species',
              pt.size = 0.3,       
              raster = FALSE) + 
  ggtitle("MG0 Cross-Species Subsets") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoAxes()

ggsave('results/figure/species_dimplot.pdf', p1, width = 9, height = 3)
ggsave('results/figure/species_dimplot.png', p1, width = 9, height = 3)

# Fig.02 ---------------
p2 <- DimPlot(integrated_mg, 
              reduction = "umap", 
              group.by = "Species",     
              cols = species_colors,    
              pt.size = 0.2,            
              raster = FALSE) + 
  ggtitle("Species") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoAxes()                             

ggsave('results/figure/species_celltype.pdf', p2, width = 4, height = 3)
ggsave('results/figure/species_celltype.png', p2, width = 4, height = 3)
  }))
))

# F3B  Monkey OLM-like versus other MG expression of MOG, MBP, RNF220, TRIM2
# Kind: computational
mg_register_panel("F3B", list(
  list(source="source/project/code_upload/MultiSpecies.r", ranges="63-70;103-122", requested_ranges="63-70;103-122", role="", code=quote({
monkey_hip_mg <- qread('data/processed_data/monkey_hip_mg.qs')
monkey_integrated_barcodes <- colnames(integrated_mg)[integrated_mg$Species == "Monkey"]
monkey_celltypes <- integrated_mg$new_celltype[integrated_mg$Species == "Monkey"]
orig_monkey_barcodes <- gsub("^Monkey_", "", monkey_integrated_barcodes)

monkey_anno_df <- data.frame(new_celltype = monkey_celltypes, row.names = orig_monkey_barcodes)
monkey_hip_mg <- AddMetaData(monkey_hip_mg, metadata = monkey_anno_df)
qsave(monkey_hip_mg, 'data/processed_data/monkey_hip_mg_v2.qs')

# Fig.04 ---------------
DefaultAssay(monkey_hip_mg) <- 'RNA'
monkey_hip_mg$new_celltype <- factor(monkey_hip_mg$new_celltype[, 1], levels = c("Monkey_MG0_like", "Other_MG"))

compare_groups_monkey <- setdiff(levels(monkey_hip_mg$new_celltype), "Monkey_MG0_like")
my_comparisons_monkey <- lapply(compare_groups_monkey, function(x) c(x, "Monkey_MG0_like"))

p_violin_monkey <- FeatureStatPlot(
  srt         = monkey_hip_mg, 
  group.by    = "new_celltype",
  stat.by     = c("MOG", "MBP", "RNF220", "TRIM2"), 
  add_box     = FALSE, 
  comparisons = my_comparisons_monkey,
  bg.by       = "new_celltype",
  palcolor    = custom_colors,
  bg_palcolor = custom_colors,
  ncol        = 4
)

ggsave(filename = 'results/figure/Monkey_Violin_Stat_Markers.pdf', plot = p_violin_monkey, width = 16, height = 4)
  }))
))

# F3C  Mouse OLM-like versus other MG expression of Mog, Mbp, Rnf220, Trim2
# Kind: computational
mg_register_panel("F3C", list(
  list(source="source/project/code_upload/MultiSpecies.r", ranges="73-101", requested_ranges="73-101", role="", code=quote({
mouse_wm_mg <- qread('data/processed_data/mouse_wm_mg.qs')
mouse_integrated_barcodes <- colnames(integrated_mg)[integrated_mg$Species == "Mouse"]
mouse_celltypes <- integrated_mg$new_celltype[integrated_mg$Species == "Mouse"]
orig_mouse_barcodes <- gsub("^Mouse_", "", mouse_integrated_barcodes)

mouse_anno_df <- data.frame(new_celltype = mouse_celltypes, row.names = orig_mouse_barcodes)
mouse_wm_mg <- AddMetaData(mouse_wm_mg, metadata = mouse_anno_df)
qsave(mouse_wm_mg, 'data/processed_data/mouse_wm_mg_v2.qs')

# Fig.03 ---------------
DefaultAssay(mouse_wm_mg) <- 'RNA'
mouse_wm_mg$new_celltype <- factor(mouse_wm_mg$new_celltype[, 1], levels = c("Mouse_MG0_like", "Other_MG"))

compare_groups_mouse <- setdiff(levels(mouse_wm_mg$new_celltype), "Mouse_MG0_like")
my_comparisons_mouse <- lapply(compare_groups_mouse, function(x) c(x, "Mouse_MG0_like"))

p_violin_mouse <- FeatureStatPlot(
  srt         = mouse_wm_mg, 
  group.by    = "new_celltype",
  stat.by     = c("Mog", "Mbp", "Rnf220", "Trim2"), 
  add_box     = FALSE, 
  comparisons = my_comparisons_mouse,
  bg.by       = "new_celltype",
  palcolor    = custom_colors,
  bg_palcolor = custom_colors,
  ncol        = 4
)

ggsave(filename = 'results/figure/Mouse_Violin_Stat_Markers.pdf', plot = p_violin_mouse, width = 16, height = 4)
  }))
))

# F3D  Monkey-3M Iba1/Trim2/DAPI immunofluorescence, 10um
# Kind: experimental_or_manual
mg_register_panel("F3D", NULL)

# F3E  C57-CC-LPS mouse corpus callosum Iba1/Trim2/DAPI immunofluorescence, 10um
# Kind: experimental_or_manual
mg_register_panel("F3E", NULL)

# F3F  C57-CC-LPS mouse MBP/Iba1/DAPI immunofluorescence, 2um
# Kind: experimental_or_manual
mg_register_panel("F3F", NULL)

# F3G  C57-CC-LPS mouse RNF220/Iba1/DAPI immunofluorescence, 10um
# Kind: experimental_or_manual
mg_register_panel("F3G", NULL)

# F3H  C57-LPS dual immunogold TEM, Trim2 10nm and Iba1 4nm, overview and zooms
# Kind: experimental_or_manual
mg_register_panel("F3H", NULL)
