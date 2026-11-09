# Verification copy of the original cross-species plotting script.
# Only server paths, output names, and the missing object assignment below are
# repaired; plotting logic and labels are retained verbatim.
setwd(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies"))


library(qs)
library(Seurat)

manual_colors <- c(
  "MG2" = "#377EB8", "MG5" = "#8B0000", "MG0" = "#4DAF4A", 
  "MG3" = "#984EA3", "MG4" = "#E6AB02", "MG1" = "#A65628"
)


combine <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/02.cross_species_integrated_noref_Human_Mouse_Monkey.qs"))
# The original file reads this object as `combine` but uses `integrated_mg`
# below without creating it.  The figure code requires this alias.
integrated_mg <- combine
DimPlot(combine,group.by = 'CellType1', cols = manual_colors, split.by = 'Species')
DimPlot(combine,group.by = 'integrated_snn_res.0.1')


integrated_mg$new_celltype <- as.character(integrated_mg$CellType1)
target_cells <- integrated_mg$integrated_snn_res.0.1 == '1' & integrated_mg$Species == 'Mouse'
integrated_mg$new_celltype[target_cells] <- "Mouse_MG0_like"
Idents(integrated_mg) <- "new_celltype"
table(integrated_mg$new_celltype)

target_cells <- integrated_mg$integrated_snn_res.0.1 == '1' & integrated_mg$Species == 'Monkey'
integrated_mg$new_celltype[target_cells] <- "Monkey_MG0_like"
Idents(integrated_mg) <- "new_celltype"
table(integrated_mg$new_celltype)


#-------------
keep_types <- c("MG0", "Mouse_MG0_like", "Monkey_MG0_like")
integrated_mg$new_celltype <- as.character(integrated_mg$new_celltype)
integrated_mg$new_celltype <- ifelse(
  integrated_mg$new_celltype %in% keep_types,
  integrated_mg$new_celltype, 
  "Other_MG"
)
Idents(integrated_mg) <- "new_celltype"
table(integrated_mg$new_celltype, useNA = "ifany")

DimPlot(integrated_mg,group.by = 'new_celltype',split.by = 'Species')


#-----------------
# 1. 创建命名颜色向量 (Named Vector)
custom_colors <- c(
  "MG0"             = "#4DAF4A",
  "Mouse_MG0_like"  = "#33A02C",  
  "Monkey_MG0_like" = "#B2DF8A",  
  "Other_MG"        = "#D9D9D9"   
)

library(ggplot2)

(p1 <- DimPlot(integrated_mg, 
        reduction = "umap", 
        group.by = "new_celltype", 
        cols = custom_colors,
        split.by = 'Species',
        pt.size = 0.3,       # 由于 Other_MG 细胞很多，适当调整点的大小
        raster = FALSE) + 
  ggtitle("MG0 Cross-Species Subsets") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoAxes()
)
ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/figure/species_dimplot_figureblock_reproduced_24server.pdf"),p1,width = 9,height = 3)
ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/figure/species_dimplot_figureblock_reproduced_24server.png"),p1,width = 9,height = 3)

# 1. 为三个物种分配 3 种区分度高的经典颜色
species_colors <- c(
  "Human"  = "#E69F00",  
  "Mouse"  = "#56B4E9",  
  "Monkey" = "#009E73"   
)

# 2. 画在同一个 UMAP 当中
(p2 <- DimPlot(integrated_mg, 
              reduction = "umap", 
              group.by = "Species",     # 确保这里的拼写是 Species
              cols = species_colors,    # 应用上面定义的 3 个颜色
              pt.size = 0.2,            # 调整点的大小
              raster = FALSE) + 
  ggtitle("Species") +
  theme(plot.title = element_text(hjust = 0.5)) +
  NoAxes()                              # 取消横纵坐标及刻度
)

ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/figure/species_celltype_figureblock_reproduced_24server.pdf"),p2,width = 4,height = 3)
ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/results/figure/species_celltype_figureblock_reproduced_24server.png"),p2,width = 4,height = 3)

message('SUCCESS: original cross-species figure block reproduced.')
