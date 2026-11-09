library(qs)
library(Seurat)

# Mouse findallmarkers----------------
mouse_wm_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/mouse_wm_mg_v2.qs"))
DefaultAssay(mouse_wm_mg) <- "RNA"
Idents(mouse_wm_mg) <- "new_celltype"
mouse_mg_markers <- FindAllMarkers(mouse_wm_mg, 
                                   only.pos = TRUE, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.25)
qsave(mouse_mg_markers,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/mouse_mg_markers.qs"))


# Monkey findallmarkers----------------
monkey_wm_mg <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/monkey_hip_mg_v2.qs"))
DefaultAssay(monkey_wm_mg) <- "RNA"
Idents(monkey_wm_mg) <- "new_celltype"
monkey_mg_markers <- FindAllMarkers(monkey_wm_mg, 
                                   only.pos = TRUE, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.25)
qsave(monkey_mg_markers,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/MultiSpecies/data/processed_data/monkey_mg_markers.qs"))
