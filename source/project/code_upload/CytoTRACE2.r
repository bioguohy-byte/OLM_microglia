# CytoTRACE2 ---------------------
library(CytoTRACE2) 
library(qs)
library(Seurat)
library(ggplot2)

data <- qread('data/processed_data/mg_cca_meta_v2.qs')

cytotrace2_res <- cytotrace2(data, 
                             is_seurat = TRUE,
                             slot_type = "counts", 
                             species = 'human')

qsave(cytotrace2_res, 'scripts/CytoTRACE2/data/processed_data/cytotrace2_result.qs')

# Fig.01 ---------------
p1 <- FeaturePlot(cytotrace2_res, features = 'CytoTRACE2_Score') +
  scale_color_distiller(palette = "Blues", direction = 1)

ggsave('scripts/CytoTRACE2/results/fig/score.pdf', plot = p1, width = 4, height = 3.5)

# Fig.02 ---------------
p2 <- VlnPlot(cytotrace2_res, features = "CytoTRACE2_Score", group.by = "CellType1", pt.size = 0)

ggsave('scripts/CytoTRACE2/results/fig/score_vln.pdf', plot = p2, width = 6, height = 4)