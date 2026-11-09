library(CytoTRACE2, lib.loc=Sys.getenv("MG_EXTERNAL_018")) # loading
library(qs)

# load rds
data <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))

DefaultAssay(data)

cytotrace2_res <- cytotrace2(data, #seurat对象
                             is_seurat = TRUE,
                             slot_type = "counts", #counts和data都可以
                             species = 'human')#物种要选择，默认是小鼠
class(cytotrace2_res)


qsave(cytotrace2_res,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/CytoTRACE2/data/processed_data/cytotrace2_result.qs"))

cytotrace2_res <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/CytoTRACE2/data/processed_data/cytotrace2_result.qs"))


p1 <- FeaturePlot(cytotrace2_res, features = 'CytoTRACE2_Score') +
  scale_color_distiller(palette = "Blues", direction = 1)
ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/CytoTRACE2/results/fig/score.pdf"),p1,width=4,height=3.5)


VlnPlot(cytotrace2_res, features = "CytoTRACE2_Score", group.by = "CellType1", pt.size = 0)

# library(CytoTRACE2, lib.loc=Sys.getenv("MG_EXTERNAL_018")) # loading
# library(qs)
# 
# # load rds
# data <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/data/processed_data/mg_cca_meta_v2.qs"))
# 
# DefaultAssay(data)
# 
# cytotrace2_res <- cytotrace2(data, #seurat对象
#                              is_seurat = TRUE, 
#                              slot_type = "data", #counts和data都可以
#                              species = 'human')#物种要选择，默认是小鼠
# class(cytotrace2_res)
# 
# 
# qsave(cytotrace2_res,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/CytoTRACE2/data/processed_data/cytotrace2_result_data.qs"))
# 
# cytotrace2_res <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/CytoTRACE2/data/processed_data/cytotrace2_result_data.qs"))
# 
# 
# DimPlot(cytotrace2_res,group.by = 'CytoTRACE2_Potency')
# FeaturePlot(cytotrace2_res,features = 'CytoTRACE2_Score')
# 
# VlnPlot(cytotrace2_res, features = "CytoTRACE2_Score", group.by = "CellType1", pt.size = 0)
