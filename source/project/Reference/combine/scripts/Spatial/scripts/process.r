library(qs)

merfish <- readRDS(Sys.getenv("MG_EXTERNAL_026"))

merfish$plot_celltype <- as.character(merfish$cell.type_manual)
merfish$plot_celltype[merfish$plot_celltype %in% c('Myelin_MG')] <- 'ST_MG0'  
merfish$plot_celltype[merfish$plot_celltype %in% c('MG_0')] <- 'ST_MG1'  
merfish$plot_celltype[merfish$plot_celltype %in% c('MG_1')] <- 'ST_MG2'  
merfish$plot_celltype[merfish$plot_celltype %in% c('MG_2')] <- 'ST_MG3'  
merfish$plot_celltype[merfish$plot_celltype %in% c('MG_3')] <- 'ST_MG4'  
merfish$plot_celltype[merfish$plot_celltype %in% c('MG_5')] <- 'ST_MG5'  

saveRDS(merfish,file = Sys.getenv("MG_EXTERNAL_027"))
merfish <- readRDS(Sys.getenv("MG_EXTERNAL_027"))
merfish_combined <- readRDS(Sys.getenv("MG_EXTERNAL_028"))

merfish_combined$cell.type_manual2 = merfish_combined$cell.type_manual
merfish_combined$cell.type_manual2[which(merfish_combined$cell.type_manual2 == 'MG')] = merfish$plot_celltype

unique(merfish_combined$cell.type_manual2)

saveRDS(merfish_combined,Sys.getenv("MG_EXTERNAL_022"))

#----------------------------
merfish <- readRDS(Sys.getenv("MG_EXTERNAL_022"))
merfish_s2r1 <- readRDS(Sys.getenv("MG_EXTERNAL_023"))

# 把细胞类型从整合的数据映射回s2r1单独的数据中---------------
all_celltype <- as.data.frame(merfish$cell.type_manual2)
colnames(all_celltype) <- 'CellType'
s2r1_celltype <- all_celltype[rownames(all_celltype) %in% colnames(merfish_s2r1), , drop = FALSE]
merfish_s2r1$CellType <- s2r1_celltype
saveRDS(merfish_s2r1,'data/save_data/seperate/merfish_s2r1_v2.rds')
