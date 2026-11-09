cell_colors <- c(
  "Neuron"      = "#E41A1C",  # 纯正红
  "OLIG"        = "#377EB8",  # 强力蓝
  "OPC"         = "#4DAF4A",  # 鲜绿色
  "Astrocyte"   = "#00CED1",  # 深青色
  "Endothelial" = "#984EA3",  # 亮紫色
  "Immune"      = "#FF7F00",  # 亮橙色
  "VLMN"        = "#A65628",  # 深棕色
  "Ependymal"   = "#E6AB02"   # 暗金/土黄色
)


library(qs)
library(dplyr)
library(magrittr)


# cpz mouse 改一下细胞大群的名字-------------------
mouse <- qread(Sys.getenv("MG_EXTERNAL_017"))
mouse$combine_celltype<-as.character(mouse$CellType1)
mouse@meta.data %<>%
  mutate(combine_celltype=case_when(CellType1%in%c("Endothelial","VLMN")~"Vascular",
                            CellType1%in%c("Immune")~"T cell",
                            TRUE ~ combine_celltype
  ))
my_order <- c("Microglia","Neuron", "Astrocyte", "OLIG", "OPC",  "Ependymal", "Vascular", "T cell")
mouse$combine_celltype <- factor(mouse$combine_celltype, levels = my_order)
qsave(mouse,Sys.getenv("MG_EXTERNAL_017"))


# yao 改细胞大群名字--------------------
yao <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Yao_et_al./data/processed_data/Yao_v2.qs"))
yao$combine_celltype<-as.character(yao$CellType)
yao@meta.data %<>%
  mutate(combine_celltype=case_when(combine_celltype%in%c("Vasular")~"Vascular",
                                    combine_celltype%in%c("Micro")~"Microglia",
                                    combine_celltype%in%c("Neurons")~"Neuron",
                                    combine_celltype%in%c("Astro","Unknown1")~"Astrocyte",
                                    combine_celltype%in%c("Oli")~"OLIG",
                                    combine_celltype%in%c("Unknown2")~"Ependymal",
                                    combine_celltype%in%c("T cells")~"T cell",
                                    TRUE ~ combine_celltype
  ))
my_order <- c("Microglia","Neuron", "Astrocyte", "OLIG", "OPC",  "Ependymal", "Vascular", "T cell")
yao$combine_celltype <- factor(yao$combine_celltype, levels = my_order)
qsave(yao,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Yao_et_al./data/processed_data/Yao_v2.qs"))


yao_stroke <- subset(yao,subset = Sample %in% c('l_48y'))
qsave(yao_stroke,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Yao_et_al./data/processed_data/Yao_stroke.qs"))
