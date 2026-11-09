
Fullard <- qread(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Fullard_et_al./data/processed_data/Fullard.qs"))
Fullard_mg <- subset(Fullard,subset = celltype == 'Mic')
qsave(Fullard_mg,file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/Fullard_et_al./data/processed_data/Fullard_mg.qs"))

DimPlot(Fullard_mg)
FeaturePlot(Fullard_mg,features = 'PLP1')
