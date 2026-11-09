# Validation-only copy of:
#   combine/scripts/CellChat/scripts/major_celltype.r (Yao section, lines 31-50)
# Original files and objects are never overwritten.

library(qs)
library(dplyr)
library(magrittr)

project_root <- file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference")
validation_root <- file.path(project_root, "combine/results/stroke_cellchat_script_validation_20260917")
dir.create(validation_root, recursive = TRUE, showWarnings = FALSE)

yao_input <- file.path(project_root, "Yao_et_al./data/processed_data/Yao_v2.qs")
yao_output <- file.path(validation_root, "Yao_v2_with_combine_celltype_reproduced.qs")
stroke_output <- file.path(validation_root, "Yao_stroke_reproduced.qs")

# The following label normalization and sample selection are unchanged from the
# Yao block in the original script. The unrelated mouse block is omitted.
yao <- qread(yao_input)
yao$combine_celltype <- as.character(yao$CellType)
yao@meta.data %<>%
  mutate(combine_celltype = case_when(
    combine_celltype %in% c("Vasular") ~ "Vascular",
    combine_celltype %in% c("Micro") ~ "Microglia",
    combine_celltype %in% c("Neurons") ~ "Neuron",
    combine_celltype %in% c("Astro", "Unknown1") ~ "Astrocyte",
    combine_celltype %in% c("Oli") ~ "OLIG",
    combine_celltype %in% c("Unknown2") ~ "Ependymal",
    combine_celltype %in% c("T cells") ~ "T cell",
    TRUE ~ combine_celltype
  ))

my_order <- c("Microglia", "Neuron", "Astrocyte", "OLIG", "OPC",
              "Ependymal", "Vascular", "T cell")
yao$combine_celltype <- factor(yao$combine_celltype, levels = my_order)
qsave(yao, yao_output)

yao_stroke <- subset(yao, subset = Sample %in% c("l_48y"))
qsave(yao_stroke, stroke_output)

cat("Reproduced Yao Stroke object:\n", stroke_output, "\n", sep = "")
print(table(yao_stroke$combine_celltype, useNA = "ifany"))
