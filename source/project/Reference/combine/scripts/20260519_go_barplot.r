library(readxl)
library(ggplot2)
library(dplyr)

go_data <- read_excel(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/scripts/BULK/data/saved_data/20260519_go.xlsx"))

p <- ggplot(go_data, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
  geom_col(width = 0.65) +
  scale_fill_gradient(low = "#2E86AB", high = "#A2BCE0") +
  labs(title = "GO Enrichment Bar Plot", x = "Gene Count", y = "GO Term", fill = "Adjusted p-value") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10), panel.grid = element_blank())

ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/results/figure/20260519_go_barplot.pdf"), plot = p, width = 8, height = 4, dpi = 300)
ggsave(file.path(Sys.getenv("MG_PROJECT_ROOT"), "Reference/combine/results/figure/20260519_go_barplot.png"), plot = p, width = 8, height = 4, dpi = 300)
