# OLM_microglia

Integrative single-nucleus analysis identifies an oligodendrocyte-like microglial state associated with myelin repair

This repository provides the analysis code for the figures in the manuscript "Integrative single-nucleus analysis identifies an oligodendrocyte-like microglial state associated with myelin repair".

Microglia are established regulators of myelin homeostasis through myelin debris clearance and inflammatory microenvironment modulation, yet whether additional specialized states contribute to myelin formation or maintenance remains unknown. Here, single-cell transcriptomics of human brains across six neurological conditions and multiple brain regions identified a previously unrecognized microglial state—oligodendrocyte-like microglia (OLM)—that expresses oligodendrocyte signature genes (MBP, MOG, ST18) and the OLM-enriched markers TRIM2 and RNF220 while retaining microglial identity. OLM were also identified in macaque and mouse brains, indicating evolutionary conservation, and exhibited a distinctive needle-like morphology. Spatial transcriptomics and ligand-receptor modeling further revealed that OLM closely interact and colocalize with oligodendrocytes. Functionally, myelin debris was sufficient to drive HMC3 cells toward an OLM-like state in vitro, and transplantation of these cells into the prematurely aged SAMP8 mouse brain recruited oligodendrocytes and rescued demyelination. These findings uncover a transdifferentiation trajectory from microglia to an oligodendrocyte-like subtype and suggest that OLM may represent an endogenous cellular resource for myelin repair.

The figure-ordered code is in `analysis/`: `Figure_1.R` to `Figure_6.R` and `Supplementary_Figure_1.R` to `Supplementary_Figure_8.R`. Small analysis tables are in `data/`, and individual result plots are in `figure/pdf/` and `figure/png/`.

To view the code for a panel in R:

```r
source("analysis/Figure_1.R")
mg_list_panels()
mg_show_panel("F1A")
```

Each chapter contains the panel code. Scientific sections use their corresponding input objects and R packages; input-location settings are listed in `config/`.
