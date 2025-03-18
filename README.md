# Analysis of scRNA-seq data from nuclear transfer and in vitro fertilized Xenopus Laevis embryos #
[![DOI](https://zenodo.org/badge/328938212.svg)](https://doi.org/10.5281/zenodo.14982290)

This repository contains the code associated to the manuscript [Zikmund*, Fiorentino*, et al, Differentiation Success of Reprogrammed Cells is Heterogenous In Vivo and Modulated by Somatic Cell Identity Memory, Stem Cell Reports, 2025](https://doi.org/10.1016/j.stemcr.2025.102447).

Specifically, we provide [Script of Script (SoS)](https://vatlab.github.io/sos-docs/notebook.html) Jupyter notebooks, R scripts and Python Jupyter notebooks to perform:
1. Batch integration using Seurat, robustness analysis for cell clustering and cell clustering using the Louvain algorithm (SoS)
2. Outlier cells identification in cluster 6 (SoS)
3. Cell type composition analysis using a generalized linear model (SoS)
4. Subclustering analysis of the "Mixed States" cluster (SoS)
5. Differential expression analysis between NT and IVF cells globally (R script)
6. Differential expression analysis between NT and IVF cells for each cell cluster (R script)
7. Definition and analysis of ON- and OFF- memory genes at cluster level (R script)
8. RNA velocity analysis using the dynamical model from scvelo (Python Jupyter Notebook)
9. Fate mapping analysis based on the RNA velocity using CellRank (Python Jupyter Notebook)
10. Double bar heatmap of cluster marker genes (SoS)
