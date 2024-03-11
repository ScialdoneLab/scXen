# scXen
## Analysis of scRNA-seq data from nuclear transfer Xenopus Laevis embryos ##

This repository contains the code associated to the manuscript [Zikmund*, Fiorentino*, et al, TITLE, biorxiv, 2024](...).

Specifically, we provide [Script of Script (SoS)](https://vatlab.github.io/sos-docs/notebook.html) Jupyter notebooks, R scripts and Python Jupyter notebooks to perform:
1. Batch integration using Seurat, robustness analysis for cell clustering and cell clustering using the Louvain algorithm (SoS)
2. Outlier cells identification in cluster 6 (SoS)
3. Cell type composition analysis (SoS)
4. Computation of specific cluster markers (SoS)
5. Differential expression analysis between NT and IVF cells globally (R script)
6. Differential expression analysis between NT and IVF cells for each cluster (R script)
7. Definition and analysis of ON- and OFF- memory genes at cluster level (R script)
8. RNA velocity analysis (Python Jupyter Notebook)
9. CellRank analysis (Python Jupyter Notebook)
10. Subclustering analysis of the "Mixed States" cluster (SoS)

The /data_utils/ folder contains files needed to run the codes.

The UMI count matrices can be downloaded from the Gene Expression Omnibus (GEO) with accession number XXX.

If you use our code and/or data, please cite:
ADD CITATION
