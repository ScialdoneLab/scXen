# scXen
## Analysis of scRNA-seq data from nuclear transfer Xenopus embryos ##

This repository contains SoS Jupyter notebooks to perform:
1. Batch integration using Seurat, robustness analysis for clustering and cell clustering using the Louvain algorithm
2. Outlier identification in cluster 6
3. Cell type composition analysis (using chi square test or a generalised linear model)
4. Computation of specific cluster markers
5. Deifnition and analysis of ON- and OFF- memory genes at cell-type level
6. Comparison with bulk RNA-seq data
7. Definition of a memory index and computation of the entropy of batch mixing (this has not been done yet)
8. Computation of UMAPs per experiment
9. RNA velocity analysis for each batch (this is a Python3 Jupyter notebook)

The data (both the raw counts oer batch or the RDS object with the processed data) can be downloaded from the IES server (folder /nas_storage/jonathan.fiorentino/scXen_data/). The folder contains also the bulk RNA-seq data and other lists of genes used in the analyses.
