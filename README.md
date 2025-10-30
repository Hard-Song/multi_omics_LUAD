# Multi-Omics-LUAD

This code performs single-cell RNA sequencing (scRNA-seq) data analysis for the GSE117570 dataset.
Core functionalities include:
Loading and merging multiple single-cell datasets
Quality control filtering based on feature counts and mitochondrial gene percentage
Data normalization and identification of variable features
Dimensionality reduction (PCA and UMAP) with batch effect correction using Harmony
Cell clustering and visualization of clusters
Identification of marker genes for each cluster
The analysis utilizes the Seurat package for key scRNA-seq processing steps, following standard workflows for data integration, dimensionality reduction, and cell type characterization.
