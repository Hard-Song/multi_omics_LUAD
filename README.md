# Multi-Omics-LUAD
This code implements a comprehensive single-cell RNA sequencing (scRNA-seq) data analysis pipeline for processing the GSE117570 dataset. It follows established workflows for scRNA-seq data processing to facilitate the identification and characterization of cell populations.

## Core Workflow

The analysis proceeds through several key stages:
- Data Loading & Merging: The code loads multiple filtered scRNA-seq datasets from the GSE117570 cohort and merges them into a single integrated dataset, preserving sample-specific identifiers for downstream analysis.

- Quality Control: Quality metrics (including the number of detected features, total RNA counts, and percentage of mitochondrial genes) are calculated and used to filter low-quality cells, ensuring only cells with robust profiles (500–6000 features and <10% mitochondrial content) are retained.

- Normalization & Feature Selection: The filtered data undergoes normalization to account for technical variability, followed by identification of 2000 highly variable features to focus downstream analyses on biologically informative genes.

- Dimensionality Reduction & Batch Correction: Principal Component Analysis (PCA) is performed to reduce data dimensionality, with Harmony applied to correct for batch effects across samples. This is followed by UMAP projection to visualize cellular relationships in two dimensions.

- Clustering & Marker Identification: Cells are grouped into clusters using shared nearest neighbor algorithms, and cluster-specific marker genes are identified to characterize the biological identity of each cell population.

This pipeline leverages the Seurat package for core scRNA-seq processing, with visualization and statistical steps to support the interpretation of cellular heterogeneity within the dataset.
