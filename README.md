# 10X Genomics Data Integration and Visualization Shiny App

This Shiny app allows users to analyze and visualize single-cell RNA sequencing data from 10X Genomics. It supports data from two datasets (Control and Treated), performs quality control (QC) checks, integrates the datasets, and visualizes the results using various plots, including violin plots, UMAP plots, clustering results, and a cell cycle heatmap.

## **Features**

- **QC Violin Plots:** Visualizes key QC metrics, such as the number of features, total count, mitochondrial percentage, and ribosomal protein percentage for the control and treated datasets.
- **UMAP Before and After Integration:** Displays UMAP visualizations of the data before and after dataset integration.
- **Cluster Distribution:** Visualizes clusters of cells after integration and clustering.
- **Cell Cycle Heatmap:** Plots a heatmap of cell cycle gene expression across the cells, grouped by their sample identity.

## **Requirements**

- R (version >= 4.0.0)
- Shiny package
- Seurat package (version >= 4.0.0)
- ggplot2 package
- patchwork package

### **Install Required Packages**
You can install the required R packages using the following commands:

```r
install.packages("shiny")
install.packages("ggplot2")
install.packages("patchwork")
install.packages("Seurat")
