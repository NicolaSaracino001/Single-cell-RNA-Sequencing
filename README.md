# Human Pancreatic Islets Single-Cell RNA-Seq Analysis

![R](https://img.shields.io/badge/R-4.x-blue.svg)
![Seurat](https://img.shields.io/badge/Seurat-v5-green.svg)
![Bioinformatics](https://img.shields.io/badge/Domain-Bioinformatics-orange.svg)

## 📌 Project Overview
This project presents an end-to-end single-cell RNA sequencing (scRNA-seq) analysis pipeline executed on human pancreatic islet dataset using R and the **Seurat** framework. 

The goal of this study is to profile cell heterogeneity within human pancreatic islets, separate endocrinal and non-endocrinal cell types, identify cell-type-specific marker genes, and assign biological identities to distinct single-cell clusters.

---

## 🛠️ Computational Pipeline

The analysis was structured into 7 distinct phases to ensure methodological rigor and reproducibility:

1. **Data Loading & Preprocessing**: Importing sparse gene expression matrices into Seurat objects.
2. **Quality Control (QC) & Filtering**: Removal of low-quality cells, doublets, and damaged cells using mitochondrial gene expression thresholds (`percent.mt`) and UMI counts (`nFeature_RNA`, `nCount_RNA`).
3. **Normalization & Feature Selection**: Data normalization using `LogNormalize` and identification of top highly variable features (HVGs).
4. **Dimensionality Reduction (PCA)**: Principal Component Analysis, evaluation of cell-cycle phase influences (G1, S, G2M), and determination of significant PCs using Elbow Plotting.
5. **Clustering & Parameter Optimization**: Benchmarking multiple parameter settings (15 PCs at Res 0.5 vs. 20 PCs at Res 0.8) to prevent over-clustering.
6. **Marker Gene Discovery**: Identification of differentially expressed genes (DEGs) per cluster using non-parametric Wilcoxon rank-sum tests (`FindAllMarkers`).
7. **Cell Type Annotation**: Mapping identified markers to established pancreatic lineage markers to annotate clusters.

---

## 📊 Key Results & Visualizations

### 1. Quality Control & Filtering
Filtering out damaged cells and high-mitochondrial noise ensured clean transcriptomic profiles downstream.

![Distribution of nFeature_RNA, nCount_RNA, percent.mt] (qc_violin_plots.png) 
> - `qc_scatter_plots.png` (Correlation between UMI counts and mitochondrial percentage)

### 2. PCA & Dimension Selection
The Elbow Plot indicated that the majority of biological variation is captured within the first 15 Principal Components.
> - `elbow_plot.png` (Standard Deviation vs. PC numbers)

### 3. Clustering Benchmark
Comparing two parameter settings demonstrated that **15 PCs with Resolution 0.5** provided a biologically sound balance compared to over-clustered alternatives:
> - `umap_comparison_15vs20.png` (Side-by-side UMAP comparison)

### 4. Marker Gene Expression & Profiling
Cluster marker specificity was confirmed using multi-gene heatmaps and dot plots:
> - `heatmap_top10_markers.png` (Top 10 differentially expressed genes per cluster)
> - `dotplot_markers.png` (Expression intensity and percentage per cluster)

### 5. Biological Validation with Canonical Markers
Canonical markers were evaluated across all clusters:
- **`INS`** (Beta cells)
- **`GCG`** (Alpha cells)
- **`SST`** (Delta cells)
- **`KRT19`** (Ductal cells)
- **`PECAM1`** (Endothelial cells)
- **`COL1A1`** (Stellate / Stromal cells)

> - `featureplot_canonical_markers.png` (UMAP expression overlay of canonical markers)
> - `violinplot_ins.png` (INS expression across clusters)

### 6. Final Cell Type Annotation
The final UMAP map delineates all major endocrinal and stromal cell populations present in the tissue sample:
> - `final_annotated_umap.png` (Final annotated UMAP cluster map)

---

## 💻 Dependencies & Requirements

To run this pipeline, install the following R packages:

```r
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))