---
layout: single
classes:
  - landing
author: Hugo Chenel
title: "Preprocessing scRNA-seq data with Seurat"
date: 2024-11-27
categories: []
tags: [Science & Tech]
excerpt: "Tutorial - Single cell series 1.1"
permalink: /blog/preprocessing-scrna-seq-data/
header:
  overlay_color: "#93b874"
  overlay_filter: "0.005"
show_date: true
pagination: false
show_taxonomy: false
read_time: true
toc: true
toc_sticky: true
toc_label: "Table of contents"
toc_icon: "cog"
sidebar:
  nav: "sc_seurat"
---

# Preprocessing scRNA-seq data with Seurat

**Author:** Hugo Chenel  
**Purpose:** This advanced tutorial guides researchers through preprocessing single-cell RNA-seq (scRNA-seq) data using Seurat, a powerful R package for single-cell analysis. Topics include quality control, normalization, clustering, and multimodal (RNA + ADT) analysis.
<br>Talk about single cell and CITEseq, CLL and my data.

---

## Prerequisites
Make sure the required R packages are installed and loaded:

```r
install.packages(c("Seurat", "ggplot2", "patchwork", "cowplot", "dplyr", "ComplexHeatmap"))
BiocManager::install(c("SingleR", "celldex"))
```

```r
# Load libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
library(SingleR)
library(celldex)
library(ComplexHeatmap)
```
## Load CITE-seq data
CITE-seq data contains RNA and protein expression profiles. We will load the raw count matrices stored in HDF5 files using ```Read10X_h5```. Each dataset corresponds to a different experimental condition or time point.

```r
# Load raw count matrices
cll.filtered.matrix_J1 <- Read10X_h5("path_to_file_J1.h5")
cll.filtered.matrix_J4 <- Read10X_h5("path_to_file_J4.h5")
cll.filtered.matrix_J8 <- Read10X_h5("path_to_file_J8.h5")
cll.filtered.matrix_J11 <- Read10X_h5("path_to_file_J11.h5")
cll.filtered.matrix_J14 <- Read10X_h5("path_to_file_J14.h5")
```

## Create Seurat objects for RNA and ADT data
We create Seurat objects for RNA data (```Gene Expression```) and add protein data (```Antibody Capture```) as separate assays. The datasets are then merged into a single Seurat object for downstream analyses.

```r
# Create Seurat objects for RNA
cll.J1.rna <- CreateSeuratObject(counts = cll.filtered.matrix_J1$`Gene Expression`, project = "J1")
cll.J4.rna <- CreateSeuratObject(counts = cll.filtered.matrix_J4$`Gene Expression`, project = "J4")
# Repeat for other time points...

# Add protein data (ADT)
cll.J1.rna[["ADT"]] <- CreateAssayObject(counts = cll.filtered.matrix_J1$`Antibody Capture`)
# Repeat for other time points...

# Merge datasets
cll.combined <- merge(cll.J1.rna, y = c(cll.J4.rna, cll.J8.rna, cll.J11.rna, cll.J14.rna), 
                      add.cell.ids = c("J1", "J4", "J8", "J11", "J14"), project = "CLL_Combined")
```

## Perform quality control and normalization
Quality control removes low-quality cells (e.g., those with excessive mitochondrial content). Normalization adjusts for sequencing depth, ensuring comparability across cells. Identifying variable features highlights the most informative genes.

### Quality control

```r
# Add mitochondrial gene percentage
cll.combined[["percent.mt"]] <- PercentageFeatureSet(cll.combined, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(cll.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells
cll.combined <- subset(cll.combined, subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 5)
```
### Normalization and variable features

```r
# Normalize data
cll.combined <- NormalizeData(cll.combined)

# Identify highly variable genes
cll.combined <- FindVariableFeatures(cll.combined, selection.method = "vst", nfeatures = 2000)
```
🖼️ Add violin plots and scatter plots here



















## Session information
SessionInfo()
