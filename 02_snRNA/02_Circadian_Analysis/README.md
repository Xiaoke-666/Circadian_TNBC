# Circadian Analysis (Single-cell RNA-seq)

This repository provides a reproducible pipeline for identifying and analyzing circadian gene programs from single-cell RNA-seq data.

The workflow includes pseudobulk generation, circadian rhythm detection using DiscoRhythm, repeated replicate-based analysis, visualization, and circadian gene module discovery.

---

## Overview

This pipeline consists of **five major steps**:

### Step 1. Pseudobulk generation

Generate pseudobulk expression matrices by aggregating cells.

The pipeline performs the following steps:

1. Load a processed `AnnData` object (`.h5ad`)
2. Create combined group labels (`Subclass × TimePoint`)
3. Assign pseudo-replicates within each group
4. Aggregate counts (pseudobulk)
5. Normalize and log-transform expression
6. Export expression matrix and metadata

### Step 2. Circadian gene detection (DiscoRhythm)

Run DiscoRhythm on pseudobulk expression data to identify rhythmic genes.

This step includes:

* metadata formatting
* gene filtering
* outlier detection (correlation + PCA)
* replicate consistency filtering (ANOVA)
* period detection (24h)
* oscillation detection (cosinor model)

### Step 3. Repeated DiscoRhythm analysis (robustness)

Run circadian analysis across multiple pseudobulk replicates (default: 10 runs).

For each replicate:

* subset metadata per cell type
* run DiscoRhythm independently
* store results per replicate and cell type

### Step 4. Visualization of circadian genes

Visualize circadian gene expression patterns for each cell type.

Features:

* spline smoothing
* z-score normalization
* heatmaps ordered by acrophase
* 8-phase circadian annotation

### Step 5. Circadian gene module analysis

Identify higher-order temporal structure of circadian genes.

Includes:

* clustering circadian genes into modules
* UMAP visualization
* module correlation analysis
* cross-cell-type comparison
