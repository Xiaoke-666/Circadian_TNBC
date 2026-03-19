# ArchR analysis pipeline for snATAC-seq

This directory contains the ArchR-based preprocessing, integration, and downstream analysis pipeline for snATAC-seq data.

## Directory structure

<pre class="overflow-visible! px-0!" data-start="261" data-end="344"><div class="relative w-full mt-4 mb-1"><div class=""><div class="relative"><div class="h-full min-h-0 min-w-0"><div class="h-full min-h-0 min-w-0"><div class="border border-token-border-light border-radius-3xl corner-superellipse/1.1 rounded-3xl"><div class="h-full w-full border-radius-3xl bg-token-bg-elevated-secondary corner-superellipse/1.1 overflow-clip rounded-3xl lxnfua_clipPathFallback"><div class="pointer-events-none absolute inset-x-4 top-12 bottom-4"><div class="pointer-events-none sticky z-40 shrink-0 z-1!"><div class="sticky bg-token-border-light"></div></div></div><div class=""><div class="relative z-0 flex max-w-full"><div id="code-block-viewer" dir="ltr" class="q9tKkq_viewer cm-editor z-10 light:cm-light dark:cm-light flex h-full w-full flex-col items-stretch ͼd ͼr"><div class="cm-scroller"><div class="cm-content q9tKkq_readonly"><span>01_ArchR/</span><br/><span>├── 01_ArchR_integration.R</span><br/><span>└── 02_ArchR_downstream_analysis.R</span></div></div></div></div></div></div></div></div></div><div class=""><div class=""></div></div></div></div></div></pre>

## Overview

The workflow is organized into two major steps:

* **Step 1: `01_ArchR_integration.R`**

  Preprocess raw snATAC-seq fragment files, generate Arrow files, remove doublets, perform dimensionality reduction and clustering, call an initial peak set, and integrate refined UMAP coordinates and cell annotations from external h5ad-based analysis.
* **Step 2: `02_ArchR_downstream_analysis.R`**

  Perform downstream analyses on the annotated ArchR project, including UMAP visualization, marker gene and marker peak analysis, motif deviation analysis, RNA-ATAC integration, co-accessibility, peak-to-gene linkage, grouped matrix export, and bigWig generation.

---

## Requirements

### R packages

The pipeline uses the following R packages:

<pre class="overflow-visible! px-0!" data-start="1100" data-end="1182"><div class="relative w-full mt-4 mb-1"><div class=""><div class="relative"><div class="h-full min-h-0 min-w-0"><div class="h-full min-h-0 min-w-0"><div class="border border-token-border-light border-radius-3xl corner-superellipse/1.1 rounded-3xl"><div class="h-full w-full border-radius-3xl bg-token-bg-elevated-secondary corner-superellipse/1.1 overflow-clip rounded-3xl lxnfua_clipPathFallback"><div class="pointer-events-none absolute inset-x-4 top-12 bottom-4"><div class="pointer-events-none sticky z-40 shrink-0 z-1!"><div class="sticky bg-token-border-light"></div></div></div><div class=""><div class="relative z-0 flex max-w-full"><div id="code-block-viewer" dir="ltr" class="q9tKkq_viewer cm-editor z-10 light:cm-light dark:cm-light flex h-full w-full flex-col items-stretch ͼd ͼr"><div class="cm-scroller"><div class="cm-content q9tKkq_readonly"><span class="ͼm">ArchR</span><br/><span class="ͼm">GenomicRanges</span><br/><span class="ͼm">stringr</span><br/><span class="ͼm">dplyr</span><br/><span class="ͼm">data.table</span><br/><span class="ͼm">BSgenome.Mmusculus.UCSC.mm10</span></div></div></div></div></div></div></div></div></div><div class=""><div class=""></div></div></div></div></div></pre>

### External tools

* `MACS2` must be installed and available in `$PATH` for peak calling.

### Genome

This workflow is configured for:

* **mouse genome** : `mm10`

---

## Step 1. ArchR preprocessing and integration

### Script

<pre class="overflow-visible! px-0!" data-start="1415" data-end="1458"><div class="relative w-full mt-4 mb-1"><div class=""><div class="relative"><div class="h-full min-h-0 min-w-0"><div class="h-full min-h-0 min-w-0"><div class="border border-token-border-light border-radius-3xl corner-superellipse/1.1 rounded-3xl"><div class="h-full w-full border-radius-3xl bg-token-bg-elevated-secondary corner-superellipse/1.1 overflow-clip rounded-3xl lxnfua_clipPathFallback"><div class="pointer-events-none absolute inset-x-4 top-12 bottom-4"><div class="pointer-events-none sticky z-40 shrink-0 z-1!"><div class="sticky bg-token-border-light"></div></div></div><div class=""><div class="relative z-0 flex max-w-full"><div id="code-block-viewer" dir="ltr" class="q9tKkq_viewer cm-editor z-10 light:cm-light dark:cm-light flex h-full w-full flex-col items-stretch ͼd ͼr"><div class="cm-scroller"><div class="cm-content q9tKkq_readonly"><span>01_ArchR/01_ArchR_integration.R</span></div></div></div></div></div></div></div></div></div><div class=""><div class=""></div></div></div></div></div></pre>

### Purpose

This script performs the initial ArchR preprocessing workflow starting from raw fragment files and produces an annotated ArchR project for downstream analysis.

### Main procedures

1. **Create Arrow files**
   * Read raw fragment files from the input directory
   * Generate Arrow files with QC filtering
   * Add TileMatrix and GeneScoreMatrix
2. **Infer doublets**
   * Compute doublet scores using ArchR
   * Remove putative doublets
3. **Create ArchR project**
   * Build an ArchR project from Arrow files
   * Save intermediate project objects and cell metadata
4. **Dimensionality reduction**
   * Run iterative LSI on the TileMatrix
5. **Clustering and UMAP**
   * Perform clustering
   * Compute UMAP embedding
6. **Initial peak calling**
   * Run `addGroupCoverages`
   * Run `addReproduciblePeakSet`
   * Add PeakMatrix
7. **Integrate external results**
   * Import external UMAP coordinates generated from h5ad/muon analysis
   * Import refined cell-level metadata from external annotation tables
   * Update the ArchR project with final cell type annotations
8. **Save annotated project**
   * Write the final annotated ArchR project for downstream analyses

### Notes

* Peak calling in this step is performed in two stages conceptually:
  * First, an **initial cluster-based peak set** is generated using default ArchR clustering.
  * This initial peak set is then used as the feature space for downstream h5ad/muon-based refinement.
  * Refined UMAP coordinates and cell annotations are imported back into ArchR to create the final annotated project.

### Expected outputs

Typical outputs from this step include:

* Arrow files
* QC summaries
* initial ArchR project
* doublet-filtered ArchR project
* initial cluster-based peak matrix
* annotated ArchR project with imported UMAP and cell annotations

---

## Step 2. Downstream ArchR analysis

### Script

<pre class="overflow-visible! px-0!" data-start="3355" data-end="3406"><div class="relative w-full mt-4 mb-1"><div class=""><div class="relative"><div class="h-full min-h-0 min-w-0"><div class="h-full min-h-0 min-w-0"><div class="border border-token-border-light border-radius-3xl corner-superellipse/1.1 rounded-3xl"><div class="h-full w-full border-radius-3xl bg-token-bg-elevated-secondary corner-superellipse/1.1 overflow-clip rounded-3xl lxnfua_clipPathFallback"><div class="pointer-events-none absolute inset-x-4 top-12 bottom-4"><div class="pointer-events-none sticky z-40 shrink-0 z-1!"><div class="sticky bg-token-border-light"></div></div></div><div class=""><div class="relative z-0 flex max-w-full"><div id="code-block-viewer" dir="ltr" class="q9tKkq_viewer cm-editor z-10 light:cm-light dark:cm-light flex h-full w-full flex-col items-stretch ͼd ͼr"><div class="cm-scroller"><div class="cm-content q9tKkq_readonly"><span>01_ArchR/02_ArchR_downstream_analysis.R</span></div></div></div></div></div></div></div></div></div><div class=""><div class=""></div></div></div></div></div></pre>

### Purpose

This script performs downstream analyses using the annotated ArchR project generated in Step 1.

### Main procedures

1. **Load annotated ArchR project**
   * Load the final project created in Step 1
2. **UMAP visualization**
   * Plot UMAP embeddings colored by annotated ATAC cell types
3. **Marker gene analysis**
   * Compute marker genes from the GeneScoreMatrix
   * Export marker gene tables by cell type
4. **Peak calling by annotated cell type**
   * Recompute group coverages using finalized ATAC cell type annotations
   * Call reproducible peaks for each annotated cell type
   * Add PeakMatrix
   * Export peak sets in TSV and BED format
5. **Marker peak analysis**
   * Identify marker peaks from the PeakMatrix
6. **Motif analysis**
   * Add motif annotations
   * Compute chromVAR deviation scores
   * Export group-level motif deviation matrices
7. **RNA-ATAC integration**
   * Load matched RNA reference
   * Perform cross-platform integration using `addGeneIntegrationMatrix`
   * Visualize predicted group assignments on UMAP
8. **Co-accessibility analysis**
   * Compute co-accessible peak pairs
   * Export co-accessibility results
9. **Peak-to-gene linkage**
   * Compute peak-to-gene links
   * Export linked peak-gene pairs in TSV and BED format
10. **Grouped matrix export**

* Export grouped GeneScoreMatrix by:
  * ATAC cell type
  * cell type × time point
* Export grouped PeakMatrix by:
  * ATAC cell type
  * cell type × time point

11. **bigWig export**

* Generate group-level bigWig tracks for visualization

### Notes

* In this step, peak calling is repeated using  **final annotated ATAC cell types** , producing a refined and biologically interpretable peak set for downstream analyses.
* This script assumes that external RNA reference data and annotation files are already prepared and available.

### Expected outputs

Typical outputs from this step include:

* UMAP PDFs
* marker gene tables
* reproducible peak sets
* marker peak results
* chromVAR deviation matrices
* RNA-ATAC integration outputs
* co-accessibility results
* peak-to-gene links
* grouped gene score and peak matrices
* bigWig tracks

---

## Output summary

After both steps, the workflow produces:

* filtered and annotated ArchR projects
* refined peak sets
* marker genes and marker peaks
* motif deviation results
* RNA-ATAC integration results
* co-accessibility and peak-to-gene linkage tables
* grouped accessibility and gene score matrices
* bigWig tracks for genome browser visualization
