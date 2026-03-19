# Muon-based snATAC-seq processing pipeline

This directory contains the Muon-based workflow for generating peak-by-cell `.h5ad` files from raw fragment files and performing downstream integration and annotation transfer.

## Directory structure

<pre class="overflow-visible! px-0!" data-start="308" data-end="417"><div class="relative w-full mt-4 mb-1"><div class=""><div class="relative"><div class="h-full min-h-0 min-w-0"><div class="h-full min-h-0 min-w-0"><div class="border border-token-border-light border-radius-3xl corner-superellipse/1.1 rounded-3xl"><div class="h-full w-full border-radius-3xl bg-token-bg-elevated-secondary corner-superellipse/1.1 overflow-clip rounded-3xl lxnfua_clipPathFallback"><div class="pointer-events-none absolute inset-x-4 top-12 bottom-4"><div class="pointer-events-none sticky z-40 shrink-0 z-1!"><div class="sticky bg-token-border-light"></div></div></div><div class=""><div class="relative z-0 flex max-w-full"><div id="code-block-viewer" dir="ltr" class="q9tKkq_viewer cm-editor z-10 light:cm-light dark:cm-light flex h-full w-full flex-col items-stretch ͼd ͼr"><div class="cm-scroller"><div class="cm-content q9tKkq_readonly"><span>02_Muon/</span><br/><span>├── 01_generate_atac_h5ad.py</span><br/><span>├── 02_run_generate_atac_h5ad.sh</span><br/><span>└── 03_muon_integration.py</span></div></div></div></div></div></div></div></div></div><div class=""><div class=""></div></div></div></div></div></pre>

## Overview

This workflow is organized into three steps:

* **Step 1: `01_generate_atac_h5ad.py`**

  Generate per-sample snATAC-seq `.h5ad` files from raw fragment files and a predefined peak set.
* **Step 2: `02_run_generate_atac_h5ad.sh`**

  Submit or run the h5ad generation step on an HPC cluster using a batch script.
* **Step 3: `03_muon_integration.py`**

  Merge sample-level `.h5ad` files, harmonize barcodes, filter cells using ArchR QC output, add matched RNA annotations, and run downstream snATAC integration.

---

### Input data

The workflow assumes the following inputs are available:

* raw fragment files (`*.atac_fragments.tsv.gz`)
* a predefined peak set file
* ArchR-filtered cell metadata table
* matched RNA metadata table

---

# Step 1. Generate per-sample h5ad files

## Script

<pre class="overflow-visible! px-0!" data-start="1595" data-end="1639"><div class="relative w-full mt-4 mb-1"><div class=""><div class="relative"><div class="h-full min-h-0 min-w-0"><div class="h-full min-h-0 min-w-0"><div class="border border-token-border-light border-radius-3xl corner-superellipse/1.1 rounded-3xl"><div class="h-full w-full border-radius-3xl bg-token-bg-elevated-secondary corner-superellipse/1.1 overflow-clip rounded-3xl lxnfua_clipPathFallback"><div class="pointer-events-none absolute inset-x-4 top-12 bottom-4"><div class="pointer-events-none sticky z-40 shrink-0 z-1!"><div class="sticky bg-token-border-light"></div></div></div><div class=""><div class="relative z-0 flex max-w-full"><div id="code-block-viewer" dir="ltr" class="q9tKkq_viewer cm-editor z-10 light:cm-light dark:cm-light flex h-full w-full flex-col items-stretch ͼd ͼr"><div class="cm-scroller"><div class="cm-content q9tKkq_readonly"><span>02_Muon/01_generate_atac_h5ad.py</span></div></div></div></div></div></div></div></div></div><div class=""><div class=""></div></div></div></div></div></pre>

## Purpose

This script creates a per-sample snATAC-seq `.h5ad` object from:

* a fragment file
* a predefined set of peak regions

It counts fragments overlapping peaks and stores the result as an AnnData object for downstream Muon/Scanpy analysis.

## Main procedures

1. Read the fragment file for a given sample
2. Summarize total fragment counts per barcode
3. Initialize an AnnData object with cell metadata
4. Read the peak file
5. Restrict peaks to standard chromosomes
6. Locate fragments and count overlaps in peaks
7. Convert the count matrix to integer format
8. Perform basic cell filtering
9. Write the per-sample `.h5ad` output

## Inputs

* `--indir`: directory containing fragment files
* `--sample`: sample name
* `--peak`: peak file with chromosome, start, and end columns
* `--outdir`: output directory

## Example

<pre class="overflow-visible! px-0!" data-start="2477" data-end="2639"><div class="relative w-full mt-4 mb-1"><div class=""><div class="relative"><div class="h-full min-h-0 min-w-0"><div class="h-full min-h-0 min-w-0"><div class="border border-token-border-light border-radius-3xl corner-superellipse/1.1 rounded-3xl"><div class="h-full w-full border-radius-3xl bg-token-bg-elevated-secondary corner-superellipse/1.1 overflow-clip rounded-3xl lxnfua_clipPathFallback"><div class="pointer-events-none absolute inset-x-4 top-12 bottom-4"><div class="pointer-events-none sticky z-40 shrink-0 z-1!"><div class="sticky bg-token-border-light"></div></div></div><div class=""><div class="relative z-0 flex max-w-full"><div id="code-block-viewer" dir="ltr" class="q9tKkq_viewer cm-editor z-10 light:cm-light dark:cm-light flex h-full w-full flex-col items-stretch ͼd ͼr"><div class="cm-scroller"><div class="cm-content q9tKkq_readonly"><span>python 01_generate_atac_h5ad.py \</span><br/><span></span><span class="ͼn">--indir</span><span> /path/to/fragments \</span><br/><span></span><span class="ͼn">--sample</span><span> Sample1 \</span><br/><span></span><span class="ͼn">--peak</span><span> /path/to/peaks.tsv \</span><br/><span></span><span class="ͼn">--outdir</span><span> /path/to/output</span></div></div></div></div></div></div></div></div></div><div class=""><div class=""></div></div></div></div></div></pre>

## Output

This step produces one `.h5ad` file per sample:

<pre class="overflow-visible! px-0!" data-start="2701" data-end="2735"><div class="relative w-full mt-4 mb-1"><div class=""><div class="relative"><div class="h-full min-h-0 min-w-0"><div class="h-full min-h-0 min-w-0"><div class="border border-token-border-light border-radius-3xl corner-superellipse/1.1 rounded-3xl"><div class="h-full w-full border-radius-3xl bg-token-bg-elevated-secondary corner-superellipse/1.1 overflow-clip rounded-3xl lxnfua_clipPathFallback"><div class="pointer-events-none absolute inset-x-4 top-12 bottom-4"><div class="pointer-events-none sticky z-40 shrink-0 z-1!"><div class="sticky bg-token-border-light"></div></div></div><div class=""><div class="relative z-0 flex max-w-full"><div id="code-block-viewer" dir="ltr" class="q9tKkq_viewer cm-editor z-10 light:cm-light dark:cm-light flex h-full w-full flex-col items-stretch ͼd ͼr"><div class="cm-scroller"><div class="cm-content q9tKkq_readonly"><span><outdir>/<sample>.h5ad</span></div></div></div></div></div></div></div></div></div><div class=""><div class=""></div></div></div></div></div></pre>

# Step 2. Run h5ad generation on cluster

## Script

<pre class="overflow-visible! px-0!" data-start="2983" data-end="3031"><div class="relative w-full mt-4 mb-1"><div class=""><div class="relative"><div class="h-full min-h-0 min-w-0"><div class="h-full min-h-0 min-w-0"><div class="border border-token-border-light border-radius-3xl corner-superellipse/1.1 rounded-3xl"><div class="h-full w-full border-radius-3xl bg-token-bg-elevated-secondary corner-superellipse/1.1 overflow-clip rounded-3xl lxnfua_clipPathFallback"><div class="pointer-events-none absolute inset-x-4 top-12 bottom-4"><div class="pointer-events-none sticky z-40 shrink-0 z-1!"><div class="sticky bg-token-border-light"></div></div></div><div class=""><div class="relative z-0 flex max-w-full"><div id="code-block-viewer" dir="ltr" class="q9tKkq_viewer cm-editor z-10 light:cm-light dark:cm-light flex h-full w-full flex-col items-stretch ͼd ͼr"><div class="cm-scroller"><div class="cm-content q9tKkq_readonly"><span>02_Muon/02_run_generate_atac_h5ad.sh</span></div></div></div></div></div></div></div></div></div><div class=""><div class=""></div></div></div></div></div></pre>

## Purpose

This shell script provides an HPC-friendly wrapper for running `01_generate_atac_h5ad.py` on a compute cluster.

It is intended for generating per-sample `.h5ad` files in batch.

## Main procedures

1. Set SLURM job parameters
2. Activate the appropriate conda environment
3. Define input, output, and peak file paths
4. Run the Python script for a specified sample

# Step 3. Merge and integrate Muon h5ad files

## Script

<pre class="overflow-visible! px-0!" data-start="4194" data-end="4236"><div class="relative w-full mt-4 mb-1"><div class=""><div class="relative"><div class="h-full min-h-0 min-w-0"><div class="h-full min-h-0 min-w-0"><div class="border border-token-border-light border-radius-3xl corner-superellipse/1.1 rounded-3xl"><div class="h-full w-full border-radius-3xl bg-token-bg-elevated-secondary corner-superellipse/1.1 overflow-clip rounded-3xl lxnfua_clipPathFallback"><div class="pointer-events-none absolute inset-x-4 top-12 bottom-4"><div class="pointer-events-none sticky z-40 shrink-0 z-1!"><div class="sticky bg-token-border-light"></div></div></div><div class=""><div class="relative z-0 flex max-w-full"><div id="code-block-viewer" dir="ltr" class="q9tKkq_viewer cm-editor z-10 light:cm-light dark:cm-light flex h-full w-full flex-col items-stretch ͼd ͼr"><div class="cm-scroller"><div class="cm-content q9tKkq_readonly"><span>02_Muon/03_muon_integration.py</span></div></div></div></div></div></div></div></div></div><div class=""><div class=""></div></div></div></div></div></pre>

## Purpose

This script merges the sample-level `.h5ad` files generated in Step 1, filters them using ArchR-derived cell metadata, adds matched RNA annotations, and runs downstream snATAC integration.

## Main procedures

1. Read all `.h5ad` files in the working directory
2. Concatenate them into a single AnnData object
3. Standardize cell barcodes
4. Filter cells using the ArchR filtered cell metadata table
5. Filter to cells overlapping with RNA metadata
6. Add RNA-derived annotations into `adata.obs`
7. Save an intermediate merged `.h5ad`
8. Call `snATAC_integrate.py` for downstream integration

## Output

This step produces:

* an intermediate merged `.h5ad` file
* an integrated `.h5ad` file generated through `snATAC_integrate.py`

---

## Relationship to the ArchR pipeline

This Muon workflow complements the ArchR workflow:

* **ArchR** is used for fragment QC, doublet filtering, initial peak calling, and project-level downstream chromatin analyses
* **Muon/Scanpy** is used here to generate peak-by-cell matrices, refine embeddings, and support downstream integration and annotation transfer
