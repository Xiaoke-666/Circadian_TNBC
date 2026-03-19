#!/usr/bin/env python3

import os
import sys
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
from muon import atac as ac



# -----------------------------
# Global settings
# -----------------------------
WORKDIR = "/path/to/working_directory"
ARCHR_CELL_FILE = "/path/to/archr_cell_metadata.tsv"
RNA_OBS_FILE = "/path/to/rna_obs.tsv.gz"

# intermediate output before integration
OUTPUT_FILE = "snATAC1.h5ad"

# integration settings
INTEGRATE_SCRIPT = Path(__file__).resolve().parent / "snATAC_integrate.py"
INTEGRATE_PREFIX = "snATAC1.noCorr"
BATCH_COLUMN = "Sample"

os.chdir(WORKDIR)

sc.figure_size = (8, 8)
sc.set_figure_params(
    scanpy=True,
    dpi=300,
    dpi_save=300,
    fontsize=4,
    color_map="viridis",
)
pd.options.display.max_rows = 999


# -----------------------------
# Plotting helpers
# -----------------------------
def pl_umap(adata, columns, prefix):
    for col in columns:
        sc.pl.umap(
            adata,
            color=col,
            legend_loc="on data",
            legend_fontsize=4,
            save=f".{prefix}.{col}.ondata.png",
        )
        sc.pl.umap(
            adata,
            color=col,
            save=f".{prefix}.{col}.png",
        )

# -----------------------------
# Read and merge h5ad files
# -----------------------------
def read_and_merge_h5ad_files(workdir):
    filenames = sorted([f for f in os.listdir(workdir) if f.endswith(".h5ad")])

    if len(filenames) == 0:
        raise FileNotFoundError("No .h5ad files were found in the working directory.")

    print(f"Found {len(filenames)} h5ad files")

    adatas = []
    for filename in filenames:
        print(f"Reading {filename}")
        adata = sc.read(filename)
        # adata.var_names_make_unique()
        adatas.append(adata)

    filenames_clean = [os.path.splitext(filename)[0] for filename in filenames]

    if len(adatas) == 1:
        merged = adatas[0].copy()
    else:
        merged = adatas[0].concatenate(
            adatas[1:],
            batch_categories=filenames_clean,
        )

    return merged


# -----------------------------
# Standardize barcodes
# -----------------------------
def standardize_barcodes(adata):
    adata = adata.copy()
    adata.obs["Sample_barcode"] = adata.obs.index.str.replace(r"-1-.*", "-1", regex=True)
    adata.obs.index = adata.obs["Sample_barcode"]
    adata.obs.index.name = None
    return adata


# -----------------------------
# Filter by ArchR cells
# -----------------------------
def filter_by_archr_cells(adata, archr_cell_file):
    archr_cells = pd.read_csv(
        archr_cell_file,
        sep="\t",
        header=0,
        index_col=0,
    )

    filtered = adata[adata.obs.index.isin(archr_cells.index), :].copy()
    print(f"After ArchR cell filtering: {filtered.n_obs} cells")
    return filtered


# -----------------------------
# Add matched RNA annotations
# -----------------------------
def add_rna_annotations(adata, rna_obs_file):
    rna = pd.read_csv(
        rna_obs_file,
        sep="\t",
        header=0,
        index_col=0,
    )

    adata = adata[adata.obs.index.isin(rna.index), :].copy()
    print(f"After RNA overlap filtering: {adata.n_obs} cells")

    for col in rna.columns:
        if col not in adata.obs.columns:
            adata.obs[col] = adata.obs.index.map(rna[col])

    return adata


# -----------------------------
# Save output
# -----------------------------
def save_h5ad(adata, output_file):
    print(f"Writing intermediate h5ad: {output_file}")
    adata.write_h5ad(output_file)


# -----------------------------
# Run snATAC integration
# -----------------------------
def run_integration(input_file, integrate_script, prefix, batch_column, output_dir="."):
    integrate_script = Path(integrate_script)

    if not integrate_script.exists():
        raise FileNotFoundError(f"Integration script not found: {integrate_script}")

    cmd = [
        sys.executable,
        str(integrate_script),
        "--input", input_file,
        "--output", output_dir,
        "--prefix", prefix,
        "--batch", batch_column,
    ]

    print("Running integration command:")
    print(" ".join(cmd))

    subprocess.run(cmd, check=True)


# -----------------------------
# Main
# -----------------------------
def main():
    print("Step 1: Read and merge h5ad files")
    adata = read_and_merge_h5ad_files(WORKDIR)

    print("Step 2: Standardize cell barcodes")
    adata = standardize_barcodes(adata)

    print("Step 3: Filter by ArchR cells")
    atac1 = filter_by_archr_cells(adata, ARCHR_CELL_FILE)

    print("Step 4: Add matched RNA annotations")
    atac1 = add_rna_annotations(atac1, RNA_OBS_FILE)

    print("Step 5: Save intermediate h5ad")
    save_h5ad(atac1, OUTPUT_FILE)

    print("Step 6: Run snATAC integration")
    run_integration(
        input_file=OUTPUT_FILE,
        integrate_script=INTEGRATE_SCRIPT,
        prefix=INTEGRATE_PREFIX,
        batch_column=BATCH_COLUMN,
        output_dir="./",
    )

    print("Pipeline completed successfully.")


if __name__ == "__main__":
    main()