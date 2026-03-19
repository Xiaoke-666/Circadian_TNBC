#!/usr/bin/env python3

import os
import sys
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


# -----------------------------
# Global settings
# -----------------------------
WORKDIR = "./Input"
OUTPUT_DIR = "."
ATAC_FILE = "./snATAC7.h5ad"

# external integration script
SNRNA_INTEGRATE_SCRIPT = Path(__file__).resolve().parent / "snRNA_integrate.py"

# sample annotation dictionary
SAMPLE_DICT = {
    "Pool141_4": "ZT4.CD45-.rep1",
    "Pool141_5": "ZT4.CD45-.rep2",
    "Pool141_6": "ZT10.CD45-.rep1",
    "Pool141_7": "ZT10.CD45-.rep2",
    "Pool141_8": "ZT16.CD45-.rep1",
    "Pool141_9": "ZT16.CD45-.rep2",
    "Pool141_10": "ZT22.CD45-.rep1",
    "Pool141_11": "ZT22.CD45-.rep2",
    "Pool141_12": "ZT4.CD45+.rep1",
    "Pool141_13": "ZT4.CD45+.rep2",
    "Pool141_14": "ZT10.CD45+.rep1",
    "Pool141_15": "ZT10.CD45+.rep2",
    "Pool141_16": "ZT16.CD45+.rep1",
    "Pool141_17": "ZT16.CD45+.rep2",
    "Pool141_18": "ZT22.CD45+.rep1",
    "Pool141_19": "ZT22.CD45+.rep2",
}

MAJOR_DISC = {
    "T-cells": ["0", "5", "8", "9"],
    "B-cells": ["4", "12"],
    "Epithelial": ["1", "3", "6", "7", "16"],
    "Mesenchymal": ["10", "18"],
    "Myeloid": ["2", "11", "13", "14", "15"],
    "Plasmablasts": ["17"],
}

sc.figure_size = (8, 8)
sc.set_figure_params(
    scanpy=True,
    dpi=300,
    dpi_save=300,
    fontsize=8,
    color_map="viridis",
)


# -----------------------------
# Plotting helper
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
# Read and merge sample h5ad files
# -----------------------------
def read_and_merge_h5ad_files(workdir):
    filenames = sorted(
        [f for f in os.listdir(workdir) if f.endswith(".raw.cellbender_filtered.h5ad")]
    )

    if len(filenames) == 0:
        raise FileNotFoundError("No .raw.cellbender_filtered.h5ad files found.")

    print(f"Found {len(filenames)} files")

    adatas = []
    for filename in filenames:
        print(f"Reading {filename}")
        adata = sc.read(os.path.join(workdir, filename))
        adata.var_names_make_unique()

        sample_name = filename.replace(".raw.cellbender_filtered.h5ad", "")
        adata.obs["Sample"] = sample_name
        adata.obs.index = adata.obs["Sample"] + "#" + adata.obs.index.astype(str)

        adatas.append(adata)

    filenames_clean = [f.replace(".raw.cellbender_filtered.h5ad", "") for f in filenames]

    if len(adatas) == 1:
        merged = adatas[0].copy()
    else:
        merged = adatas[0].concatenate(adatas[1:], batch_categories=filenames_clean)

    return merged


# -----------------------------
# Annotate metadata
# -----------------------------
def annotate_samples(adata):
    adata = adata.copy()

    adata.obs["Sample_barcode"] = adata.obs.index.str.replace(r"-1-.*", "-1", regex=True)
    adata.obs.index = adata.obs["Sample_barcode"]
    adata.obs.index.name = None

    adata.obs["Sample"] = adata.obs["Sample"].map(SAMPLE_DICT)
    adata.obs["TimePoint"] = adata.obs["Sample"].str.split(".").str[0]
    adata.obs["CD45"] = adata.obs["Sample"].str.split(".").str[1]
    adata.obs["Replicate"] = adata.obs["Sample"].str.split(".").str[2]

    # keep first 3 columns in var as in original script
    adata.var = adata.var.iloc[:, 0:3].copy()
    adata.var.columns = ["gene_ids", "feature_types", "genome"]

    return adata


# -----------------------------
# QC and doublet detection
# -----------------------------
def run_qc_and_filter(adata):
    rna = adata.copy()

    rna.var["mt"] = rna.var_names.str.startswith("mt-")
    rna.var["ribo"] = rna.var_names.str.startswith(("RPS", "RPL", "Rps", "Rpl"))

    sc.pp.calculate_qc_metrics(
        rna,
        qc_vars=["mt", "ribo"],
        inplace=True,
        percent_top=[],
        log1p=True,
    )

    # scrublet
    sc.external.pp.scrublet(rna, expected_doublet_rate=0.05, batch_key="batch")
    rna.write_h5ad("rna1.raw.h5ad", compression="gzip", compression_opts=6)

    # basic filtering
    sc.pp.filter_cells(rna, min_genes=200)
    rna2 = rna[rna.obs.pct_counts_mt < 10, :].copy()
    rna2 = rna2[rna2.obs["predicted_doublet"] == False, :].copy()

    rna2.write_h5ad("rna2.h5ad", compression="gzip", compression_opts=6)
    return rna2


# -----------------------------
# Remove doublet-enriched clusters
# -----------------------------
def subset_clean_clusters(rna2):
    rna3 = rna2[
        (
            rna2.obs["leiden_clust_15"].isin(["1", "2", "5", "10", "14", "17", "15,1", "15,2", "15,5"])
            & rna2.obs["CD45"].isin(["CD45-"])
        )
        | (
            rna2.obs["leiden_clust_15"].isin(["0", "3", "4", "6", "7", "8", "9", "11", "12", "13", "16", "15,0", "15,3", "15,4"])
            & rna2.obs["CD45"].isin(["CD45+"])
        ),
        :
    ].copy()

    rna3.write_h5ad("rna3.h5ad", compression="gzip", compression_opts=6)
    return rna3


# -----------------------------
# Annotate major cell types
# -----------------------------
def annotate_major_celltypes(rna3):
    rna3 = rna3.copy()
    rna3.obs["CellType"] = "Unknown"

    for celltype, clusters in MAJOR_DISC.items():
        rna3.obs["CellType"] = np.where(
            rna3.obs["leiden_clust"].isin(clusters),
            celltype,
            rna3.obs["CellType"],
        )

    pl_umap(rna3, ["CellType"], "rna3.Cellbender.noCorr")

    rna3.write_h5ad("rna3.Cellbender.noCorr.h5ad", compression="gzip", compression_opts=6)

    for major in rna3.obs["CellType"].unique():
        print(f"Writing subset: {major}")
        tmp = rna3[rna3.obs["CellType"] == major, :].copy()
        tmp.write_h5ad(f"rna3.Cellbender.{major}.h5ad", compression="gzip", compression_opts=6)

    return rna3


# -----------------------------
# Marker gene analysis
# -----------------------------
def run_marker_analysis(rna3, groupby="CellCellType"):
    rna3 = rna3.copy()

    rna3.layers["raw"] = rna3.X.copy()
    sc.pp.normalize_total(rna3, target_sum=1e4)
    sc.pp.log1p(rna3)
    sc.tl.rank_genes_groups(rna3, groupby, method="wilcoxon")

    result = rna3.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    marker_tables = []
    for group in groups:
        print(f"Exporting markers for {group}")
        tmp = pd.DataFrame(
            {
                "names": result["names"][group],
                "logfoldchanges": result["logfoldchanges"][group],
                "pvals_adj": result["pvals_adj"][group],
            }
        )
        tmp["leiden"] = group
        tmp = tmp[(tmp["pvals_adj"] < 0.05) & (tmp["logfoldchanges"] > 0.5)].copy()

        tmp.to_csv(
            f"rna3.CellType.{group}.padj005.logfc05.csv",
            index=False,
        )
        marker_tables.append(tmp)

    if marker_tables:
        marker_df = pd.concat(marker_tables, axis=0, ignore_index=True)
        marker_df.to_csv("rna3.all_markers.csv", index=False)

    rna3.write_h5ad("rna3.h5ad", compression="gzip", compression_opts=6)
    return rna3


# -----------------------------
# Match RNA and ATAC cells
# -----------------------------
def generate_rna_atac_overlap(rna_file, atac_file, output_file):
    rna3 = sc.read(rna_file)
    atac = sc.read(atac_file)

    intersected_cells = np.intersect1d(rna3.obs.index.tolist(), atac.obs.index.tolist()).tolist()
    rna5 = rna3[rna3.obs.index.isin(intersected_cells), :].copy()

    print(f"RNA-ATAC overlapping cells: {rna5.n_obs}")
    rna5.write_h5ad(output_file, compression="gzip", compression_opts=6)

    return rna5


# -----------------------------
# Run external integration
# -----------------------------
def run_snrna_integration(input_file, prefix, batch_column="Sample", output_dir="."):
    if not SNRNA_INTEGRATE_SCRIPT.exists():
        raise FileNotFoundError(f"snRNA integration script not found: {SNRNA_INTEGRATE_SCRIPT}")

    cmd = [
        sys.executable,
        str(SNRNA_INTEGRATE_SCRIPT),
        "--input", input_file,
        "--output", output_dir,
        "--prefix", prefix,
        "--batch", batch_column,
    ]

    print("Running snRNA integration:")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


# -----------------------------
# Main
# -----------------------------
def main():
    print("Step 1: Read and merge raw CellBender-filtered h5ad files")
    adata = read_and_merge_h5ad_files(WORKDIR)

    print("Step 2: Annotate sample-level metadata")
    adata = annotate_samples(adata)

    print("Step 3: QC and doublet filtering")
    rna2 = run_qc_and_filter(adata)

    print("Step 4: Plot basic UMAP metadata")
    pl_umap(rna2, ["TimePoint", "CD45", "Replicate"], "rna2")

    print("Step 5: Remove doublet-enriched clusters")
    rna3 = subset_clean_clusters(rna2)

    print("Step 6: Annotate major cell types")
    rna3 = annotate_major_celltypes(rna3)
    # After the annotation of major cell types, we seperated the h5ad for each major cell type, which will be used for the downstream sub-clustering and annotation. 
    # We did a second round of doublet removal and quality control by removing the clusters that are enriched with doublets whithin each major cell type.
    # Cell subtype annotation was then performed for each major cell type separately based on the sub-clustering and marker gene analysis.
    # The overall h5ad file was then filtered to keep only the cells in the clean clusters, which will be used for the downstream anlysis.

    print("Step 7: Save obs table")
    rna3.obs.to_csv("rna3.Cellbender.noCorr.obs.csv")

    print("Step 8: Marker gene analysis")
    rna3 = run_marker_analysis(rna3, groupby="CellType")

    print("Step 9: Generate RNA-ATAC overlap object")
    rna5 = generate_rna_atac_overlap(
        rna_file="rna3.h5ad",
        atac_file=ATAC_FILE,
        output_file="rna5.Cellbender.overlap.h5ad",
    )

    print("Step 10: Run snRNA integration")
    run_snrna_integration(
        input_file="rna5.Cellbender.overlap.h5ad",
        prefix="rna5.Cellbender.overlap",
        batch_column="Sample",
        output_dir=OUTPUT_DIR,
    )

    print("Pipeline completed successfully.")


if __name__ == "__main__":
    main()