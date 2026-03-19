#!/usr/bin/env python3

"""
Circadian gene module analysis

This script:
1. Subsets RNA data to circadian genes
2. Performs clustering (Leiden) to define modules
3. Annotates circadian genes with modules
4. Computes module-level correlations
5. Visualizes UMAP and acrophase patterns
"""

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Paths (relative, GitHub-safe)
# -----------------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_DIR / "data"
RESULT_DIR = PROJECT_DIR / "results"

RNA_FILE = DATA_DIR / "rna.h5ad"
CIRCADIAN_GENE_FILE = DATA_DIR / "circadian_genes.tsv"

OUTPUT_DIR = RESULT_DIR / "circadian_modules"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def plot_umap_clean(adata, keys, prefix):
    """
    Plot clean UMAP (no axes/legend) for figure usage.
    """
    for key in keys:
        print(f"Plotting {key}")
        sc.pl.umap(
            adata,
            color=key,
            frameon=False,
            cmap="Spectral_r",
            vmin="p20",
            vmax="p95",
            show=False,
        )
        plt.savefig(
            OUTPUT_DIR / f"UMAP.{prefix}.{key}.png",
            dpi=300,
            bbox_inches="tight",
            transparent=True,
        )
        plt.close()


# -----------------------------------------------------------------------------
# Load and preprocess
# -----------------------------------------------------------------------------
def load_data():
    print("Loading RNA data...")
    rna = sc.read(RNA_FILE)

    # use raw counts
    if "raw" in rna.layers:
        rna.X = rna.layers["raw"].copy()
        del rna.layers["raw"]

    # transpose to gene × cell
    rna = rna.T

    circadian_genes = pd.read_csv(CIRCADIAN_GENE_FILE, sep="\t")

    return rna, circadian_genes


def subset_circadian_genes(rna, circadian_genes):
    print("Subsetting circadian genes...")
    genes = circadian_genes["GeneID"].unique().tolist()
    rna_sub = rna[rna.obs.index.isin(genes)].copy()
    return rna_sub


# -----------------------------------------------------------------------------
# Clustering (module definition)
# -----------------------------------------------------------------------------
def compute_modules(rna):
    print("Running UMAP + Leiden clustering...")
    sc.tl.umap(rna, min_dist=0.01)
    sc.tl.leiden(rna, resolution=1.2, key_added="leiden")

    # rename clusters to M1, M2, ...
    rna.obs["module"] = ["M" + str(int(x) + 1) for x in rna.obs["leiden"]]
    return rna


# -----------------------------------------------------------------------------
# Export gene modules
# -----------------------------------------------------------------------------
def export_module_gene_lists(rna):
    print("Exporting gene lists per module...")

    data = {}
    max_len = 0

    modules = sorted(rna.obs["module"].unique())

    for m in modules:
        genes = rna.obs[rna.obs["module"] == m].index.tolist()
        data[m] = genes
        max_len = max(max_len, len(genes))

    # pad
    for k in data:
        data[k] += [np.nan] * (max_len - len(data[k]))

    df = pd.DataFrame(data)
    df.to_csv(OUTPUT_DIR / "CircadianGeneModules.tsv", sep="\t", index=False)


# -----------------------------------------------------------------------------
# Correlation analysis
# -----------------------------------------------------------------------------
def compute_module_correlation(rna):
    print("Computing module correlation...")

    df = pd.DataFrame(rna.obsm["X_pca"])
    df["module"] = rna.obs["module"].values

    df_mean = df.groupby("module").mean()
    corr = df_mean.T.corr()

    corr.to_csv(OUTPUT_DIR / "Module_correlation.tsv", sep="\t")


# -----------------------------------------------------------------------------
# Map acrophase from DiscoRhythm
# -----------------------------------------------------------------------------
def map_acrophase(rna, circadian_dir):
    print("Mapping acrophase...")

    celltypes = [
        "CD4_T", "CD8_T", "B", "Cancer_Epi",
        "Macro_Mono", "NKT", "Neu", "CAFs"
    ]

    for ct in celltypes:
        file = circadian_dir / f"{ct}.discoODAres.tsv"
        if not file.exists():
            continue

        df = pd.read_csv(file, sep="\t")
        mapping = df.set_index("GeneID")["acrophase"].to_dict()

        rna.obs[f"acrophase.{ct}"] = rna.obs.index.map(mapping)

    return rna


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    rna, circadian_genes = load_data()

    rna = subset_circadian_genes(rna, circadian_genes)

    rna.write(OUTPUT_DIR / "rna_circadian.h5ad")

    rna = compute_modules(rna)

    export_module_gene_lists(rna)

    compute_module_correlation(rna)

    # optional: map acrophase
    circadian_dir = RESULT_DIR / "CircadianGenes"
    rna = map_acrophase(rna, circadian_dir)

    plot_umap_clean(
        rna,
        ["module"],
        prefix="CircadianModules"
    )

    print("Done.")


if __name__ == "__main__":
    main()