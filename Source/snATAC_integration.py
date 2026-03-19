#!/usr/bin/env python3

"""
Workflow for snATAC-seq data analysis:
generate integrated h5ad file based on muon / scanpy.
"""

from optparse import OptionParser
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
from muon import atac as ac


def integrate(
    adata,
    output=None,
    prefix="",
    batch=None,
    min_n_cells_by_counts=20,
    hvg=0,
    use_harmony=False,
    max_iter_harmony=50,
    use_bbknn=False,
    resolution=1.0,
    min_dist=0.5,
    leiden="leiden_clust",
    plot=None,
    compression=6,
    **kwargs,
):
    """
    Integrate ATAC h5ad data, run TF-IDF, LSI, optional Harmony/BBKNN,
    build neighbors/UMAP, Leiden clustering, and save results.
    """

    sc.settings.set_figure_params(dpi=350)

    print(f"Reading h5ad: {adata}")
    atac = sc.read_h5ad(adata)

    # ensure output path exists
    if output is None:
        output = "./"
    os.makedirs(output, exist_ok=True)

    # validate batch column
    if batch is not None and batch not in atac.obs.columns:
        print(f"Warning: batch column '{batch}' not found in atac.obs. batch=None will be used.")
        batch = None

    # remove batches with too few cells
    if batch is not None:
        print("Removing samples with fewer than 5 cells in the given batch column")
        count_of_cells = atac.obs.groupby(batch).size()
        valid_batches = count_of_cells[count_of_cells >= 5].index.tolist()
        atac = atac[atac.obs[batch].isin(valid_batches)].copy()

    # preserve raw counts
    if "raw" not in atac.layers:
        print('Copying .X to atac.layers["raw"]')
        atac.layers["raw"] = atac.X.copy()
    else:
        print('Copying atac.layers["raw"] to .X')
        atac.X = atac.layers["raw"].copy()

    print("Re-calculating QC metrics (skipped in current script)")
    # sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=True, inplace=True, layer="raw")

    # filter peaks by detection
    if min_n_cells_by_counts > 0 and "n_cells_by_counts" in atac.var.columns:
        print("Filtering peaks with low detection")
        mu.pp.filter_var(atac, "n_cells_by_counts", lambda x: x >= min_n_cells_by_counts)

    # TF-IDF normalization
    print("Running TF-IDF normalization")
    ac.pp.tfidf(atac, scale_factor=1e4)

    # LSI
    print("Running LSI")
    ac.tl.lsi(atac)

    # restore raw matrix to X
    atac.X = atac.layers["raw"].copy()
    atac.X = atac.X.astype(np.float32)
    del atac.layers["raw"]

    # remove first LSI component
    print("Removing the first LSI component")
    atac.obsm["X_lsi"] = atac.obsm["X_lsi"][:, 1:]
    atac.varm["LSI"] = atac.varm["LSI"][:, 1:]
    atac.uns["lsi"]["stdev"] = atac.uns["lsi"]["stdev"][1:]

    use_rep = "X_lsi"
    n_pcs = len(atac.uns["lsi"]["stdev"])

    # optional Harmony
    if batch is not None and use_harmony:
        print("Running Harmony")
        use_rep_adj = f"{use_rep}_harmony"
        sc.external.pp.harmony_integrate(
            atac,
            batch,
            basis=use_rep,
            adjusted_basis=use_rep_adj,
            max_iter_harmony=max_iter_harmony,
        )
        use_rep = use_rep_adj

    # optional BBKNN
    if batch is not None and use_bbknn:
        print("Running BBKNN")
        sc.external.pp.bbknn(atac, batch_key=batch, use_rep=use_rep, n_pcs=n_pcs)
    else:
        print(f"Running neighbors with use_rep={use_rep}")
        sc.pp.neighbors(atac, use_rep=use_rep, n_pcs=n_pcs)

    # UMAP
    print("Computing UMAP")
    sc.tl.umap(atac, min_dist=min_dist)

    # Leiden clustering
    print("Running Leiden clustering")
    sc.tl.leiden(atac, key_added=leiden, resolution=resolution)

    # optionally prefix leiden labels
    if prefix:
        atac.obs[leiden] = [f"{prefix}{v}" for v in atac.obs[leiden].astype(str)]

    # default plotting columns
    default_plots = [
        "Chr", "TSSEnrichment", "ReadsInTSS", "ReadsInPromoter",
        "ReadsInBlacklist", "PromoterRatio", "PassQC", "NucleosomeRatio",
        "nMultiFrags", "nMonoFrags", "nFrags", "nDiFrags", "DoubletScore",
        "DoubletEnrichment", "BlacklistRatio", "n_genes_by_counts",
        "log1p_n_genes_by_counts", "total_counts", "log1p_total_counts",
        "total_counts_Distal", "log1p_total_counts_Distal", "pct_counts_Distal",
        "total_counts_Promoter", "log1p_total_counts_Promoter", "pct_counts_Promoter",
        "total_counts_Exonic", "log1p_total_counts_Exonic", "pct_counts_Exonic",
        "region", "frag"
    ]

    if plot is not None:
        plot = np.union1d(default_plots, plot)
    else:
        plot = default_plots

    print("Generating UMAP plots")
    for col in plot:
        if col in atac.obs.columns:
            print(f"Plotting {col}")
            sc.pl.umap(atac, color=col, show=False, save=f".{prefix}.{col}.png")
            sc.pl.umap(
                atac,
                color=col,
                legend_loc="on data",
                legend_fontsize=5,
                show=False,
                save=f".{prefix}.{col}.ondata.png",
            )

    print(f"Plotting {leiden}")
    sc.pl.umap(atac, color=leiden, show=False, save=f"_{leiden}_beside.png")
    sc.pl.umap(
        atac,
        color=leiden,
        legend_loc="on data",
        legend_fontsize=4,
        show=False,
        save=f"_{leiden}_ondata.png",
    )

    # save result
    output_file = os.path.join(output, f"{prefix}.h5ad" if prefix else "integrated.h5ad")
    print(f"Writing h5ad to: {output_file}")
    atac.write_h5ad(output_file, compression="gzip")

    print("Done.")


if __name__ == "__main__":
    usage = "\nWorkflow for snATAC-seq data analysis - generating h5ad file (based on muon, scanpy)"
    parser = OptionParser(usage, version="Version 1.0")

    parser.add_option("-i", "--input", dest="h5ad", type="string", help="Input h5ad file")
    parser.add_option("-o", "--output", dest="output", type="string", help="Output directory")
    parser.add_option(
        "-p", "--prefix",
        dest="prefix",
        type="string",
        default="",
        help='Prefix added to Leiden cluster labels'
    )
    parser.add_option("-b", "--batch", default=None, type="str", help="Batch column name")
    parser.add_option(
        "--min_n_cells_by_counts",
        dest="min_n_cells_by_counts",
        type="int",
        default=20,
        help="Minimum number of cells by counts for filtering peaks"
    )
    parser.add_option("--hvg", default=0, type="int")
    parser.add_option("--use_harmony", dest="use_harmony", default=False, action="store_true")
    parser.add_option("--max_iter_harmony", default=50, type="int")
    parser.add_option("--use_bbknn", dest="use_bbknn", default=False, action="store_true")
    parser.add_option("-r", "--resolution", default=1.0, type="float", help="Resolution for Leiden")
    parser.add_option("--min_dist", default=0.5, type="float", help="Minimum distance for UMAP")
    parser.add_option("-l", "--leiden", default="leiden_clust", type="str", help="Leiden clustering column")
    parser.add_option("--plot", default=None, type="str", action="append")
    parser.add_option("--compression", default=6, type="int")

    options, args = parser.parse_args()

    if not options.h5ad:
        parser.error("Please provide --input")
    if not options.output:
        parser.error("Please provide --output")

    integrate(
        adata=options.h5ad,
        output=options.output,
        prefix=options.prefix,
        batch=options.batch,
        min_n_cells_by_counts=options.min_n_cells_by_counts,
        hvg=options.hvg,
        use_harmony=options.use_harmony,
        max_iter_harmony=options.max_iter_harmony,
        use_bbknn=options.use_bbknn,
        resolution=options.resolution,
        min_dist=options.min_dist,
        leiden=options.leiden,
        plot=options.plot,
        compression=options.compression,
    )