#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np
import glob
import argparse
import os


def compute_tf_connectivity(
    h5ad_file,
    tsv_pattern,
    output_file,
    min_mean=0.02,
    max_mean=4,
    min_disp=0.5
):
    """
    Compute TF connectivity network from chromVAR deviations matrix.

    Parameters
    ----------
    h5ad_file : str
        Path to input .h5ad file
    tsv_pattern : str
        Pattern for circadian gene .tsv files (e.g., "*.tsv")
    output_file : str
        Output CSV file for edges
    """

    # Read TF matrix
    tf_mat = sc.read_h5ad(h5ad_file)

    # Read circadian gene files
    files = glob.glob(tsv_pattern)
    if len(files) == 0:
        raise ValueError(f"No TSV files found with pattern: {tsv_pattern}")

    df_list = [pd.read_table(f, sep="\t") for f in files]
    df_comb = pd.concat(df_list, ignore_index=True)

    if "GeneID" not in df_comb.columns:
        raise ValueError("TSV files must contain 'GeneID' column")

    genes = df_comb["GeneID"].unique()

    # Subset TF matrix
    sub_tf = tf_mat[tf_mat.obs_names.isin(genes)].copy()
    sub_tf = sub_tf[~sub_tf.obs_names.duplicated()].copy()

    if sub_tf.shape[0] == 0:
        raise ValueError("No matching genes found in TF matrix.")

    # Compute neighbors
    sc.pp.highly_variable_genes(
        sub_tf,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp
    )
    sc.pp.pca(sub_tf)
    sc.pp.neighbors(sub_tf, use_rep="X_pca")

    # Extract connectivity and distance matrices
    conn = sub_tf.obsp["connectivities"]
    dis = sub_tf.obsp["distances"]

    # Convert to COO format
    coo_conn = conn.tocoo()
    coo_dist = dis.tocoo()

    # Build edges dataframe
    edges = pd.DataFrame({
        "source": [sub_tf.obs_names[i] for i in coo_conn.row],
        "target": [sub_tf.obs_names[j] for j in coo_conn.col],
        "connectivity": coo_conn.data
    })

    # Distance lookup
    dist_dict = {
        (i, j): d for i, j, d in zip(coo_dist.row, coo_dist.col, coo_dist.data)
    }

    edges["distance"] = [
        dist_dict.get((i, j), np.nan)
        for i, j in zip(coo_conn.row, coo_conn.col)
    ]

    # Remove self-loops and duplicates
    edges = edges[edges["source"] != edges["target"]]
    edges = edges.drop_duplicates(subset=["source", "target"])

    # Save output
    edges.to_csv(output_file, index=False)

    print(f"Saved connectivity network to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute TF connectivity network from chromVAR deviations"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input h5ad file"
    )

    parser.add_argument(
        "-p", "--pattern",
        default="*.tsv",
        help="Pattern for circadian gene TSV files (default: *.tsv)"
    )

    parser.add_argument(
        "-o", "--output",
        default="TF_connectivity.csv",
        help="Output CSV file (default: TF_connectivity.csv)"
    )

    parser.add_argument("--min_mean", type=float, default=0.02)
    parser.add_argument("--max_mean", type=float, default=4)
    parser.add_argument("--min_disp", type=float, default=0.5)

    args = parser.parse_args()

    compute_tf_connectivity(
        h5ad_file=args.input,
        tsv_pattern=args.pattern,
        output_file=args.output,
        min_mean=args.min_mean,
        max_mean=args.max_mean,
        min_disp=args.min_disp
    )


if __name__ == "__main__":
    main()

