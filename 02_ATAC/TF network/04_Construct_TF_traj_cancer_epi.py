#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import os

from matplotlib import rc_context
import matplotlib.pyplot as plt
from adjustText import adjust_text
from matplotlib.collections import LineCollection


def run_cancer_tf_paga(
    h5ad_file,
    cancer_tsv,
    tf_rank_file,
    output_prefix,
    min_mean=0.02,
    max_mean=4,
    min_disp=0.5
):

    # ---------------------------
    # Load data
    # ---------------------------
    tf_mat = sc.read_h5ad(h5ad_file)

    df_cancer = pd.read_table(cancer_tsv)
    tf_filtered = pd.read_csv(tf_rank_file)

    # Filter Cancer_Epi TFs
    tf_filtered = tf_filtered.loc[tf_filtered["celltype"] == "Cancer_Epi", :]
    df_cancer = df_cancer.loc[df_cancer["GeneID"].isin(tf_filtered["motif_name"]), :]

    TFs = df_cancer["GeneID"].unique()

    if len(TFs) == 0:
        raise ValueError("No TFs found after filtering.")

    # ---------------------------
    # Subset AnnData
    # ---------------------------
    tf_mat = tf_mat[~tf_mat.obs_names.duplicated()].copy()

    if "TF" not in tf_mat.obs.columns:
        raise ValueError("Column 'TF' not found in obs.")

    if "ATAC.Subclass.June18_2025" not in tf_mat.var.columns:
        raise ValueError("Column 'ATAC.Subclass.June18_2025' not found in var.")

    sub = tf_mat[tf_mat.obs["TF"].isin(TFs), :].copy()
    sub = sub[:, sub.var["ATAC.Subclass.June18_2025"].isin(["Cancer_Epi"])].copy()

    if sub.shape[0] == 0 or sub.shape[1] == 0:
        raise ValueError("Subset resulted in empty matrix.")

    # ---------------------------
    # Compute neighbors
    # ---------------------------
    sc.pp.highly_variable_genes(sub, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
    sc.pp.pca(sub)
    sc.pp.neighbors(sub, use_rep="X_pca")

    # ---------------------------
    # Build edge list
    # ---------------------------
    conn = sub.obsp["connectivities"]
    dis = sub.obsp["distances"]

    coo_conn = conn.tocoo()
    coo_dist = dis.tocoo()

    edges = pd.DataFrame({
        "source": [sub.obs_names[i] for i in coo_conn.row],
        "target": [sub.obs_names[j] for j in coo_conn.col],
        "connectivity": coo_conn.data
    })

    dist_dict = {(i, j): d for i, j, d in zip(coo_dist.row, coo_dist.col, coo_dist.data)}

    edges["distance"] = [
        dist_dict.get((i, j), np.nan)
        for i, j in zip(coo_conn.row, coo_conn.col)
    ]

    edges = edges[edges["source"] != edges["target"]]
    edges = edges.drop_duplicates(subset=["source", "target"])

    # ---------------------------
    # Annotate edges
    # ---------------------------
    edges_anno = edges.merge(df_cancer, left_on="source", right_on="GeneID", how="left")

    edges_anno = edges_anno[
        ["source", "target", "connectivity", "distance", "acrophase", "pvalue"]
    ]

    # Discretize phase
    bins = [0, 3, 6, 9, 12, 15, 18, 21, 24]
    labels = ["P1","P2","P3","P4","P5","P6","P7","P8"]

    edges_anno["phase"] = pd.cut(
        edges_anno["acrophase"],
        bins=bins,
        labels=labels,
        right=False
    )

    # ---------------------------
    # Map phase back to AnnData
    # ---------------------------
    tf_phase = edges_anno[["source","phase"]].drop_duplicates()
    phase_map = dict(zip(tf_phase["source"], tf_phase["phase"]))

    sub.obs["phase"] = sub.obs["TF"].map(phase_map)

    # ---------------------------
    # PAGA + UMAP
    # ---------------------------
    sc.tl.paga(sub, groups="phase")
    sc.pl.paga(sub, plot=False)

    sc.tl.umap(sub, init_pos="paga", spread=1.2)

    # Save AnnData
    out_h5ad = f"{output_prefix}.h5ad"
    sub.write_h5ad(out_h5ad, compression="gzip")

    # ---------------------------
    # Plot
    # ---------------------------
    sub.uns["phase_colors"] = [
        "#FF4545","#FF9C73","#FBD288","#FCF596",
        "#DFF2EB","#B9E5E8","#7AB2D3","#4A628A"
    ]

    TFs_label = [
        "Jun","Junb","Jund","Fos","Fosb","Fosl1",
        "Tead3","Tead4","Smad2","Smad3",
        "Ppara","Ppard","Pparg","Thra","Thrb","Rarg",
        "Nfkb1","Rel","Rela","Relb","Bcl6",
        "Bach1","Bach2","Nfe2l2","Nr1d1","Nr1d2",
        "Rora","Rorb"
    ]

    with rc_context({"figure.figsize": (5, 4)}):
        sc.pl.paga_compare(
            sub,
            edges=True,
            color="phase",
            frameon=False,
            legend_loc="none",
            edge_width_scale=0.3,
            show=False,
            node_size_scale=3.0,
            size=150.0
        )

    fig = plt.gcf()
    axes_list = fig.axes
    ax_umap = axes_list[-2]
    ax_paga = axes_list[-1]

    # Improve aesthetics
    for coll in ax_umap.collections:
        coll.set_linewidth(0.1)
        coll.set_edgecolor("#C0C0C0")

    for artist in ax_paga.collections:
        if isinstance(artist, LineCollection):
            artist.set_color("gray")
            artist.set_alpha(0.8)

    # Label TFs
    coords = sub.obsm["X_umap"]
    obs_names = np.array(sub.obs_names)

    texts = []
    for tf in TFs_label:
        if tf in obs_names:
            i = np.where(obs_names == tf)[0][0]
            x, y = coords[i, 0], coords[i, 1]
            txt = ax_umap.text(x, y, tf, fontsize=10, weight="bold")
            texts.append(txt)

    adjust_text(
        texts,
        ax=ax_umap,
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
        force_points=0.5,
        force_text=0.5
    )

    out_pdf = f"{output_prefix}.pdf"
    plt.savefig(out_pdf, bbox_inches="tight")

    print(f"Saved: {out_h5ad}")
    print(f"Saved: {out_pdf}")


# ---------------------------
# CLI
# ---------------------------
def main():
    parser = argparse.ArgumentParser(description="Cancer_Epi TF circadian PAGA analysis")

    parser.add_argument("-i", "--input", required=True, help="Input h5ad file")
    parser.add_argument("-c", "--cancer", required=True, help="Cancer circadian TF TSV")
    parser.add_argument("-t", "--tf_rank", required=True, help="Filtered cancer TF CSV")
    parser.add_argument("-o", "--output", default="Cancer_circadian")

    parser.add_argument("--min_mean", type=float, default=0.02)
    parser.add_argument("--max_mean", type=float, default=4)
    parser.add_argument("--min_disp", type=float, default=0.5)

    args = parser.parse_args()

    run_cancer_tf_paga(
        h5ad_file=args.input,
        cancer_tsv=args.cancer,
        tf_rank_file=args.tf_rank,
        output_prefix=args.output,
        min_mean=args.min_mean,
        max_mean=args.max_mean,
        min_disp=args.min_disp
    )


if __name__ == "__main__":
    main()