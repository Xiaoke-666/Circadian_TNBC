#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"

TF_DIR = RESULTS_DIR / "circadian_tf"
COACTIVITY_DIR = TF_DIR / "TF_motif_coactivity"
COACTIVITY_DIR.mkdir(parents=True, exist_ok=True)

CHROMVAR_FILE = DATA_DIR / "ChromVar_Deviations.h5ad"
RNA_METADATA_FILE = DATA_DIR / "rna_metadata.tsv.gz"

# Example cell type
CELLTYPE = "CD4_T"
CELLTYPE_MODEL_FILE = TF_DIR / f"ChromVar_Deviations.{CELLTYPE}.Circadian.h5ad"
CIRCADIAN_TF_FILE = TF_DIR / f"{CELLTYPE}.circadianTF.significant.tsv"


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
def assign_phase_group(acrophase_series: pd.Series) -> pd.Series:
    """
    Assign circadian phase groups based on acrophase.

    Phase bins:
    P1: [0,3), P2: [3,6), ..., P8: [21,24]
    """
    return pd.cut(
        acrophase_series,
        bins=[-1, 3, 6, 9, 12, 15, 18, 21, 24],
        labels=["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"]
    )


def prepare_clean_chromvar_h5ad(
    chromvar_file: Path,
    rna_metadata_file: Path,
    output_file: Path,
) -> sc.AnnData:
    """
    Load chromVAR deviation matrix, transfer selected metadata from RNA annotations,
    and save a cleaned AnnData object.
    """
    print("Loading chromVAR deviation matrix...")
    adata = sc.read(chromvar_file)

    print("Loading RNA metadata...")
    rna_meta = pd.read_csv(rna_metadata_file, sep="\t", index_col=0)

    for col in ["ATAC_CellType", "ATAC_Subtype"]:
        if col in rna_meta.columns:
            adata.obs[col] = adata.obs.index.map(rna_meta[col])

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    print(f"Writing cleaned chromVAR object: {output_file}")
    adata.write(output_file, compression="gzip")

    return adata


def transpose_and_save(adata: sc.AnnData, output_file: Path) -> sc.AnnData:
    """
    Transpose AnnData so TFs become observations for downstream TF-space modeling.
    """
    adata_t = adata.T
    adata_t.obs_names_make_unique()
    adata_t.var_names_make_unique()

    print(f"Writing transposed chromVAR object: {output_file}")
    adata_t.write(output_file, compression="gzip")
    return adata_t


def run_tf_manifold(
    adata: sc.AnnData,
    min_dist: float = 0.6,
    leiden_key: str = "leiden",
    resolution: float = 3.0,
) -> sc.AnnData:
    """
    Run PCA, neighbors, UMAP, tSNE, and Leiden clustering in TF space.

    Notes
    -----
    The input AnnData is expected to be transposed such that each observation
    represents a TF rather than a cell.
    """
    adata = adata.copy()

    print("Running PCA...")
    sc.pp.pca(adata)

    print("Running neighbors...")
    sc.pp.neighbors(adata, use_rep="X_pca")

    print("Running UMAP...")
    sc.tl.umap(adata, min_dist=min_dist)

    print("Running Leiden clustering...")
    sc.tl.leiden(adata, resolution=resolution, key_added=leiden_key)

    return adata


def build_connectivity_distance_table(adata: sc.AnnData) -> pd.DataFrame:
    """
    Extract connectivities and distances from AnnData graph and return long-format table.
    """
    print("Extracting connectivities...")
    connectivities = pd.DataFrame(
        adata.obsp["connectivities"].todense(),
        index=adata.obs.index,
        columns=adata.obs.index,
    )
    connectivities_long = connectivities.stack().reset_index()
    connectivities_long.columns = ["cell1", "cell2", "connectivity"]

    print("Extracting distances...")
    distances = pd.DataFrame(
        adata.obsp["distances"].todense(),
        index=adata.obs.index,
        columns=adata.obs.index,
    )
    distances_long = distances.stack().reset_index()
    distances_long.columns = ["cell1", "cell2", "distance"]

    merged = connectivities_long.merge(
        distances_long,
        on=["cell1", "cell2"],
        how="inner",
    )
    return merged


def annotate_with_circadian_tf_stats(
    edge_df: pd.DataFrame,
    circadian_tf_df: pd.DataFrame,
    key_col: str = "cell1",
) -> pd.DataFrame:
    """
    Map circadian TF rhythmicity statistics onto the edge table.
    Assumes node IDs in `key_col` match `GeneID` in circadian_tf_df.
    """
    circadian_tf_df = circadian_tf_df.copy()
    circadian_tf_df["Phase_group"] = assign_phase_group(circadian_tf_df["acrophase"])

    stat_cols = [
        "acrophase", "amplitude", "pvalue", "qvalue",
        "Rsq", "mesor", "sincoef", "coscoef", "Phase_group"
    ]

    circadian_tf_df = circadian_tf_df.set_index("GeneID")

    for col in stat_cols:
        edge_df[col] = edge_df[key_col].map(circadian_tf_df[col])

    return edge_df


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------
def main():
    # -------------------------------------------------------------------------
    # Step 1. Prepare cleaned chromVAR object
    # -------------------------------------------------------------------------
    clean_h5ad = COACTIVITY_DIR / "ChromVar_Deviations.Clean.h5ad"
    clean_t_h5ad = COACTIVITY_DIR / "ChromVar_Deviations.Clean.T.h5ad"
    tf_manifold_file = COACTIVITY_DIR / "ChromVar_Deviations.Clean.T.manifold.h5ad"

    if not clean_h5ad.exists():
        adata = prepare_clean_chromvar_h5ad(
            chromvar_file=CHROMVAR_FILE,
            rna_metadata_file=RNA_METADATA_FILE,
            output_file=clean_h5ad,
        )
    else:
        adata = sc.read(clean_h5ad)

    # -------------------------------------------------------------------------
    # Step 2. Transpose to TF space
    # -------------------------------------------------------------------------
    if not clean_t_h5ad.exists():
        adata_t = transpose_and_save(adata, clean_t_h5ad)
    else:
        adata_t = sc.read(clean_t_h5ad)

    # -------------------------------------------------------------------------
    # Step 3. Build TF manifold model (PCA / UMAP / Leiden)
    # -------------------------------------------------------------------------
    if not tf_manifold_file.exists():
        print("Building TF manifold model...")
        adata_t_model = run_tf_manifold(
            adata=adata_t,
            hvg=0,
            min_dist=0.6,
            leiden_key="leiden",
            resolution=3.0,
        )
        print(f"Writing TF manifold model: {tf_manifold_file}")
        adata_t_model.write(tf_manifold_file, compression="gzip")
    else:
        print("Loading existing TF manifold model...")
        adata_t_model = sc.read(tf_manifold_file)

    # -------------------------------------------------------------------------
    # Step 4. Load integrated circadian TF model for one example cell type
    # -------------------------------------------------------------------------
    print(f"Loading circadian TF model for {CELLTYPE}...")
    tf_model = sc.read(CELLTYPE_MODEL_FILE)

    # -------------------------------------------------------------------------
    # Step 5. Export connectivity and distance matrix
    # -------------------------------------------------------------------------
    graph_df = build_connectivity_distance_table(tf_model)

    out_full = COACTIVITY_DIR / f"{CELLTYPE}.connectivity_distance.csv"
    graph_df.to_csv(out_full, index=False)

    graph_df = graph_df[graph_df["connectivity"] > 0].copy()

    # -------------------------------------------------------------------------
    # Step 6. Load significant circadian TFs
    # -------------------------------------------------------------------------
    print("Loading circadian TF statistics...")
    circadian_tfs = pd.read_csv(CIRCADIAN_TF_FILE, sep="\t")

    # -------------------------------------------------------------------------
    # Step 7. Map rhythmicity stats onto graph edges
    # -------------------------------------------------------------------------
    graph_df = annotate_with_circadian_tf_stats(
        edge_df=graph_df,
        circadian_tf_df=circadian_tfs,
        key_col="cell1",
    )

    out_filtered = COACTIVITY_DIR / f"{CELLTYPE}.connectivity_distance.filtered.csv"
    graph_df.to_csv(out_filtered, index=False)

    print("Done.")


if __name__ == "__main__":
    main()