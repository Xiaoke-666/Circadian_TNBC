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
OUTPUT_DIR = PROJECT_DIR / "results" / "circadian_tf"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TF_H5AD = DATA_DIR / "ChromVar_Deviations.h5ad"
RNA_METADATA_H5AD = DATA_DIR / "rna_metadata_reference.h5ad"

N_REPEATS = 10
N_PSEUDO_REPS = 3


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
def pseudobulk_normalized(adata: sc.AnnData, group_col: str, method: str = "mean") -> sc.AnnData:
    """
    Generate pseudobulk profiles from TF deviation scores.

    Parameters
    ----------
    adata : AnnData
        Input AnnData containing TF deviation scores.
    group_col : str
        Column in `adata.obs` used to define pseudobulk groups.
    method : str
        Aggregation method: 'mean' or 'median'.

    Returns
    -------
    AnnData
        Pseudobulk AnnData object.
    """
    assert method in ["mean", "median"]

    if "raw" not in adata.layers:
        print("Storing current .X as .layers['raw']")
        adata.layers["raw"] = adata.X.copy()

    groups = adata.obs[group_col].astype(str).unique()
    pseudobulks = []

    for group in groups:
        print(f"Processing group: {group}")
        subset = adata[adata.obs[group_col] == group].copy()

        if method == "mean":
            values = subset.X.mean(axis=0)
        else:
            if hasattr(subset.X, "toarray"):
                values = np.median(subset.X.toarray(), axis=0)
            else:
                values = np.median(subset.X, axis=0)

        pb = sc.AnnData(X=np.asarray(values).reshape(1, -1))
        pb.obs["group"] = group
        pb.var = subset.var.copy()
        pseudobulks.append(pb)

    pb_adata = sc.concat(pseudobulks)
    pb_adata.obs.index = pb_adata.obs["group"]
    pb_adata.obs.index.name = None

    return pb_adata


def assign_pseudo_replicates(obs_df: pd.DataFrame, group_col: str, n_reps: int, seed: int) -> pd.DataFrame:
    """
    Randomly assign cells within each group to pseudo-replicates.
    """
    rng = np.random.default_rng(seed)
    obs_df = obs_df.copy()

    assigned = []
    for group, group_df in obs_df.groupby(group_col, sort=False):
        group_df = group_df.copy()
        group_df["pseudo_replicate"] = rng.integers(0, n_reps, size=group_df.shape[0])
        group_df[f"{group_col}.pseudo_replicate"] = (
            group_df[group_col].astype(str) + "." + group_df["pseudo_replicate"].astype(str)
        )
        assigned.append(group_df)

    return pd.concat(assigned, axis=0)


def build_metadata(
    adata: sc.AnnData,
    celltype_col: str,
    timepoint_col: str,
    pseudorep_col: str,
    sample_col: str,
) -> pd.DataFrame:
    """
    Build metadata table for pseudobulk outputs.
    """
    metadata = adata.obs[
        [celltype_col, timepoint_col, pseudorep_col, sample_col]
    ].drop_duplicates().copy()

    metadata.columns = ["Celltype", "TimePoint", "Replicate", "SampleName"]
    metadata["Replicate"] = (
        pd.to_numeric(metadata["Replicate"], errors="coerce")
        .fillna(0)
        .astype(int) + 1
    )

    return metadata


def export_pseudobulk(
    pb: sc.AnnData,
    metadata: pd.DataFrame,
    matrix_file: Path,
    metadata_file: Path,
) -> None:
    """
    Export pseudobulk matrix and metadata.
    """
    pb_df = pd.DataFrame(pb.X.T, index=pb.var.index, columns=pb.obs.index)
    pb_df.columns = metadata["SampleName"]
    pb_df.index.name = "GeneID"

    if str(matrix_file).endswith(".gz"):
        pb_df.to_csv(matrix_file, sep="\t", compression="gzip")
    else:
        pb_df.to_csv(matrix_file, sep="\t")

    metadata.to_csv(metadata_file, sep="\t", index=False)


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------
def main():
    print("Loading TF deviation matrix...")
    adata = sc.read(TF_H5AD)

    print("Loading RNA reference metadata...")
    rna = sc.read(RNA_METADATA_H5AD)

    # Keep only cells with valid CellType_TimePoint labels
    adata = adata[adata.obs["CellType_TimePoint"].notna()].copy()

    # Transfer metadata from RNA object
    for col in rna.obs.columns:
        if col not in adata.obs.columns:
            adata.obs[col] = adata.obs.index.map(rna.obs[col])

    # Standardize TF names
    adata.var["TF"] = adata.var["name"].astype(str).str.split("_").str[0]
    adata.var.index = adata.var["TF"]
    adata.var.index.name = None

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    # Save cleaned input object
    adata.write(OUTPUT_DIR / "ChromVar_Deviations.cleaned.h5ad", compression="gzip")

    # -------------------------------------------------------------------------
    # Repeated pseudobulk generation
    # -------------------------------------------------------------------------
    print("Generating repeated pseudobulk TF matrices...")

    for i in range(N_REPEATS):
        print(f"\n========== Repeat {i + 1}/{N_REPEATS} ==========")

        assigned_obs = assign_pseudo_replicates(
            adata.obs,
            group_col="CellType_TimePoint",
            n_reps=N_PSEUDO_REPS,
            seed=i,
        )

        adata.obs["pseudo_replicate"] = assigned_obs.loc[adata.obs.index, "pseudo_replicate"]
        adata.obs["CellType_TimePoint.pseudo_replicate"] = assigned_obs.loc[
            adata.obs.index, "CellType_TimePoint.pseudo_replicate"
        ]

        print(adata.obs["CellType_TimePoint.pseudo_replicate"].value_counts())

        pb = pseudobulk_normalized(
            adata,
            group_col="CellType_TimePoint.pseudo_replicate",
            method="mean",
        )

        metadata = build_metadata(
            adata=adata,
            celltype_col="ATAC_CellType",
            timepoint_col="TimePoint",
            pseudorep_col="pseudo_replicate",
            sample_col="CellType_TimePoint.pseudo_replicate",
        )

        export_pseudobulk(
            pb=pb,
            metadata=metadata,
            matrix_file=OUTPUT_DIR / f"CellType_TimePoint.pseudo_replicate_{i}.tsv.gz",
            metadata_file=OUTPUT_DIR / f"CellType_TimePoint.pseudo_replicate_metadata_{i}.tsv",
        )

    # -------------------------------------------------------------------------
    # One overall pseudobulk export
    # -------------------------------------------------------------------------
    print("\nGenerating one overall TF pseudobulk matrix...")

    assigned_obs = assign_pseudo_replicates(
        adata.obs,
        group_col="CellType_TimePoint",
        n_reps=N_PSEUDO_REPS,
        seed=0,
    )

    adata.obs["pseudo_replicate"] = assigned_obs.loc[adata.obs.index, "pseudo_replicate"]
    adata.obs["CellType_TimePoint.pseudo_replicate"] = assigned_obs.loc[
        adata.obs.index, "CellType_TimePoint.pseudo_replicate"
    ]

    pb = pseudobulk_normalized(
        adata,
        group_col="CellType_TimePoint.pseudo_replicate",
        method="mean",
    )

    metadata = build_metadata(
        adata=adata,
        celltype_col="ATAC_CellType",
        timepoint_col="TimePoint",
        pseudorep_col="pseudo_replicate",
        sample_col="CellType_TimePoint.pseudo_replicate",
    )

    export_pseudobulk(
        pb=pb,
        metadata=metadata,
        matrix_file=OUTPUT_DIR / "CellType_TimePoint.pseudo_replicate.tsv",
        metadata_file=OUTPUT_DIR / "CellType_TimePoint.pseudo_replicate_metadata.tsv",
    )

    print("\nPseudo-replicate TF matrices generated successfully.")


if __name__ == "__main__":
    main()