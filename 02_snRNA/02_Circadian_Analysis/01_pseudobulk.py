from pathlib import Path
import os

import numpy as np
import pandas as pd
import scanpy as sc


# -----------------------------------------------------------------------------
# Global settings
# -----------------------------------------------------------------------------
sc.settings.set_figure_params(
    scanpy=True,
    dpi=300,
    dpi_save=300,
    fontsize=6,
    color_map="viridis",
)
pd.options.display.max_rows = 999

RANDOM_SEED = 42


# -----------------------------------------------------------------------------
# Core functions
# -----------------------------------------------------------------------------
def pseudobulk_by_group(adata, group_col):
    """
    Generate pseudobulk expression profiles by summing counts within each group.

    Parameters
    ----------
    adata : AnnData
        Input single-cell AnnData object.
    group_col : str
        Column in `adata.obs` defining groups to aggregate.

    Returns
    -------
    AnnData
        Pseudobulk AnnData with one sample per group.
    """
    adata = adata.copy()

    if "raw" in adata.layers:
        print('Copying adata.layers["raw"] to adata.X')
        adata.X = adata.layers["raw"].copy()
    else:
        print('Copying adata.X to adata.layers["raw"]')
        adata.layers["raw"] = adata.X.copy()

    groups = adata.obs[group_col].astype(str).unique()
    pseudobulks = []

    for group in groups:
        print(f"Aggregating: {group}")
        subset = adata[adata.obs[group_col] == group].copy()

        pb = sc.AnnData(
            X=subset.X.sum(axis=0),
            var=subset.var.copy(),
        )
        pb.obs["group"] = group
        pseudobulks.append(pb)

    pb_adata = sc.concat(pseudobulks)
    pb_adata.obs.index = pb_adata.obs["group"]
    pb_adata.obs.index.name = None

    pb_adata.layers["raw"] = pb_adata.X.copy()
    sc.pp.normalize_total(pb_adata, target_sum=1e4)
    sc.pp.log1p(pb_adata)

    return pb_adata


def sanitize_strings(series):
    """
    Sanitize strings for output sample names and cell type labels.
    """
    return (
        series.astype(str)
        .str.replace("/", "_", regex=False)
        .str.replace(" ", "_", regex=False)
        .str.replace("+", "", regex=False)
    )


def assign_pseudo_replicates(obs_df, group_col, n_reps=3, seed=42):
    """
    Randomly assign cells within each group to pseudo-replicates.

    Parameters
    ----------
    obs_df : pandas.DataFrame
        Cell metadata table.
    group_col : str
        Column defining grouping (e.g. CellType_TimePoint).
    n_reps : int
        Number of pseudo-replicates per group.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    pandas.DataFrame
        Updated metadata with pseudo-replicate assignments.
    """
    rng = np.random.default_rng(seed)
    obs_df = obs_df.copy()

    pseudo_rep_list = []
    for group, group_df in obs_df.groupby(group_col, sort=False):
        print(f"Assigning pseudo-replicates: {group}")
        group_df = group_df.copy()
        group_df["pseudo_replicate"] = rng.integers(0, n_reps, size=group_df.shape[0])
        pseudo_rep_list.append(group_df)

    return pd.concat(pseudo_rep_list, axis=0)


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------
def main():
    workdir = Path(
        "./"
    )
    os.chdir(workdir)

    # Read input AnnData object
    rna_overlap = sc.read("RNA_overlap.h5ad")

    # Create combined label: CellType.TimePoint
    rna_overlap.obs["CellType_TimePoint"] = (
        rna_overlap.obs["CellType"].astype(str)
        + "."
        + rna_overlap.obs["TimePoint"].astype(str)
    )

    # Assign pseudo-replicates within each CellType-TimePoint group
    obs_with_reps = assign_pseudo_replicates(
        rna_overlap.obs,
        group_col="CellType_TimePoint",
        n_reps=3,
        seed=RANDOM_SEED,
    )

    rna_overlap.obs["pseudo_replicate"] = obs_with_reps.loc[
        rna_overlap.obs.index, "pseudo_replicate"
    ].astype(str)

    rna_overlap.obs["CellType_TimePoint_pseudo_replicate"] = (
        rna_overlap.obs["CellType_TimePoint"].astype(str)
        + "."
        + rna_overlap.obs["pseudo_replicate"].astype(str)
    )

    # Generate pseudobulk matrix
    pb = pseudobulk_by_group(rna_overlap, "CellType_TimePoint_pseudo_replicate")
    pb_df = pd.DataFrame(pb.X.T, index=pb.var.index, columns=pb.obs.index)

    # Build metadata table
    metadata = (
        rna_overlap.obs[
            [
                "CellType",
                "TimePoint",
                "CellType_TimePoint_pseudo_replicate",
            ]
        ]
        .drop_duplicates()
        .copy()
    )

    metadata.columns = ["Celltype", "Time_point", "SampleName"]
    metadata["Replicate"] = metadata["SampleName"].str.split(".").str[-1].astype(int) + 1

    metadata = metadata.sort_values(
        by=["Celltype", "Time_point", "Replicate"]
    ).reset_index(drop=True)

    # Reorder columns to match metadata
    pb_df = pb_df[metadata["SampleName"]]

    # Clean output names
    metadata["SampleName"] = sanitize_strings(metadata["SampleName"])
    metadata["Celltype"] = sanitize_strings(metadata["Celltype"])

    pb_df.columns = metadata["SampleName"]
    pb_df.index.name = "GeneID"

    # Save outputs
    pb_df.to_csv(
        "CellType_TimePoint.pseudo_replicate.tsv.gz",
        sep="\t",
        compression="gzip",
    )
    metadata.to_csv(
        "CellType_TimePoint.pseudo_replicate.metadata.tsv",
        sep="\t",
        index=False,
    )

    print("Done.")


if __name__ == "__main__":
    main()