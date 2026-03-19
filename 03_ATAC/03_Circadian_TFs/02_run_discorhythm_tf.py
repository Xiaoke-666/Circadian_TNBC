#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import subprocess

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_DIR / "results" / "circadian_tf"
OUTPUT_DIR = DATA_DIR / "CircadianTF_Subclass"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DISCORHYTHM_SCRIPT = PROJECT_DIR / "scripts" / "02_run_discorhythm.R"

N_REPEATS = 10

CELLTYPES = [
    "CD4_T", "CD8_T", "B",
    "Cancer_Epi", "DC", "Endo", "Macro_Mono", "NK", "NKT", "Neu",
    "Normal_epi", "Plasmablasts", "CAFs", "PVLs"
]

# -----------------------------------------------------------------------------
# Run DiscoRhythm for each cell type and replicate
# -----------------------------------------------------------------------------
def run_discorhythm_for_celltype(celltype, replicate_idx):
    print(f"\n[INFO] Processing {celltype} | replicate {replicate_idx}")

    metadata_file = DATA_DIR / f"CellType_TimePoint.pseudo_replicate_metadata_{replicate_idx}.tsv"
    matrix_file = DATA_DIR / f"CellType_TimePoint.pseudo_replicate_{replicate_idx}.tsv.gz"

    metadata = pd.read_csv(metadata_file, sep="\t")

    # subset metadata for one cell type
    metadata_ct = metadata[metadata["Celltype"] == celltype]

    out_dir = OUTPUT_DIR / f"Replicate_{replicate_idx+1}" / celltype
    out_dir.mkdir(parents=True, exist_ok=True)

    metadata_ct_file = out_dir / f"{celltype}.metadata.tsv"
    metadata_ct.to_csv(metadata_ct_file, sep="\t", index=False)

    # run R script (DiscoRhythm)
    cmd = [
        "Rscript",
        str(DISCORHYTHM_SCRIPT),
        "-i", str(out_dir),
        "-o", str(out_dir),
        "-d", str(matrix_file),
        "-m", str(metadata_ct_file),
        "-p", celltype
    ]

    print("[INFO] Running DiscoRhythm...")
    subprocess.run(cmd, check=True)

    print(f"[INFO] Done: {celltype} | replicate {replicate_idx}")


# -----------------------------------------------------------------------------
# Integration across replicates
# -----------------------------------------------------------------------------
def integrate_results(celltype):
    print(f"\n[INFO] Integrating results for {celltype}")

    all_df = None

    for i in range(1, N_REPEATS):
        file = OUTPUT_DIR / f"Replicate_{i}" / celltype / f"{celltype}.discoODAres_CS_ANNOATION.tsv"

        tmp = pd.read_csv(file, sep="\t")

        tmp = tmp[
            ["GeneID", "acrophase", "amplitude", "pvalue", "qvalue", "Rsq", "mesor", "sincoef", "coscoef"]
        ]

        tmp.columns = [
            "GeneID",
            f"acrophase_{i}",
            f"amplitude_{i}",
            f"pvalue_{i}",
            f"qvalue_{i}",
            f"Rsq_{i}",
            f"mesor_{i}",
            f"sincoef_{i}",
            f"coscoef_{i}",
        ]

        if all_df is None:
            all_df = tmp
        else:
            all_df = pd.merge(all_df, tmp, on="GeneID", how="outer")

    # -------------------------------------------------------------------------
    # Aggregate statistics
    # -------------------------------------------------------------------------
    idx = list(range(1, N_REPEATS))

    all_df["acrophase"] = all_df[[f"acrophase_{i}" for i in idx]].median(axis=1)
    all_df["amplitude"] = all_df[[f"amplitude_{i}" for i in idx]].median(axis=1)

    all_df["pvalue"] = all_df[[f"pvalue_{i}" for i in idx]].mean(axis=1)
    all_df["qvalue"] = all_df[[f"qvalue_{i}" for i in idx]].mean(axis=1)
    all_df["Rsq"] = all_df[[f"Rsq_{i}" for i in idx]].mean(axis=1)
    all_df["mesor"] = all_df[[f"mesor_{i}" for i in idx]].mean(axis=1)
    all_df["sincoef"] = all_df[[f"sincoef_{i}" for i in idx]].mean(axis=1)
    all_df["coscoef"] = all_df[[f"coscoef_{i}" for i in idx]].mean(axis=1)

    # -------------------------------------------------------------------------
    # Save full results
    # -------------------------------------------------------------------------
    out_all = OUTPUT_DIR / f"{celltype}.circadianTF.all.tsv"
    all_df.to_csv(out_all, sep="\t", index=False)

    # -------------------------------------------------------------------------
    # Filter significant TFs
    # -------------------------------------------------------------------------
    sig = all_df[all_df["pvalue"] < 0.02].copy()
    sig = sig.sort_values("acrophase")

    out_sig = OUTPUT_DIR / f"{celltype}.circadianTF.significant.tsv"
    sig.to_csv(out_sig, sep="\t", index=False)

    print(f"[INFO] Saved: {celltype}")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():

    # Step 1: Run DiscoRhythm
    for celltype in CELLTYPES:
        for i in range(N_REPEATS):
            run_discorhythm_for_celltype(celltype, i)

    # Step 2: Integration
    integrate_celltypes = [
        "CD4_T", "CD8_T", "B",
        "Cancer_Epi", "Macro_Mono", "NKT", "Neu", "CAFs"
    ]

    for celltype in integrate_celltypes:
        integrate_results(celltype)


if __name__ == "__main__":
    main()