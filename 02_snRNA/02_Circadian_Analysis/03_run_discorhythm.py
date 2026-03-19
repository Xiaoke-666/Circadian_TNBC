#!/usr/bin/env python3

from pathlib import Path
import subprocess
import pandas as pd


CELLTYPES = [
    "CD4_T", "CD8_T", "B", "Cancer_Epi", "DC", "Endo", "Macro_Mono",
    "NK", "NKT", "Neu", "Normal_epi", "Plasmablasts", "CAFs", "PVLs"
]

# Project-relative paths
PROJECT_DIR = Path(__file__).resolve().parents[1]
PSEUDOBULK_DIR = PROJECT_DIR / "data" / "pseudobulk"
OUTPUT_BASE_DIR = PROJECT_DIR / "results" / "CircadianGenes"
DISCORHYTHM_SCRIPT = PROJECT_DIR / "scripts" / "run_discorhythm.R"
RSCRIPT_BIN = "Rscript"


def run_discorhythm_for_one_celltype(
    replicate_id: int,
    celltype: str,
    pseudobulk_dir: Path,
    output_base_dir: Path,
    rscript_bin: str,
    discorhythm_script: Path,
) -> None:
    """
    Run DiscoRhythm for one cell type and one pseudobulk replicate.
    """
    metadata_file = pseudobulk_dir / f"CellType_TimePoint.pseudo_replicate_metadata_{replicate_id - 1}.tsv"
    data_file = pseudobulk_dir / f"CellType_TimePoint.pseudo_replicate_{replicate_id - 1}.tsv.gz"

    replicate_dir = output_base_dir / f"Replicate_{replicate_id}"
    celltype_dir = replicate_dir / celltype
    celltype_dir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(metadata_file, sep="\t")
    tmp_metadata = metadata[metadata["Celltype"] == celltype].copy()
    tmp_metadata_file = celltype_dir / f"{celltype}.metadata.tsv"
    tmp_metadata.to_csv(tmp_metadata_file, sep="\t", index=False)

    print(f"[Replicate {replicate_id}] Running DiscoRhythm for {celltype}")

    cmd = [
        rscript_bin,
        str(discorhythm_script),
        "-i", str(celltype_dir),
        "-o", str(celltype_dir),
        "-d", str(data_file),
        "-m", str(tmp_metadata_file),
        "-p", celltype,
    ]

    subprocess.run(cmd, check=True)
    print(f"[Replicate {replicate_id}] {celltype} done")


def main():
    for replicate_id in range(1, 11):
        print(f"\n========== Processing replicate {replicate_id} ==========")
        for celltype in CELLTYPES:
            run_discorhythm_for_one_celltype(
                replicate_id=replicate_id,
                celltype=celltype,
                pseudobulk_dir=PSEUDOBULK_DIR,
                output_base_dir=OUTPUT_BASE_DIR,
                rscript_bin=RSCRIPT_BIN,
                discorhythm_script=DISCORHYTHM_SCRIPT,
            )

    print("\nAll replicates completed.")


if __name__ == "__main__":
    main()