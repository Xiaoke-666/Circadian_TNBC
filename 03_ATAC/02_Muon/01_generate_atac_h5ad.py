#!/usr/bin/env python3

from pathlib import Path
from optparse import OptionParser

import anndata
import muon as mu
import numpy as np
import pandas as pd
from muon import atac as ac


def generate_h5ad(indir, sample, peaks, outdir, **kwargs):
    """
    Generate an AnnData (.h5ad) object for snATAC-seq from fragment files and peak regions.

    Parameters
    ----------
    indir : str
        Input directory containing fragment files.
    sample : str
        Sample name. Expected fragment file: <sample>.atac_fragments.tsv.gz
    peaks : str
        Peak BED/TSV file with at least 3 columns: chromosome, start, end.
    outdir : str
        Output directory for the generated h5ad file.
    """
    indir = Path(indir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fragment_file = indir / f"{sample}.atac_fragments.tsv.gz"
    peak_file = Path(peaks)
    output_file = outdir / f"{sample}.h5ad"

    print("Reading fragments...")
    fragment_df = pd.read_csv(
        fragment_file,
        sep="\t",
        header=None,
        index_col=None,
        engine="c",
        comment="#",
    )
    fragment_df.columns = ["seqnames", "start", "end", "barcode", "count"]

    # Summarize total fragments per barcode
    cell_meta = (
        fragment_df.groupby("barcode")["count"]
        .sum()
        .reset_index()
        .rename(columns={"count": "total"})
    )
    cell_meta["Sample"] = sample
    cell_meta["Sample_barcode"] = cell_meta["Sample"] + "#" + cell_meta["barcode"]
    cell_meta.index = cell_meta["barcode"]
    cell_meta.index.name = None

    print("Finished reading fragments.")

    # Initialize AnnData object
    atac = anndata.AnnData(obs=cell_meta)

    # Read peaks
    print("Reading peak file...")
    peak_df = pd.read_csv(peak_file, sep="\t", header=0)
    peak_df = peak_df.iloc[:, [0, 1, 2]].copy()
    peak_df.columns = ["seqnames", "start", "end"]

    # Keep standard chromosomes only
    valid_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    peak_df = peak_df[peak_df["seqnames"].isin(valid_chroms)].copy()

    peak_df.index = [
        f"{chrom}:{start}-{end}"
        for chrom, start, end in zip(peak_df["seqnames"], peak_df["start"], peak_df["end"])
    ]
    peak_df["interval"] = peak_df.index.astype(str)

    print("Locating fragments...")
    ac.tl.locate_fragments(atac, str(fragment_file))

    print("Counting fragments in peaks...")
    atac = ac.tl.count_fragments_features(
        atac,
        peak_df.rename(
            {"seqnames": "Chromosome", "start": "Start", "end": "End"},
            axis=1,
        ),
        extend_upstream=0,
        extend_downstream=0,
    )

    print("Converting count matrix to int16...")
    atac.X = atac.X.astype(np.int16)

    # Rename obs index
    atac.obs.index = atac.obs["Sample_barcode"]
    atac.obs.index.name = None

    # Basic filtering
    print(f"Before filtering: {atac.n_obs} cells")
    mu.pp.filter_obs(atac, "total", lambda x: (x >= 50) & (x <= 500000))
    print(f"After filtering: {atac.n_obs} cells")

    # Doublet detection placeholder
    print("Skipping doublet calculation (not enabled in this script).")

    print(f"Writing h5ad file to: {output_file}")
    atac.write_h5ad(output_file, compression="gzip", compression_opts=5)

    print("Done!")


if __name__ == "__main__":
    usage = "\nWorkflow for snATAC-seq data analysis: generate h5ad file (based on muon/scanpy)"
    parser = OptionParser(usage, version="Version 1.0")

    parser.add_option(
        "-i", "--indir",
        dest="indir",
        type="string",
        help="Input directory containing fragment files"
    )
    parser.add_option(
        "-s", "--sample",
        dest="sample",
        type="string",
        help="Sample name"
    )
    parser.add_option(
        "-p", "--peak",
        dest="peaks",
        type="string",
        help="Peak file (BED/TSV with chromosome, start, end)",
    )
    parser.add_option(
        "-o", "--outdir",
        dest="outdir",
        type="string",
        help="Output directory"
    )

    options, _ = parser.parse_args()

    required_args = ["indir", "sample", "peaks", "outdir"]
    missing_args = [arg for arg in required_args if getattr(options, arg) is None]

    if missing_args:
        parser.error("Missing required arguments: " + ", ".join(missing_args))

    generate_h5ad(**vars(options))