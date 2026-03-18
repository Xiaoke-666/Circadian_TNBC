"""
Plot core circadian gene expression vs inferred phase from a Fit_Output_*.csv. 
"""

import argparse
from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def log(message, verbose=False):
    if verbose:
        print(message)


def get_gene_column(df):
    for candidate in ["Gene_Symbol", "gene", "Gene", "symbol"]:
        if candidate in df.columns:
            return candidate

    for column in df.columns:
        if "gene" in str(column).lower():
            return column

    return df.columns[0]


def get_sample_column(df):
    for candidate in ["ID", "Sample_ID", "Sample"]:
        if candidate in df.columns:
            return candidate
    raise ValueError("Fit_Output file must contain one of: ID, Sample_ID, Sample.")


def get_phase_column(df):
    for candidate in ["Phase_MA", "Phases_MA", "Phase", "Predicted_Phase_Hours"]:
        if candidate in df.columns:
            return candidate
    raise ValueError(
        "Fit_Output file must contain one of: Phase_MA, Phases_MA, Phase, Predicted_Phase_Hours."
    )


def phase_series_to_hours(series, phase_col):
    values = pd.to_numeric(series, errors="coerce")
    if phase_col == "Predicted_Phase_Hours":
        return values % 24.0
    return (values % (2 * np.pi)) * (24.0 / (2 * np.pi))


def candidate_dataset_names(result_dataset):
    names = []

    def add(value):
        if value and value not in names:
            names.append(value)

    add(result_dataset)
    add(result_dataset.removeprefix("data_"))

    stripped_timestamp = re.sub(r"_\d{8}_\d{6}$", "", result_dataset)
    add(stripped_timestamp)
    add(stripped_timestamp.removeprefix("data_"))

    return names


DEFAULT_GENES = [
    "ARNTL", "CLOCK", "NPAS2", "NR1D1", "BHLHE41", "NR1D2", "DBP",
    "CIART", "PER1", "PER2", "PER3", "TEF", "HLF", "CRY1", "CRY2",
    "RORC", "NFIL3"
]


def find_repo_root(start_path: Path) -> Path:
    """Walk upward from a results path until the repository root is found."""

    for candidate in [start_path, *start_path.parents]:
        if (candidate / "src").exists() and (
            (candidate / "scripts").exists() or (candidate / "workflows").exists()
        ):
            return candidate
    raise ValueError(
        "Could not infer repository root from the results path. "
        "Please pass --expression explicitly."
    )


def find_fit_file(base: Path) -> Path:
    """Given a CYCLOPS results directory, pick a Fit_Output*.csv file.

    If base itself is a Fit_Output*.csv, just return it.
    Otherwise, search under base for files starting with "Fit_Output" and
    choose the lexicographically latest one (usually the most recent).
    Supports both Fit_Output.csv and Fit_Output_*.csv formats.
    """

    base = base.resolve()
    if base.is_file() and base.name.startswith("Fit_Output") and base.suffix == ".csv":
        return base

    if not base.is_dir():
        raise FileNotFoundError(
            "Results path was not found or is not a directory. "
            "Pass a results directory or a Fit_Output*.csv file."
        )

    # Try both Fit_Output_*.csv and Fit_Output.csv
    candidates = sorted(
        p for p in base.glob("Fit_Output*.csv") if p.is_file()
    )
    if not candidates:
        # Allow searching one level down in case user points to higher directory
        candidates = sorted(
            p for p in base.rglob("Fit_Output*.csv") if p.is_file()
        )
    if not candidates:
        raise FileNotFoundError(
            "No Fit_Output*.csv file was found in the provided results location."
        )
    return candidates[-1]


def map_fit_to_expression(fit_path: Path) -> Path:
    """Map a Fit_Output path to the corresponding expression.csv."""

    fit_path = fit_path.resolve()

    parts = list(fit_path.parts)
    if "results" not in parts:
        raise ValueError(
            "Fit path does not contain 'results'; cannot infer dataset name. "
            "Please pass --expression explicitly."
        )

    project_root = find_repo_root(fit_path.parent)
    res_idx = parts.index("results")
    if res_idx + 1 >= len(parts):
        raise ValueError("Cannot determine dataset directory from fit path")

    result_dataset = parts[res_idx + 1]
    data_root = project_root / "data"

    tried = []
    for dataset_name in candidate_dataset_names(result_dataset):
        expr_candidate = data_root / dataset_name / "expression.csv"
        tried.append(f"data/{dataset_name}/expression.csv")
        if expr_candidate.exists():
            return expr_candidate

    for expr_candidate in data_root.rglob("expression.csv"):
        relative_parent = expr_candidate.parent.relative_to(data_root)
        flattened = "_".join(relative_parent.parts)
        if flattened in candidate_dataset_names(result_dataset):
            return expr_candidate

    raise FileNotFoundError(
        "Could not map Fit_Output to expression.csv from the repository layout. Tried: "
        + ", ".join(tried)
        + ". Use --expression to specify it."
    )


def get_celltype_mapping(expr_file: Path, sample_ids):
    """Extract celltype information from Celltype_D row if it exists.
    
    Returns:
        dict mapping sample_id -> celltype, or None if Celltype_D row doesn't exist
    """
    df = pd.read_csv(expr_file, low_memory=False)
    gene_col = get_gene_column(df)
    
    # Some datasets encode cell type as a pseudo-row in expression.csv rather
    # than in a separate metadata table. This helper supports that layout.
    celltype_rows = df[df[gene_col].astype(str).str.upper() == "CELLTYPE_D"]
    if celltype_rows.empty:
        return None
    
    celltype_row = celltype_rows.iloc[0]
    celltype_map = {}
    for sid in sample_ids:
        if sid in celltype_row.index:
            celltype_map[sid] = celltype_row[sid]
    
    return celltype_map if celltype_map else None


def load_expression_for_genes(expr_file: Path, genes_upper, sample_ids):
    """Load expression values for given genes and samples.

    - expr_file: path to expression.csv
    - genes_upper: list of uppercase gene symbols
    - sample_ids: ordered sample IDs to match the phases from Fit_Output
    Returns a DataFrame with index=sample_ids and columns=genes_upper.
    """

    df = pd.read_csv(expr_file, low_memory=False)
    gene_col = get_gene_column(df)

    df["GENE_UP"] = df[gene_col].astype(str).str.upper()
    
    # Exclude Celltype_D row from gene expression data (case-insensitive)
    df = df[df["GENE_UP"] != "CELLTYPE_D"]

    # Fit output defines the sample order used for phase values, so expression
    # columns are reordered to match those sample IDs before plotting.
    sample_cols = [sid for sid in sample_ids if sid in df.columns]
    if not sample_cols:
        raise ValueError(
            "None of the Fit_Output sample IDs were found as columns in "
            f"{expr_file}. Check that dataset and results correspond."
        )

    for c in sample_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    grp = df.groupby("GENE_UP")[sample_cols].mean()

    data = {}
    for g in genes_upper:
        if g in grp.index:
            row = grp.loc[g].reindex(sample_ids).values
        else:
            row = np.full(len(sample_ids), np.nan)
        data[g] = row

    return pd.DataFrame(data, index=sample_ids)


def fit_cosine(times_hours, y):
    """Fit a simple 24h cosine: y ~ a + b*cos(w t) + c*sin(w t)."""

    mask = ~np.isnan(y) & np.isfinite(times_hours)
    if mask.sum() < 3:
        return np.nan, np.nan, np.nan  # a, amplitude, phase_rad

    x = times_hours[mask]
    v = y[mask]
    theta = x * (2 * np.pi / 24.0)
    A = np.column_stack([np.ones_like(theta), np.cos(theta), np.sin(theta)])
    try:
        coeffs, *_ = np.linalg.lstsq(A, v, rcond=None)
    except Exception:
        return np.nan, np.nan, np.nan

    a, b, c = coeffs
    amp = float(np.hypot(b, c))
    phase = float(np.arctan2(c, b) % (2 * np.pi))
    return float(a), amp, phase


def plot_expression_vs_phase(expr_by_gene, phases_hours, genes, out_png, verbose=False):
    """Scatter + cosine fit for each gene across 0–24h phase."""

    genes = list(genes)
    n = len(genes)
    if n == 0:
        log("No genes to plot.", verbose)
        return

    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 2.8 * nrows))
    axes = np.array(axes).reshape(-1)

    x_grid = np.linspace(0, 24, 300)

    fit_records = []

    for i, g in enumerate(genes):
        ax = axes[i]
        if g not in expr_by_gene.columns:
            ax.text(0.5, 0.5, "missing", ha="center", transform=ax.transAxes)
            ax.set_title(g)
            ax.set_xlabel("Phase (hours)")
            continue

        y = expr_by_gene[g].values.astype(float)
        mask = ~np.isnan(y) & np.isfinite(phases_hours)
        if mask.sum() == 0:
            ax.text(0.5, 0.5, "no data", ha="center", transform=ax.transAxes)
            ax.set_title(g)
            ax.set_xlabel("Phase (hours)")
            continue

        x = phases_hours[mask]
        v = y[mask]

        # The fit is estimated on the full ordered vector so that missing
        # expression values are handled consistently with the phase ordering.
        a, amp, phase_rad = fit_cosine(phases_hours, y)

        ax.scatter(x, v, s=18, alpha=0.8)

        if np.isfinite(amp):
            theta_grid = x_grid * (2 * np.pi / 24.0)
            b = amp * np.cos(phase_rad)
            c = amp * np.sin(phase_rad)
            y_fit = a + b * np.cos(theta_grid) + c * np.sin(theta_grid)
            ax.plot(x_grid, y_fit, color="C1", linewidth=1.5)
            title = f"{g} (amp={amp:.2f}, φ={phase_rad:.2f})"
            fit_records.append({"Gene": g, "amplitude": amp, "phase_rad": phase_rad})
        else:
            title = f"{g} (fit failed)"

        ax.set_title(title)
        ax.set_xlabel("Phase (hours)")

    # hide unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"Generated: {out_png}")

    return pd.DataFrame(fit_records)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot core circadian gene expression vs inferred phase "
            "from a Fit_Output_*.csv."
        )
    )
    parser.add_argument(
        "--results_path",
        required=True,
        help=(
            "Path to a CYCLOPS results directory (containing Fit_Output_*.csv) "
            "or to a specific Fit_Output_*.csv file."
        ),
    )
    parser.add_argument(
        "--expression",
        help=(
            "Optional path to expression.csv. If omitted, the script will try "
            "to infer it from the repository layout."
        ),
    )
    parser.add_argument(
        "--genes",
        nargs="*",
        default=None,
        help=(
            "Genes to plot (default = canonical clock genes). Gene symbols "
            "are matched case-insensitively."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print additional details about file selection and cell-type grouping.",
    )

    args = parser.parse_args()
    verbose = args.verbose

    genes = args.genes if args.genes else DEFAULT_GENES
    genes_upper = [g.upper() for g in genes]

    results_path = Path(args.results_path)
    fit_file = find_fit_file(results_path)
    log(f"Using Fit file: {fit_file.name}", verbose)

    if args.expression:
        expr_file = Path(args.expression)
    else:
        expr_file = map_fit_to_expression(fit_file)
    log(f"Using expression: {expr_file.name}", verbose)

    fit_df = pd.read_csv(fit_file)
    sample_col = get_sample_column(fit_df)
    phase_col = get_phase_column(fit_df)
    log(f"Using {sample_col} and {phase_col} from fit output", verbose)

    sample_ids = fit_df[sample_col].astype(str).tolist()
    phases_hours = phase_series_to_hours(fit_df[phase_col], phase_col).values

    # If cell-type labels are available, generate one plot per cell type so the
    # manuscript figures can be compared at the same biological resolution.
    celltype_map = get_celltype_mapping(expr_file, sample_ids)
    
    if celltype_map:
        # Group samples by celltype
        celltype_groups = {}
        for sid, ctype in celltype_map.items():
            if ctype not in celltype_groups:
                celltype_groups[ctype] = []
            celltype_groups[ctype].append(sid)
        
        log(f"Found {len(celltype_groups)} cell types: {list(celltype_groups.keys())}", verbose)
        
        # Plot for each celltype
        for ctype, ctype_samples in celltype_groups.items():
            # Get indices of these samples in the original sample_ids list
            indices = [i for i, sid in enumerate(sample_ids) if sid in ctype_samples]
            if not indices:
                continue
            
            # Filter expression and phases for this celltype
            expr_by_gene = load_expression_for_genes(expr_file, genes_upper, sample_ids)
            expr_filtered = expr_by_gene.iloc[indices]
            phases_filtered = phases_hours[indices]
            
            out_png = fit_file.with_name(f"expression_vs_phase_{ctype}.png")
            log(f"Plotting {ctype} with {len(indices)} samples...", verbose)
            plot_expression_vs_phase(expr_filtered, phases_filtered, genes_upper, out_png, verbose=verbose)
    else:
        # No Celltype_D found, plot all samples together
        log("No Celltype_D found, plotting all samples together", verbose)
        expr_by_gene = load_expression_for_genes(expr_file, genes_upper, sample_ids)
        out_png = fit_file.with_name("expression_vs_phase.png")
        plot_expression_vs_phase(expr_by_gene, phases_hours, genes_upper, out_png, verbose=verbose)


if __name__ == "__main__":
    main()
