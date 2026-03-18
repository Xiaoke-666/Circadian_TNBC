# Python Downstream Analysis

This directory contains the Python scripts used after sample-level phase inference has been completed.

## Files

- `estimate_donor_internal_phase.py`: summarizes inferred sample phases into donor internal phase (DIP)
- `plot_phase_fits.py`: visualizes expression-phase relationships for selected genes

## Environment

Install Python dependencies from the repository root:

```bash
pip install -r requirements.txt
```

## 1. Estimate donor internal phase

This script reads a result directory containing a `Fit_Output*.csv` file and computes donor-level circadian phase summaries using circular statistics and iterative outlier removal.

Run:

```bash
python src/python/estimate_donor_internal_phase.py results/example_run
```

Optional verbose mode:

```bash
python src/python/estimate_donor_internal_phase.py results/example_run --verbose
```

Main outputs:

- `DIP_results.csv`
- `DIP_summary.png`
- `DIP_circular.png`

## 2. Plot expression against inferred phase

This script overlays inferred phase with expression values for selected genes and fits cosine-like trends for visualization.

Run:

```bash
python src/python/plot_phase_fits.py --results_path results/example_run --expression examples/input_template/expression.csv
```

Optional verbose mode:

```bash
python src/python/plot_phase_fits.py --results_path results/example_run --expression examples/input_template/expression.csv --verbose
```

Main outputs:

- `expression_vs_phase*.png`
- additional gene-level plots depending on the selected input and output naming

## Input assumptions

- Sample IDs in the phase output should match the sample columns in `expression.csv`.
- When metadata or encoded sample names are available, donor and cell-type labels are inferred from those identifiers.
- The plotting script uses the first gene-identifier column in `expression.csv`; `Gene_Symbol` is supported but not required.
- These scripts are designed for the final analysis workflow and prioritize clarity and reproducibility over broad general-purpose packaging.

For a minimal input example, see `examples/input_template/` in the repository root.