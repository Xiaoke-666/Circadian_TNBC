# Input Template

This directory provides a minimal format template for preparing an input dataset for the circadian phase inference workflow.

It is intended as a schema example only. The files in this directory are toy placeholders and are not meant for biological interpretation.

## Files

- `expression.csv`: gene-by-sample matrix with gene identifiers in the first column
- `metadata.csv`: optional sample metadata table with sample IDs and cell-type labels
- `seed_genes.txt`: optional reference seed list retained with the input package

## Expected conventions

- Sample IDs in `metadata.csv` should match the sample columns in `expression.csv`.
- The first column in `expression.csv` is treated as the gene identifier column by the current workflow.
- Additional metadata columns can be added as needed, but `Sample` and `CellType_D` are the relevant columns for the bundled plotting logic.

## Typical usage

From the repository root:

```bash
bash scripts/run_phase_pipeline.sh examples/input_template
```

The placeholder files in this directory are intentionally minimal. Replace them with the real expression matrix and metadata before running the full analysis.