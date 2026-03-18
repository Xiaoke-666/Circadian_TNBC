# Run inferCNV and calculate CNV score

This pipeline runs [InferCNV] (https://github.com/broadinstitute/inferCNV/wiki) and calculate the CNV score based on its output

## Overview

The workflow consists of two main steps:

1. Run InferCNV pipelines
2. Calculate and plot CNV score based on InferCNV results 

---

## 1. Run InferCNV pipelines

### Input
- A `Seurat` object in `.rds` format
- A gene/chromosome positions file for mouse (avaliable at https://data.broadinstitute.org/Trinity/CTAT/cnv/)

### Output
- An output directory from InferCNV. More details for outfiles can be found: https://github.com/broadinstitute/inferCNV/wiki/Output-Files

### How to Run
Use the code in `01_Run_infercnv.R` and modify the parameters based on the guidance from InfeCNV (https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV). 

---

## 2. Calculate CNV score

### Input
- InferCNV output (`run.final.infercnv_obj.rds`)
- A `Seurat` object contains UMAP embeddings which will be used for plottin (`.rds`). 

### Output
- Calculated InferCNV score for each cell (`.csv`)
- A `Seurat` object with calculated InferCNV in `meta.data`
- UMAP Plot for InferCNV score (`.pdf`)

### Usage
```bash
chmod +x 02_Infercnv_score.R

./02_Infercnv_score.R \
  -i ./infercnv_output/run.final.infercnv_obj \
  -s ./raw_data.rds \
  --reduction Xumap_ \
  -o ./figures/infercnv
```


## Notes

- Ensure all input files are properly formatted and consistent across steps.
- The pipeline assumes prior computation of ChromVar deviation scores and circadian statistics.
- File paths can be adjusted as needed for your environment.
