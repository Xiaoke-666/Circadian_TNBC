# Construct TF Connectivity

This pipeline performs transcription factor (TF) network construction across cell classes by integrating chromatin accessibility–derived TF deviation scores with circadian statistics. It also includes downstream steps for weighting, filtering. It also contains an example of trajectory network construction for cancer epithelial cells.

## Overview

The workflow consists of four main steps:

1. Construct TF connectivity network  
2. Weight TF deviation scores using circadian statistics  
3. Filter TFs based on gene expression  
4. Example of construct TF trajectory for cancer epithelial cells  

---

## 1. Construct TF Network

### Input
- Transposed single-cell TF deviation matrix from ChromVar (`.h5ad`)
- Circadian TF statistics for each cell class (`.tsv`)

### Output
- TF connectivity matrix (`.csv`), used for downstream analysis

### Usage
```bash
chmod +x 01_Construct_TF_connectivity.py

./01_Construct_TF_connectivity.py \
  -i ChromVar_Deviations.T.h5ad \
  -p "*.tsv" \
  -o TF_connectivity.csv
```

---

## 2. Weight TF Deviation Scores by Circadian Statistics

### Input
- Single-cell TF deviation matrix from ChromVar (`.h5ad`)
- Circadian TF statistics for each cell class (`.tsv`)
- TF connectivity matrix from Step 1 (`.csv`)

### Output
- TF connectivity matrix with weighted TF deviation scores (`.csv`)

### Usage
```bash
chmod +x 02_Weighted_TF_deviation.R

./02_Weighted_TF_deviation.R \
  -i ChromVar_Deviations.h5ad \
  -n TF_connectivity.csv \
  -d ./circadian_TF/ \
  -o network_weighted_zero_anno_celltype.csv
```

---

## 3. Filter TFs by Gene Expression

### Input
- Motif-to-gene mapping file (`.csv`) with columns:
  - `motif_name`
  - `gene_symbol`
- Pseudobulk gene expression matrix (`.csv.gz`)
- Weighted TF connectivity matrix from Step 2 (`.csv`)

### Output
- Filtered TF network (`.csv`), excluding:
  - TFs with no expression
  - TFs within the lowest 10% of expression levels

### Usage
```bash
chmod +x 03_Filter_TF_by_GEX.R

./03_Filter_TF_by_GEX.R \
  -m motif_mapped.csv \
  -g pseudobulk_GE_mtx.csv.gz \
  -n network_weighted_zero_anno_celltype.csv \
  -o network_motif_weighted_zero_celltype_filtered.csv
```

---

## 4. Construct TF Trajectory for Cancer Epithelial Cells

### Input
- Transposed single-cell TF deviation matrix from ChromVar (`.h5ad`)
- Circadian TF statistics for cancer epithelial cells (`.tsv`)
- Filtered TF list ranked by gene expression (`.csv`)

### Output
- `.h5ad`: TF trajectory with UMAP and PAGA embeddings  
- `.pdf`: Trajectory visualization  

### Usage
```bash
chmod +x 04_Construct_TF_traj_cancer_epi.py

./04_Construct_TF_traj_cancer_epi.py \
  -i ./circadian/ChromVar_Deviations.T.Cancer_epi.h5ad \
  -c ./circadian/Cancer_Epi.discoODAres_CS_ANNOTATION.outer.pvalue_0.02.tsv \
  -t ./circadian/TF_rank_GEX.csv \
  -o ./circadian/Cancer_circadian
```

---

## Notes

- Ensure all input files are properly formatted and consistent across steps.
- The pipeline assumes prior computation of ChromVar deviation scores and circadian statistics.
- File paths can be adjusted as needed for your environment.
