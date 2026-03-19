# Circadian TF Analysis

This repository provides a simple pipeline to identify circadian transcription factors (TFs) from chromVAR deviation scores.

---

## Workflow

### Step 1. Pseudobulk generation

`01_Pseudobulk.py`

- Input: chromVAR `h5ad`
- Grouping: `CellType × TimePoint × pseudo-replicate`
- Repeated random pseudo-replicates (default: 10)
- Output: pseudobulk TF matrices + metadata

---

### Step 2. Circadian TF detection

`02_run_discorhythm_tf.py`

- Run DiscoRhythm for each cell type and replicate
- Integrate results across replicates
- Identify significant circadian TFs

Output:

- `<celltype>.circadianTF.all.tsv`
- `<celltype>.circadianTF.significant.tsv`

---

### Step 3. TF module / co-activity analysis

`03_Circadian_TF_Modules.py`

- Build TF manifold (PCA / UMAP / Leiden)
- Extract TF–TF connectivity and distance
- Map circadian statistics onto TF network
