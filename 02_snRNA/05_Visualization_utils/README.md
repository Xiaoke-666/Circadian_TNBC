# Circadian Visualization Toolkit

This repository provides a collection of R functions for visualizing circadian patterns in transcriptomic and functional enrichment data. The toolkit supports pathway-level analysis, gene-level visualization, and multi-group comparisons.


## Installation

Ensure the following R packages are installed:
`install.packages(c("dplyr", "tidyr", "ggplot2", "stringr"))`

Load the script:
`source("circadian_visualization_utils.R")`

---

## Functions
### 1. Pathway Circadian Curve
Visualizes the average circadian rhythm of a gene set with variability.

### How to run
```bash
plot_pathway_circadian(
  input_file = "data/circadian.tsv", 
  gene_set_file = "data/genes.txt", 
  output_prefix = "output/pathway", 
  plot_title = "Pathway Circadian Pattern"
  )
```
### Output
- CSV of selected genes
- PDF plot of averaged circadian curve

### 2. Functional Enrichment Polar Heatmap
Generates a circular heatmap showing phase-specific enrichment.

### How to run
```bash
plot_polar_enrichment( 
  input_file = "data/enrichment.csv", 
  output_file = "output/polar_heatmap.pdf", custom_order = c("A", "B", "C"), 
  plot_title = "Functional Enrichment"
  )
```

### 3. Acrophase Radial Plot
Displays distribution of acrophase across cell types.

### How to run
```bash
plot_acrophase_radial( 
  input_files = list( CD4_T = "data/CD4_T.tsv", CD8_T = "data/CD8_T.tsv" ), 
  output_file = "output/acrophase.pdf", 
  group_order = c("CD4_T", "CD8_T"), 
  group_colors = c("CD4_T"="#c4222a", "CD8_T"="#e07582") )
```
### 4. Gene-Level Circadian Expression
Plots fitted circadian curve with observed expression.

### How to run
```bash
plot_gene_circadian_expression( 
  circadian_file = "data/circadian.tsv", 
  rna_file = "data/rna.tsv", 
  gene_name = "PER1", 
  subtype_name = "CD4_T", subtype_map = c("CD4_T" = "CD4_T") )
```

## Input Data Requirements
### Circadian result file from DiscoRhythm
Must include:
- GeneID
- mesor
- sincoef
- coscoef
- acrophase
- amplitude
### Enrichment result file from Metascape
Must include:
- Description of pathway
- LogP values
### RNA expression file

