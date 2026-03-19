#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(viridis))

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------
parser <- ArgumentParser(
  description = "Plot circadian gene heatmaps for each cell type."
)

parser$add_argument(
  "--matrix_dir",
  required = TRUE,
  help = "Directory containing single-cell circadian expression matrices"
)

parser$add_argument(
  "--circadian_dir",
  required = TRUE,
  help = "Directory containing integrated DiscoRhythm circadian gene results"
)

parser$add_argument(
  "--output_dir",
  required = TRUE,
  help = "Output directory for heatmaps"
)

parser$add_argument(
  "--celltypes",
  nargs = "+",
  default = c("Cancer_Epi", "CAFs", "CD4_T", "B", "Macro_Mono", "CD8_T", "NKT", "Neu"),
  help = "Cell types to plot"
)

args <- parser$parse_args()

matrix_dir <- args$matrix_dir
circadian_dir <- args$circadian_dir
output_dir <- args$output_dir
celltypes <- args$celltypes

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Global settings
# -----------------------------------------------------------------------------
phase_colors <- c(
  "Phase 1" = "#FF4545",
  "Phase 2" = "#FF9C73",
  "Phase 3" = "#FBD288",
  "Phase 4" = "#FCF596",
  "Phase 5" = "#DFF2EB",
  "Phase 6" = "#B9E5E8",
  "Phase 7" = "#7AB2D3",
  "Phase 8" = "#4A628A",
  "Unknown" = "gray"
)

heatmap_color <- colorRamp2(
  seq(-1.5, 1.5, length.out = 11),
  viridis(11, option = "A")
)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
assign_phase <- function(acrophase) {
  case_when(
    acrophase >= 0  & acrophase < 3  ~ "Phase 1",
    acrophase >= 3  & acrophase < 6  ~ "Phase 2",
    acrophase >= 6  & acrophase < 9  ~ "Phase 3",
    acrophase >= 9  & acrophase < 12 ~ "Phase 4",
    acrophase >= 12 & acrophase < 15 ~ "Phase 5",
    acrophase >= 15 & acrophase < 18 ~ "Phase 6",
    acrophase >= 18 & acrophase < 21 ~ "Phase 7",
    acrophase >= 21 & acrophase <= 24 ~ "Phase 8",
    TRUE ~ "Unknown"
  )
}

read_and_smooth_matrix <- function(file_path) {
  message("Reading matrix: ", file_path)

  mat <- fread(file_path, header = TRUE, data.table = FALSE, sep = "\t")
  rownames(mat) <- mat[[1]]
  mat <- mat[, -1, drop = FALSE]
  mat <- na.omit(mat)

  smoothed <- apply(mat, 1, function(x) smooth.spline(x, df = 3)$y)
  zscore <- apply(smoothed, 1, function(x) (x - mean(x)) / sd(x))
  zscore <- as.data.frame(zscore)
  rownames(zscore) <- rownames(mat)

  return(zscore)
}

read_circadian_genes <- function(file_path, pvalue_cutoff = 0.02) {
  message("Reading circadian gene table: ", file_path)

  genes <- read.table(
    file_path,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
  )

  genes <- genes[genes$pvalue < pvalue_cutoff, , drop = FALSE]
  genes <- genes[order(genes$acrophase), , drop = FALSE]
  genes$Phase <- assign_phase(genes$acrophase)

  return(genes)
}

plot_heatmap <- function(pt_matrix, phase_vector, output_file) {
  row_anno <- rowAnnotation(
    Phase = phase_vector,
    col = list(Phase = phase_colors),
    show_annotation_name = FALSE,
    width = unit(3, "mm")
  )

  ht <- Heatmap(
    pt_matrix,
    name = "z-score",
    col = heatmap_color,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 6),
    row_title_rot = 0,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    row_split = phase_vector,
    use_raster = TRUE,
    row_title = NULL
  )

  png(output_file, width = 4, height = 7, units = "in", res = 300)
  draw(row_anno + ht, heatmap_legend_side = "right")
  dev.off()
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
for (celltype in celltypes) {
  message("==================================================")
  message("Processing cell type: ", celltype)

  matrix_file <- file.path(
    matrix_dir,
    paste0("rna4.Cellbender.", celltype, ".CircadianGenes.ZT4_10_16_22.top.sorted.tsv.gz")
  )

  circadian_file <- file.path(
    circadian_dir,
    paste0(celltype, ".discoODAres_CS_ANNOATION.outer.tsv")
  )

  pt_matrix <- read_and_smooth_matrix(matrix_file)
  circadian_genes <- read_circadian_genes(circadian_file, pvalue_cutoff = 0.02)

  common_genes <- intersect(rownames(pt_matrix), rownames(circadian_genes))
  pt_matrix <- pt_matrix[common_genes, , drop = FALSE]
  row_split <- circadian_genes[common_genes, "Phase", drop = TRUE]

  n_cells <- ncol(pt_matrix) / 4
  message("Number of cells per time block: ", n_cells)
  message("Plotting heatmap...")

  output_file <- file.path(
    output_dir,
    paste0(celltype, ".CircadianGenes.heatmap.", n_cells, "cells.smooth.spline.png")
  )

  plot_heatmap(
    pt_matrix = pt_matrix,
    phase_vector = row_split,
    output_file = output_file
  )
}

message("Done.")
