#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(schard)
  library(tidyverse)
  library(argparse)
})

compute_tf_network_celltype <- function(
  h5ad_file,
  tf_network_file,
  integration_dir,
  output_file,
  exclude_celltypes = c("Normal_Epi","Endo","PVLs","Plasmablasts","DC","NK")
) {

  # ---------------------------
  # Load Seurat object
  # ---------------------------
  tf_score <- h5ad2seurat(h5ad_file)

  meta_col <- "ATAC.Subclass.June18_2025"

  if (!meta_col %in% colnames(tf_score@meta.data)) {
    stop(paste("Metadata column not found:", meta_col))
  }

  # Filter cells
  valid_cells <- rownames(tf_score@meta.data)[
    !tf_score@meta.data[[meta_col]] %in% c(exclude_celltypes, NA)
  ]

  sub <- subset(tf_score, cells = valid_cells)

  # ---------------------------
  # Build expression matrix
  # ---------------------------
  mat <- as.data.frame(t(as.matrix(sub[["RNA"]]$data)))

  mat$cell_type <- sub@meta.data[
    match(rownames(mat), rownames(sub@meta.data)),
    meta_col
  ]

  mat$cell_type <- as.character(mat$cell_type)

  # ---------------------------
  # Pseudobulk (mean per cell type)
  # ---------------------------
  num_cols <- setdiff(colnames(mat), "cell_type")

  pseudobulk <- mat %>%
    group_by(cell_type) %>%
    summarise(across(all_of(num_cols), mean), .groups = "drop") %>%
    column_to_rownames("cell_type")

  # ---------------------------
  # Load TF network
  # ---------------------------
  TF_net <- read.csv(tf_network_file)
  TFs <- unique(TF_net$source)

  # Subset pseudobulk to TFs
  pseudobulk <- as.data.frame(t(as.matrix(pseudobulk[, intersect(TFs, colnames(pseudobulk)), drop = FALSE])))

  # ---------------------------
  # Load circadian TF files
  # ---------------------------
  files <- list.files(circadian_TF, pattern = "\\.tsv$", full.names = TRUE)

  if (length(files) == 0) {
    stop("No TSV files found in circadian TF directory.")
  }

  weight <- data.frame(TF = rownames(pseudobulk))

  for (i in seq_along(files)) {
    name <- sub("\\..*", "", basename(files[i]))
    df <- read.table(files[i], sep = "\t", header = TRUE)

    if (!all(c("GeneID", "pvalue") %in% colnames(df))) {
      stop(paste("Missing required columns in:", files[i]))
    }

    df$log1_over_p <- log(1 / df$pvalue)

    weight[[name]] <- df$log1_over_p[
      match(weight$TF, df$GeneID)
    ]
  }

  weight[is.na(weight)] <- 0
  rownames(weight) <- weight$TF
  weight <- weight[, -1, drop = FALSE]

  # Align matrices
  weight <- weight[match(rownames(pseudobulk), rownames(weight)), , drop = FALSE]
  pseudobulk <- pseudobulk[, match(colnames(weight), colnames(pseudobulk)), drop = FALSE]

  # ---------------------------
  # Apply weighting
  # ---------------------------
  pseudobulk_weight <- pseudobulk * weight
  colnames(pseudobulk_weight) <- paste0(colnames(pseudobulk_weight), "_weighted")

  # Combine
  pseudo_all <- rbind(
    t(pseudobulk),
    t(pseudobulk_weight)
  )

  # ---------------------------
  # Add TF scores to network
  # ---------------------------
  for (ct in rownames(pseudo_all)) {
    tf_scores <- setNames(as.numeric(pseudo_all[ct, ]), colnames(pseudo_all))
    TF_net[[ct]] <- tf_scores[TF_net$source]
  }

  # ---------------------------
  # Remove ENS genes
  # ---------------------------
  TF_net <- TF_net[
    !grepl("ENS", TF_net$source) &
    !grepl("ENS", TF_net$target),
  ]

  # ---------------------------
  # Save output
  # ---------------------------
  write.csv(TF_net, file = output_file, quote = FALSE, row.names = FALSE)

  message("Output saved to: ", output_file)
}

# ---------------------------
# CLI interface
# ---------------------------
main <- function() {
  parser <- ArgumentParser(description = "Compute TF network with cell-type pseudobulk weighting")

  parser$add_argument("-i", "--input", required = TRUE, help = "Input h5ad file")
  parser$add_argument("-n", "--network", required = TRUE, help = "TF connectivity CSV")
  parser$add_argument("-d", "--dir", required = TRUE, help = "Directory with circadian TSV files")
  parser$add_argument("-o", "--output", default = "network_weighted_zero_anno_celltype.csv")

  args <- parser$parse_args()

  compute_tf_network_celltype(
    h5ad_file = args$input,
    tf_network_file = args$network,
    circadian_TF = args$dir,
    output_file = args$output
  )
}

if (sys.nframe() == 0) {
  main()
}