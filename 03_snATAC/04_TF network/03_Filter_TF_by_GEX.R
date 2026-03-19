#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(argparse)
})

compute_tf_motif_weighting <- function(
  motif_map_file,
  gene_exp_file,
  tf_network_file,
  output_file,
  top_fraction = 0.9
) {

  # ---------------------------
  # Load motif mapping
  # ---------------------------
  motif_df <- read.csv(motif_map_file)

  required_cols <- c("motif_name", "gene_symbol")
  if (!all(required_cols %in% colnames(motif_df))) {
    stop("motif_map_file must contain: motif_name, gene_symbol")
  }

  # Fill missing gene_symbol with motif_name
  motif_df$gene_symbol[motif_df$gene_symbol == ""] <- motif_df$motif_name[motif_df$gene_symbol == ""]

  # ---------------------------
  # Load gene expression
  # ---------------------------
  geneExp <- as.data.frame(fread(file = gene_exp_file, header = TRUE))

  if (ncol(geneExp) < 2) {
    stop("Gene expression matrix must have at least 2 columns")
  }

  colnames(geneExp)[1] <- "gene_symbol"

  # Remap unmatched gene symbols
  motif_df$gene_symbol[!motif_df$gene_symbol %in% geneExp$gene_symbol] <- 
    motif_df$motif_name[!motif_df$gene_symbol %in% geneExp$gene_symbol]

  # Filter expression
  geneExp <- geneExp[geneExp$gene_symbol %in% motif_df$gene_symbol, ]

  # ---------------------------
  # Merge motif + expression
  # ---------------------------
  motif_df <- motif_df %>%
    left_join(geneExp, by = "gene_symbol") %>%
    na.omit()

  # ---------------------------
  # Build TF expression ranking
  # ---------------------------
  df_TF_exp <- data.frame()

  expr_cols <- colnames(motif_df)[!(colnames(motif_df) %in% c("motif_name","gene_symbol"))]

  for (celltype in expr_cols) {
    df_tmp <- motif_df[, c("motif_name","gene_symbol",celltype)]

    df_tmp <- df_tmp %>%
      arrange(desc(.data[[celltype]])) %>%
      filter(.data[[celltype]] != 0)

    if (nrow(df_tmp) == 0) next

    cutoff <- round(top_fraction * nrow(df_tmp))
    df_tmp <- df_tmp[1:cutoff, ]

    colnames(df_tmp)[3] <- "Expression"
    df_tmp$celltype <- celltype

    df_TF_exp <- rbind(df_TF_exp, df_tmp)
  }

  # ---------------------------
  # Load TF network (weighted)
  # ---------------------------
  TF_network_full <- read.csv(tf_network_file)

  TF_network <- unique(
    TF_network_full[, c("source", grep("weighted", colnames(TF_network_full), value = TRUE))]
  )

  TF_network <- TF_network[TF_network$source %in% unique(df_TF_exp$motif_name), ]

  rownames(TF_network) <- TF_network$source
  TF_network <- TF_network[, -1, drop = FALSE]
  colnames(TF_network) <- gsub("_weighted", "", colnames(TF_network))

  df_TF_exp <- df_TF_exp[df_TF_exp$motif_name %in% rownames(TF_network), ]

  # ---------------------------
  # Convert to wide format
  # ---------------------------
  df_TF_exp_wide <- df_TF_exp %>%
    pivot_wider(
      names_from = celltype,
      values_from = Expression,
      values_fill = 0,
      values_fn = mean
    ) %>%
    column_to_rownames("motif_name")

  # Align matrices
  common_rows <- intersect(rownames(TF_network), rownames(df_TF_exp_wide))
  common_cols <- intersect(colnames(TF_network), colnames(df_TF_exp_wide))

  TF_network <- TF_network[common_rows, common_cols, drop = FALSE]
  df_TF_exp_wide <- df_TF_exp_wide[common_rows, common_cols, drop = FALSE]

  # Binarize expression
  df_TF_exp_wide[df_TF_exp_wide > 0] <- 1

  # ---------------------------
  # Apply weighting
  # ---------------------------
  g_w_TF_network <- TF_network * df_TF_exp_wide

  # Remove TFs with all zeros
  g_w_TF_network <- g_w_TF_network[rowSums(g_w_TF_network) > 0, , drop = FALSE]

  # ---------------------------
  # Map back to edge list
  # ---------------------------
  TF_edges <- TF_network_full[, c("source","target","connectivity","distance")]

  TF_edges <- TF_edges[
    TF_edges$source %in% rownames(g_w_TF_network) &
    TF_edges$target %in% rownames(g_w_TF_network),
  ]

  g_w_TF_network <- t(g_w_TF_network)
  rownames(g_w_TF_network) <- paste0(rownames(g_w_TF_network), "_weighted")

  # Add scores
  for (ct in rownames(g_w_TF_network)) {
    tf_scores <- setNames(as.numeric(g_w_TF_network[ct, ]), colnames(g_w_TF_network))
    TF_edges[[ct]] <- tf_scores[TF_edges$source]
  }

  # ---------------------------
  # Save output
  # ---------------------------
  write.csv(TF_edges, file = output_file, quote = FALSE, row.names = FALSE)

  message("Output saved to: ", output_file)
}

# ---------------------------
# CLI interface
# ---------------------------
main <- function() {
  parser <- ArgumentParser(description = "Motif-based TF network weighting")

  parser$add_argument("-m", "--motif_map", required = TRUE, help = "Motif mapping CSV")
  parser$add_argument("-g", "--gene_exp", required = TRUE, help = "Gene expression matrix (csv.gz)")
  parser$add_argument("-n", "--network", required = TRUE, help = "TF network CSV")
  parser$add_argument("-o", "--output", default = "network_motif_weighted_zero_celltype_filtered.csv")
  parser$add_argument("--top_fraction", type = "double", default = 0.9)

  args <- parser$parse_args()

  compute_tf_motif_weighting(
    motif_map_file = args$motif_map,
    gene_exp_file = args$gene_exp,
    tf_network_file = args$network,
    output_file = args$output,
    top_fraction = args$top_fraction
  )
}

if (sys.nframe() == 0) {
  main()
}