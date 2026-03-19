#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(Seurat)
  library(scCustomize)
  library(RColorBrewer)
  library(ggplot2)
  library(cowplot)
  library(argparse)
})

infercnv_score <- function(
  infercnv_rds,
  seurat_rds,
  output_prefix,
  reduction = "Xumap_"
) {

  # ---------------------------
  # Load inferCNV object
  # ---------------------------
  infercnv_obj <- readRDS(infercnv_rds)

  if (is.null(infercnv_obj@expr.data)) {
    stop("infercnv_obj does not contain expr.data")
  }

  exp <- infercnv_obj@expr.data

  # ---------------------------
  # Rescale (-1 to 1) and compute CNV score
  # ---------------------------
  exp <- scale(t(exp))

  tmp1 <- sweep(exp, 2, apply(exp, 2, min), "-")
  tmp2 <- apply(exp, 2, max) - apply(exp, 2, min)

  exp_1 <- t(2 * sweep(tmp1, 2, tmp2, "/") - 1)

  cnv_score <- as.data.frame(colSums(exp_1 * exp_1))
  colnames(cnv_score) <- "infercnv_score"

  # ---------------------------
  # Load Seurat object
  # ---------------------------
  seurat_obj <- readRDS(seurat_rds)

  if (!"ATAC.Subclass.June18_2025" %in% colnames(seurat_obj@meta.data)) {
    stop("Missing metadata column: ATAC.Subclass.June18_2025")
  }

  if (!"CD45" %in% colnames(seurat_obj@meta.data)) {
    stop("Missing metadata column: CD45")
  }

  Idents(seurat_obj) <- seurat_obj$ATAC.Subclass.June18_2025

  # Subset CD45-
  seurat_obj <- subset(seurat_obj, subset = CD45 == "CD45-")

  metadata <- seurat_obj@meta.data

  # Align CNV scores
  cnv_score <- cnv_score[
    match(rownames(metadata), rownames(cnv_score)),
    , drop = FALSE
  ]

  metadata$infercnv_score <- cnv_score$infercnv_score
  seurat_obj@meta.data <- metadata

  # Save updated Seurat
  out_rds <- paste0(output_prefix, "_CD45_subset.rds")
  saveRDS(seurat_obj, file = out_rds)

  # ---------------------------
  # Plotting
  # ---------------------------
  p1 <- Plot_Density_Custom(
    seurat_object = seurat_obj,
    reduction = reduction,
    pt.size = 0.1,
    features = c("infercnv_score")
  ) +
    ggtitle("Infercnv_score") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))

  Idents(seurat_obj) <- seurat_obj$ATAC.Subclass.June18_2025

  colors <- c("#a58874","#f4a85e","#e0a882","#f2c0d5","#ff9997")
  names(colors) <- levels(seurat_obj)

  p2 <- DimPlot(
    seurat_obj,
    reduction = reduction,
    label = TRUE,
    cols = colors
  ) +
    theme_void() +
    theme(legend.position = "none")

  p3 <- cowplot::plot_grid(p1, p2)

  # ---------------------------
  # Save plots
  # ---------------------------
  out_combined <- paste0(output_prefix, "_infercnv_score.pdf")

  
  ggsave(out_combined, plot = p3, width = 7, height = 3, dpi = 300)

  # ---------------------------
  # Export metadata
  # ---------------------------
  df <- seurat_obj@meta.data[, c("infercnv_score", "ATAC.Subclass.June18_2025"), drop = FALSE]

  out_csv <- paste0(output_prefix, "_infercnv_score.csv")
  write.csv(df, file = out_csv, row.names = TRUE)

  message("Saved Seurat object: ", out_rds)
  message("Saved combined figure: ", out_combined)
  message("Saved score: ", out_csv)
}

# ---------------------------
# CLI
# ---------------------------
main <- function() {
  parser <- ArgumentParser(description = "InferCNV scoring pipeline")

  parser$add_argument("-i", "--infercnv", required = TRUE, help = "infercnv RDS file")
  parser$add_argument("-s", "--seurat", required = TRUE, help = "Seurat RDS file")
  parser$add_argument("-o", "--output", default = "infercnv")

  parser$add_argument("--reduction", default = "Xumap_", help = "Reduction name (default: Xumap_)")

  args <- parser$parse_args()

  infercnv_score(
    infercnv_rds = args$infercnv,
    seurat_rds = args$seurat,
    output_prefix = args$output,
    reduction = args$reduction
  )
}

if (sys.nframe() == 0) {
  main()
}