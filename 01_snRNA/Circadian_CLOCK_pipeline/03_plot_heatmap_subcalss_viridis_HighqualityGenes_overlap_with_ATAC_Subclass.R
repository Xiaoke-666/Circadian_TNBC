library(pheatmap)
library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table) 
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)


celltypes<-c('Cancer_Epi',
 'CAFs',
 'CD4_T',
 'B',
 'Macro_Mono',
 'CD8_T',
 'NKT',
 'Neu'
)

phase_colors <- c(
  "Phase 1" = "#FF4545", "Phase 2" = "#FF9C73", "Phase 3" = "#FBD288",
  "Phase 4" = "#FCF596", "Phase 5" = "#DFF2EB", "Phase 6" = "#B9E5E8",
  "Phase 7" = "#7AB2D3", "Phase 8" = "#4A628A", "Unknown" = "gray"
)


setwd("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP/SingleCellVisualizationCircadianGenes_HighQualityGenes")
HighQualityGenes=read.csv("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/ATAC_RNA.intersected_genes.tsv",
                            sep="\t", header=TRUE)
indir="/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP/SingleCellVisualizationCircadianGenes/"



################################################
############ polished version:
celltypes<-c(
   'CAFs',
'Cancer_Epi',
 'CD4_T',
 'B',
 'Macro_Mono',
 'CD8_T',
 'NKT',
 'Neu'
)

for (celltype in celltypes) {
  message("Processing: ", celltype)

  # Load expression data
  expr_path <- file.path(indir, paste0("rna4.Cellbender.", celltype, ".CircadianGenes.ZT4_10_16_22.top.sorted.tsv.gz"))
  res <- fread(expr_path, header = TRUE, data.table = FALSE, sep = "\t")
  rownames(res) <- res$V1
  dim(res)

  res <- res[rownames(res) %in% HighQualityGenes$gene, -1]
  res <- na.omit(res)
  # print dim
  dim(res)

  # Smoothing and Z-score normalization
  pt.matrix <- apply(res, 1, function(x) smooth.spline(x, df = 3)$y)
  pt.matrix <- apply(pt.matrix, 1, function(x) (x - mean(x)) / sd(x))
  pt.matrix <- as.data.frame(pt.matrix)
  rownames(pt.matrix) <- rownames(res)


  message("Reading circadian gene annotation...")
  circadian_path <- paste0(
    "/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/",
    "01_snRNA_cellbender/01_map_snATAC_UMAP/CircadianGenes_v2/Integration/",
    celltype, ".discoODAres_CS_ANNOATION.outer.tsv"
  )
  CircadianGenes <- read.table(circadian_path, header = TRUE, sep = "\t", row.names = 1)
  CircadianGenes <- CircadianGenes[CircadianGenes$pvalue < 0.02, ]
  CircadianGenes <- CircadianGenes[order(CircadianGenes$acrophase), ]
  message("Number of raw circadian genes: ", nrow(CircadianGenes))

  CircadianGenes <- CircadianGenes[rownames(CircadianGenes) %in% rownames(pt.matrix), ]
  message("Number of circadian genes after filtering: ", nrow(CircadianGenes))

  # Assign circadian phases
  CircadianGenes$Phase <- with(CircadianGenes, case_when(
    acrophase >= 0 & acrophase < 3   ~ "Phase 1",
    acrophase >= 3 & acrophase < 6   ~ "Phase 2",
    acrophase >= 6 & acrophase < 9   ~ "Phase 3",
    acrophase >= 9 & acrophase < 12  ~ "Phase 4",
    acrophase >= 12 & acrophase < 15 ~ "Phase 5",
    acrophase >= 15 & acrophase < 18 ~ "Phase 6",
    acrophase >= 18 & acrophase < 21 ~ "Phase 7",
    acrophase >= 21 & acrophase <= 24 ~ "Phase 8",
    TRUE ~ "Unknown"
  ))

  # Filter to shared genes with annotations
  common_genes <- intersect(rownames(pt.matrix), rownames(CircadianGenes))
  pt.matrix <- pt.matrix[common_genes, ]
  row_split <- CircadianGenes[common_genes, "Phase", drop = TRUE]

  # Define annotation bar
  row_anno <- rowAnnotation(
    Phase = row_split,
    col = list(Phase = phase_colors),
    show_annotation_name = FALSE,
    width = unit(3, "mm")
  )

  # Heatmap color
  heatmap_color <- colorRamp2(seq(-1.5, 1.5, length.out = 11), viridis(11, option = "A"))
  nCells <- ncol(pt.matrix) / 4
  message("Number of cells: ", nCells)
  message("Plotting heatmap...")

  # Create heatmap object
  out <- Heatmap(
    pt.matrix,
    name = "z-score",
    col = heatmap_color,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 6),
    row_title_rot = 0,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    row_split = row_split,
    use_raster = TRUE,
    row_title = NULL
  )

  # File name prefix
  prefix <- paste0(celltype, ".CircadianGenes.heatmap.", nCells, "cells.smooth.spline")

  # Save PNGs
  png(paste0(prefix, ".png"), width = 4, height = 7, units = "in", res = 300)
  draw(row_anno + out, heatmap_legend_side = "right")
  dev.off()

  png(paste0(prefix, "2.png"), width = 4, height = 7, units = "in", res = 300)
  draw(row_anno + out, show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
  dev.off()

  # Save PDFs
  pdf(paste0(prefix, ".pdf"), width = 4, height = 7)
  draw(row_anno + out, heatmap_legend_side = "right")
  dev.off()

  pdf(paste0(prefix, "2.pdf"), width = 4, height = 7)
  draw(row_anno + out, show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
  dev.off()
}




'''
Processing: CAFs
Number of raw circadian genes: 721
Number of circadian genes after filtering: 687
Number of cells: 350

Processing: Cancer_Epi
Number of raw circadian genes: 2912
Number of circadian genes after filtering: 2755
Number of cells: 5769

Processing: CD4_T
Number of raw circadian genes: 2291
Number of circadian genes after filtering: 2221
Number of cells: 3363

Processing: B
Number of raw circadian genes: 1661
Number of circadian genes after filtering: 1624
Number of cells: 2183

Processing: Macro_Mono
Number of raw circadian genes: 1380
Number of circadian genes after filtering: 1331
Number of cells: 1006

Processing: CD8_T
Number of raw circadian genes: 786
Number of circadian genes after filtering: 763
Number of cells: 1125

Processing: NKT
Number of raw circadian genes: 391
Number of circadian genes after filtering: 376
Number of cells: 282

Processing: Neu
Number of raw circadian genes: 459
Number of circadian genes after filtering: 445
Number of cells: 395

'''