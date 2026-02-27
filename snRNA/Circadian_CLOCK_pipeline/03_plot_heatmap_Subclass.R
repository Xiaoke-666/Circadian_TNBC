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


setwd("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP/SingleCellVisualizationCircadianGenes")

for (celltype in celltypes){
  print(celltype)
  res <- fread(paste0("rna4.Cellbender.",celltype, ".CircadianGenes.ZT4_10_16_22.top.sorted.tsv.gz"),
               header = TRUE, data.table = FALSE, sep = "\t")
  rownames(res) <- res$V1
  res <- res[, -1]  # Remove the first column (row names)
  res <- na.omit(res)  # Remove rows with NA values
  # Smoothing and Z-score normalization
  pt.matrix <- apply(res, 1, function(x) smooth.spline(x, df = 3)$y)
  pt.matrix <- apply(pt.matrix, 1, function(x) (x - mean(x)) / sd(x))
  pt.matrix <- as.data.frame(pt.matrix)
  rownames(pt.matrix) <- rownames(res)
  
  print("Read circadian genes")
  CircadianGenes <- read.table(
    paste0("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP/CircadianGenes_v2/Integration/",
           celltype, ".discoODAres_CS_ANNOATION.outer.tsv"),
    header = TRUE, sep = "\t", row.names = 1
  )
  CircadianGenes <- CircadianGenes[CircadianGenes$pvalue < 0.02, ]
  CircadianGenes <- CircadianGenes[order(CircadianGenes$acrophase), ]
  # Assign circadian phases (P1–P8)
 
  CircadianGenes$Phase <- case_when(
    CircadianGenes$acrophase >= 0 & CircadianGenes$acrophase < 3  ~ "Phase 1",
    CircadianGenes$acrophase >= 3 & CircadianGenes$acrophase < 6  ~ "Phase 2",
    CircadianGenes$acrophase >= 6 & CircadianGenes$acrophase < 9  ~ "Phase 3",
    CircadianGenes$acrophase >= 9 & CircadianGenes$acrophase < 12 ~ "Phase 4",
    CircadianGenes$acrophase >= 12 & CircadianGenes$acrophase < 15 ~ "Phase 5",
    CircadianGenes$acrophase >= 15 & CircadianGenes$acrophase < 18 ~ "Phase 6",
    CircadianGenes$acrophase >= 18 & CircadianGenes$acrophase < 21 ~ "Phase 7",
    CircadianGenes$acrophase >= 21 & CircadianGenes$acrophase <= 24 ~ "Phase 8",
    TRUE ~ "Unknown"
  )

  # Filter circadian genes to keep only those with Phase annotation
  common_genes <- intersect(rownames(pt.matrix), rownames(CircadianGenes))
  pt.matrix <- pt.matrix[common_genes, ]
  row_split <- CircadianGenes[common_genes, "Phase", drop = TRUE]


  # Create row annotation bar
  row_anno <- rowAnnotation(
    Phase = row_split,
    col = list(Phase = phase_colors),
    show_annotation_name = FALSE, 
    width = unit(3, "mm")
  )


  heatmap_color <- colorRamp2(seq(-1.5, 1.5, length.out = 11), viridis(11, option = "A"))
  nCells= ncol(pt.matrix)/4

  print(paste0("Number of cells: ", nCells))
  print("Plotting heatmap for circadian genes...")
  # Generate heatmap
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

  # Save to PNG
  png(paste0(celltype, ".CircadianGenes.heatmap.", nCells, "cells.smooth.spline.png"),
      width = 4, height = 7, units = "in", res = 300)
  draw(row_anno + out, heatmap_legend_side = "right")
  dev.off()

  # without legend 
  png(paste0(celltype, ".CircadianGenes.heatmap.", nCells, "cells.smooth.spline2.png"),
      width = 4, height = 7, units = "in", res = 300)
  draw(row_anno + out, show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
  dev.off()

  #save as pdf
  pdf(paste0(celltype, ".CircadianGenes.heatmap.", nCells, "cells.smooth.spline.pdf"),
      width = 4, height = 7)
  draw(row_anno + out, heatmap_legend_side = "right")
  dev.off()

    # without legend 
  pdf(paste0(celltype, ".CircadianGenes.heatmap.", nCells, "cells.smooth.spline2.pdf"),
      width = 4, height = 7)
  draw(row_anno + out, show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
  dev.off()


}

