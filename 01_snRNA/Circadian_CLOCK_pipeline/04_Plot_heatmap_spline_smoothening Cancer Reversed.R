library(pheatmap)
library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table) 
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(circular)


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

matrixIndir="/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP/SingleCellVisualizationCircadianGenes/"

setwd("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP/SingleCellVisualizationCircadianGenes_Comparing_Subclass")

celltype1="Cancer_Epi"

celltypes<-c("CD4_T",
     'B',
      'CD8_T',
      'Macro_Mono',
      'NKT',
      'Neu'
)



##################################################################################################################################################################
######################################################  p1
##################################################################################################################################################################
setwd("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP/SingleCellVisualizationCircadianGenes_Comparing_Subclass")

for (celltype in celltypes){
  print(celltype)
  # generate a new folder
    dir.create(paste0(celltype1,"_vs_",celltype, ".Reversed.smooth.spline"), showWarnings = FALSE)
    setwd(paste0(celltype1,"_vs_",celltype, ".Reversed.smooth.spline"))
  
   print("Read circadian genes")
  CircadianGenes <- read.table(
    paste0("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP/subclass_acrophase_difference/",celltype1,"_with_",celltype, "_common_circadian_genes_withacrophasediff.csv"),
    header = TRUE, sep = ",", row.names = 1
  )
  CircadianGenes <- CircadianGenes[CircadianGenes$phase_category=="Reversed", ]

  CircadianGenes$acrophase <- CircadianGenes$ac1
  CircadianGenes <- CircadianGenes[order(CircadianGenes$acrophase), ]

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
 write.table(CircadianGenes, file = paste0(celltype1,"_with_",celltype, "_common_circadian_genes_withacrophasediff.csv"), sep = ",", row.names = TRUE, col.names = TRUE)
   print(head(CircadianGenes))

  res1 <- fread(paste0(matrixIndir,"rna4.Cellbender.",celltype1, ".CircadianGenes.ZT4_10_16_22.top.sorted.tsv.gz"),
               header = TRUE, data.table = FALSE, sep = "\t")
  rownames(res1) <- res1$V1
  res1 <- res1[, -1]  # Remove the first column (row names)
  res1 <- na.omit(res1)  # Remove rows with NA values
  # Smoothing and Z-score normalization
  pt.matrix1 <- apply(res1, 1, function(x) smooth.spline(x, df = 3)$y)
  pt.matrix1 <- apply(pt.matrix1, 1, function(x) (x - mean(x)) / sd(x))
  pt.matrix1 <- as.data.frame(pt.matrix1)
  rownames(pt.matrix1) <- rownames(res1)
  

 res2 <- fread(paste0(matrixIndir,"rna4.Cellbender.",celltype, ".CircadianGenes.ZT4_10_16_22.top.sorted.tsv.gz"),
               header = TRUE, data.table = FALSE, sep = "\t")
  rownames(res2) <- res2$V1
  res2 <- res2[, -1]  # Remove the first column (row names)
  res2 <- na.omit(res2)  # Remove rows with NA values
  # Smoothing and Z-score normalization
  pt.matrix2 <- apply(res2, 1, function(x) smooth.spline(x, df = 3)$y)
  pt.matrix2 <- apply(pt.matrix2, 1, function(x) (x - mean(x)) / sd(x))
  pt.matrix2 <- as.data.frame(pt.matrix2)
  rownames(pt.matrix2) <- rownames(res2)

  # Filter circadian genes to keep only those with Phase annotation
  common_genes <- intersect(rownames(pt.matrix1), rownames(CircadianGenes))
  common_genes <- intersect(rownames(pt.matrix2), common_genes)

  pt.matrix1 <- pt.matrix1[common_genes, ]
  pt.matrix2 <- pt.matrix2[common_genes, ]
  row_split <- CircadianGenes[common_genes, "Phase", drop = TRUE]


  # Create row annotation bar
  row_anno <- rowAnnotation(
    Phase = row_split,
    col = list(Phase = phase_colors),
    show_annotation_name = FALSE, 
    width = unit(3, "mm")
  )

  heatmap_color <- colorRamp2(seq(-1.5, 1.5, length.out = 11), viridis(11, option = "A"))

  nCells= ncol(pt.matrix1)/4
  print(paste0("Number of cells: ", nCells))
  print("Plotting heatmap for circadian genes...")
  # Generate heatmap
  out <- Heatmap(
    pt.matrix1,
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
  png(paste0(celltype1, ".CircadianGenes.heatmap.", nCells, "cells.smooth.spline.png"),
      width = 4, height = 7, units = "in", res = 300)
  draw(row_anno + out, heatmap_legend_side = "right")
  dev.off()

  # without legend 
  png(paste0(celltype1, ".CircadianGenes.heatmap.", nCells, "cells.smooth.spline2.png"),
      width = 4, height = 7, units = "in", res = 300)
  draw(row_anno + out, show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
  dev.off()
  
  setwd("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP/SingleCellVisualizationCircadianGenes_Comparing_Subclass")

}



##################################################################################################################################################################
######################################################  p2
##################################################################################################################################################################

setwd("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP/SingleCellVisualizationCircadianGenes_Comparing_Subclass")

for (celltype in celltypes){
  print(celltype)
  # generate a new folder
    dir.create(paste0(celltype1,"_vs_",celltype, ".Reversed.smooth.spline"), showWarnings = FALSE)
    setwd(paste0(celltype1,"_vs_",celltype, ".Reversed.smooth.spline"))
  
   print("Read circadian genes")
  CircadianGenes <- read.table(
    paste0("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP/subclass_acrophase_difference/",celltype1,"_with_",celltype, "_common_circadian_genes_withacrophasediff.csv"),
    header = TRUE, sep = ",", row.names = 1
  )
  
  CircadianGenes <- CircadianGenes[CircadianGenes$phase_category=="Reversed", ]

  CircadianGenes$acrophase <- CircadianGenes$ac1
  CircadianGenes <- CircadianGenes[order(CircadianGenes$acrophase), ]


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
 write.table(CircadianGenes, file = paste0(celltype1,"_with_",celltype, "_common_circadian_genes_withacrophasediff.csv"), sep = ",", row.names = TRUE, col.names = TRUE)
   print(head(CircadianGenes))

  res1 <- fread(paste0(matrixIndir,"rna4.Cellbender.",celltype1, ".CircadianGenes.ZT4_10_16_22.top.sorted.tsv.gz"),
               header = TRUE, data.table = FALSE, sep = "\t")
  rownames(res1) <- res1$V1
  res1 <- res1[, -1]  # Remove the first column (row names)
  res1 <- na.omit(res1)  # Remove rows with NA values
  # Smoothing and Z-score normalization
  pt.matrix1 <- apply(res1, 1, function(x) smooth.spline(x, df = 3)$y)
  pt.matrix1 <- apply(pt.matrix1, 1, function(x) (x - mean(x)) / sd(x))
  pt.matrix1 <- as.data.frame(pt.matrix1)
  rownames(pt.matrix1) <- rownames(res1)
  

 res2 <- fread(paste0(matrixIndir,"rna4.Cellbender.",celltype, ".CircadianGenes.ZT4_10_16_22.top.sorted.tsv.gz"),
               header = TRUE, data.table = FALSE, sep = "\t")
  rownames(res2) <- res2$V1
  res2 <- res2[, -1]  # Remove the first column (row names)
  res2 <- na.omit(res2)  # Remove rows with NA values
  # Smoothing and Z-score normalization
  pt.matrix2 <- apply(res2, 1, function(x) smooth.spline(x, df = 3)$y)
  pt.matrix2 <- apply(pt.matrix2, 1, function(x) (x - mean(x)) / sd(x))
  pt.matrix2 <- as.data.frame(pt.matrix2)
  rownames(pt.matrix2) <- rownames(res2)

  # Filter circadian genes to keep only those with Phase annotation
  common_genes <- intersect(rownames(pt.matrix1), rownames(CircadianGenes))
  common_genes <- intersect(rownames(pt.matrix2), common_genes)

  pt.matrix2 <- pt.matrix2[common_genes, ]
  row_split <- CircadianGenes[common_genes, "Phase", drop = TRUE]


  # Create row annotation bar
  row_anno <- rowAnnotation(
    Phase = row_split,
    col = list(Phase = phase_colors),
    show_annotation_name = FALSE, 
    width = unit(3, "mm")
  )

  heatmap_color <- colorRamp2(seq(-1.5, 1.5, length.out = 11), viridis(11, option = "A"))

  nCells= ncol(pt.matrix2)/4
  print(paste0("Number of cells: ", nCells))
  print("Plotting heatmap for circadian genes...")
  # Generate heatmap
  out <- Heatmap(
    pt.matrix2,
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

  setwd("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP/SingleCellVisualizationCircadianGenes_Comparing_Subclass")

}






  #save as pdf
  pdf(paste0(celltype1, ".CircadianGenes.heatmap.", nCells, "cells.smooth.spline.pdf"),
      width = 4, height = 7)
  draw(row_anno + out, heatmap_legend_side = "right")
  dev.off()

    # without legend 
  pdf(paste0(celltype1, ".CircadianGenes.heatmap.", nCells, "cells.smooth.spline2.pdf"),
      width = 4, height = 7)
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