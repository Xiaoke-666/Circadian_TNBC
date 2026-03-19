#!/usr/bin/env Rscript

# =============================================================================
# Downstream ArchR analysis for snATAC-seq
#
# This script performs:
# 1. project loading and metadata integration
# 2. UMAP visualization
# 3. marker gene analysis
# 4. peak calling by annotated cell type
# 5. peak export
# 6. marker peak and motif deviation analysis
# 7. RNA-ATAC integration
# 8. co-accessibility and peak-to-gene linkage
# 9. grouped matrix export
# 10. group-level bigWig generation
# =============================================================================

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
addArchRThreads(threads = 48)
addArchRGenome("mm10")

project_dir <- "."
workdir <- file.path(project_dir, "02_snATAC", "02_ArchR", "03_proj3")
rna_dir <- file.path(project_dir, "01_snRNA_cellbender", "01_map_snATAC_UMAP")
h5ad_dir <- file.path(project_dir, "02_snATAC", "02_h5ad_basedOnClusterPeaks", "02_overall_cellbender")

path_to_macs2 <- "macs2"

dir.create(workdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load ArchR project
# -----------------------------------------------------------------------------
message("Loading ArchR project...")
proj3 <- loadArchRProject(workdir)

saveArchRProject(
  ArchRProj = proj3,
  outputDirectory = workdir,
  load = TRUE
)

# -----------------------------------------------------------------------------
# UMAP visualization
# -----------------------------------------------------------------------------
message("Plotting UMAP by ATAC cell types...")
p1 <- plotEmbedding(
  ArchRProj = proj3,
  colorBy = "cellColData",
  name = "ATAC_CellTypes",
  embedding = "X_umap",
  size = 0.1,
  imputeWeights = getImputeWeights(proj3)
)

plotPDF(
  p1,
  name = "proj3.ATAC_CellTypes.pdf",
  ArchRProj = proj3,
  addDOC = FALSE,
  width = 5,
  height = 5
)

# -----------------------------------------------------------------------------
# Marker genes by annotated ATAC cell type
# -----------------------------------------------------------------------------
message("Computing imputation weights...")
proj3 <- addImputeWeights(proj3)

message("Computing marker genes...")
markersGS <- getMarkerFeatures(
  ArchRProj = proj3,
  useMatrix = "GeneScoreMatrix",
  groupBy = "ATAC_CellTypes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveRDS(markersGS, "proj3.ATAC_CellTypes.Marker_genes.rds")

markerList <- getMarkers(
  markersGS,
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5"
)

for (i in seq_along(markerList)) {
  write.csv(
    markerList[[i]],
    file = paste0(
      "proj3.CellType.Marker_genes.",
      names(markerList)[i],
      ".csv"
    ),
    quote = FALSE,
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# Peak calling by annotated ATAC cell type
# -----------------------------------------------------------------------------
message("Adding group coverages...")

# NOTE:
# Peak calling was performed after finalizing ATAC cell type annotations.
# Group-level coverages were generated for each annotated cell type, followed by
# reproducible peak calling to define a refined peak set for downstream analyses.

proj3 <- addGroupCoverages(
  ArchRProj = proj3,
  groupBy = "ATAC_CellTypes",
  useLabels = TRUE,
  minCells = 40,
  maxCells = 10000,
  maxFragments = 25 * 10^6,
  minReplicates = 2,
  maxReplicates = 30,
  sampleRatio = 0.8,
  kmerLength = 6,
  threads = getArchRThreads(),
  returnGroups = FALSE,
  parallelParam = NULL,
  force = TRUE,
  verbose = TRUE,
  logFile = createLogFile("addGroupCoverages")
)

message("Calling reproducible peak sets...")
proj3 <- addReproduciblePeakSet(
  ArchRProj = proj3,
  maxPeaks = 500000,
  cutOff = 0.01,
  force = TRUE,
  groupBy = "ATAC_CellTypes",
  pathToMacs2 = path_to_macs2
)

message("Adding peak matrix...")
proj3 <- addPeakMatrix(proj3)

peakSet_df <- data.frame(proj3@peakSet)
peakSet_df$Cell_type <- names(proj3@peakSet)

write.table(
  peakSet_df,
  "ATAC_CellTypes_peakSet.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

saveArchRProject(
  ArchRProj = proj3,
  outputDirectory = workdir,
  load = TRUE
)

# -----------------------------------------------------------------------------
# Export reproducible peaks as TSV and BED
# -----------------------------------------------------------------------------
peakcalls_dir <- file.path(workdir, "PeakCalls")
if (dir.exists(peakcalls_dir)) {
  message("Exporting reproducible peaks...")
  old_wd <- getwd()
  setwd(peakcalls_dir)

  peak_files <- list.files(pattern = "reproduciblePeaks.gr.rds$")
  all_peaks <- data.frame()

  for (f in peak_files) {
    tmp <- readRDS(f)
    tmp_df <- data.frame(tmp)
    cell_type <- str_replace(f, "-reproduciblePeaks.gr.rds", "")
    tmp_df$Cell_type <- cell_type

    all_peaks <- rbind(all_peaks, tmp_df)

    write.table(
      tmp_df,
      paste0(cell_type, ".tsv"),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE
    )

    write.table(
      tmp_df[, c(1, 2, 3)],
      paste0(cell_type, ".bed"),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
  }

  write.table(
    all_peaks,
    "CellType_reproduciblePeaks.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  system('for i in *.bed; do sort -k1,1 -k2,2n "$i" > "${i%.*}.srt.bed"; done')
  setwd(old_wd)
}

# -----------------------------------------------------------------------------
# Marker peaks
# -----------------------------------------------------------------------------
message("Computing marker peaks...")
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj3,
  useMatrix = "PeakMatrix",
  groupBy = "ATAC_CellTypes",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

marker_peak_list <- getMarkers(
  markersPeaks,
  cutOff = "FDR <= 0.01 & Log2FC >= 1"
)

# -----------------------------------------------------------------------------
# Motif annotations and chromVAR deviations
# -----------------------------------------------------------------------------
message("Adding motif annotations...")
proj3 <- addMotifAnnotations(
  ArchRProj = proj3,
  motifSet = "cisbp",
  name = "Motif"
)

proj3 <- addBgdPeaks(proj3)

message("Computing motif deviation matrix...")
proj3 <- addDeviationsMatrix(
  ArchRProj = proj3,
  peakAnnotation = "Motif",
  force = TRUE
)

message("Exporting group-level chromVAR z-scores...")
chromvar_z <- getGroupSE(
  ArchRProj = proj3,
  useMatrix = "MotifMatrix",
  groupBy = "ATAC_CellTypes",
  divideN = TRUE,
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)

saveRDS(chromvar_z, "ChromVar_Zscore_Mat.CellType.rds")

chromvar_z_filter <- chromvar_z[rowData(chromvar_z)$seqnames == "z", ]

motif_mat <- as.data.frame(assay(chromvar_z_filter, "MotifMatrix"))
motif_mat$ID <- rownames(motif_mat)

motif_anno <- as.data.frame(rowData(chromvar_z_filter))
motif_anno$ID <- rownames(motif_anno)

motif_export <- merge(motif_anno, motif_mat, by = "ID")

write.table(
  motif_export,
  file = "ChromVar_Zscore_Mat.CellType.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Cross-platform RNA-ATAC integration
# -----------------------------------------------------------------------------
message("Loading matched RNA reference...")
seRNA <- readRDS(
  file.path(rna_dir, "rna3.Cellbender.ATAC_UMAP.Overlap.July2_2024.sce.rds")
)

celltypes_atac <- unique(proj3$ATAC_CellTypes)
celltypes_rna <- unique(colData(seRNA)$ATAC_CellTypes)
shared_celltypes <- intersect(celltypes_atac, celltypes_rna)

groupList <- list()

for (x in shared_celltypes) {
  tmp <- SimpleList(
    ATAC = proj3$cellNames[proj3$ATAC_CellTypes %in% x],
    RNA = colnames(seRNA)[grep(x, colData(seRNA)$ATAC_CellTypes)]
  )
  message(x)
  message(length(tmp$ATAC))
  message(length(tmp$RNA))
  groupList[[x]] <- tmp
}

RNA_embedding <- as.data.frame(reducedDims(seRNA)[["UMAP"]])

message("Running gene integration...")
proj3 <- addGeneIntegrationMatrix(
  ArchRProj = proj3,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE,
  groupList = groupList,
  embeddingRNA = RNA_embedding,
  groupATAC = "ATAC_CellTypes",
  groupRNA = "ATAC_CellTypes",
  nameCell = "predictedCell_Co",
  nameGroup = "predictedGroup_Co",
  nameScore = "predictedScore_Co",
  force = TRUE
)

pal <- paletteDiscrete(values = colData(seRNA)$ATAC_CellTypes)

p_integration <- plotEmbedding(
  proj3,
  colorBy = "cellColData",
  name = "predictedGroup_Co",
  pal = pal,
  embedding = "X_umap"
)

plotPDF(
  p_integration,
  name = "Plot-UMAP-RNA-Integration.pdf",
  ArchRProj = proj3,
  addDOC = FALSE,
  width = 5,
  height = 5
)

# -----------------------------------------------------------------------------
# Co-accessibility
# -----------------------------------------------------------------------------
message("Computing co-accessibility...")
proj3 <- addCoAccessibility(
  ArchRProj = proj3,
  reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
  ArchRProj = proj3,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)

saveRDS(cA, "CoAccessibility.ATAC.CellType.rds")

cA_df <- as.data.frame(cA)
gr_df <- data.frame(metadata(cA)[[1]])

write.table(
  cA_df,
  "CoAccessibility.peak.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

write.table(
  data.frame(idx = rownames(gr_df), gr_df),
  "CoAccessibility.peak.idx.txt",
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

# -----------------------------------------------------------------------------
# Peak-to-gene links
# -----------------------------------------------------------------------------
message("Computing peak-to-gene links...")
proj3 <- addPeak2GeneLinks(
  ArchRProj = proj3,
  reducedDims = "IterativeLSI"
)

export_p2g <- function(p2g, prefix) {
  p2g_df <- as.data.frame(p2g)
  metadata_peak <- as.data.frame(metadata(p2g)[[1]])
  metadata_gene <- as.data.frame(metadata(p2g)[[2]])

  p2g_df$peak_chr <- metadata_peak[p2g_df$idxATAC, "seqnames"]
  p2g_df$peak_start <- metadata_peak[p2g_df$idxATAC, "start"]
  p2g_df$peak_end <- metadata_peak[p2g_df$idxATAC, "end"]

  p2g_df$gene_chr <- metadata_gene[p2g_df$idxRNA, "seqnames"]
  p2g_df$TSS <- metadata_gene[p2g_df$idxRNA, "start"]
  p2g_df$gene_name <- metadata_gene[p2g_df$idxRNA, "name"]

  write.table(
    p2g_df,
    paste0(prefix, ".tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  p2g_bed <- p2g_df[, c(
    "peak_chr", "peak_start", "peak_end", "Correlation",
    "FDR", "VarQATAC", "VarQRNA", "gene_chr", "TSS", "gene_name"
  )]

  write.table(
    p2g_bed,
    paste0(prefix, ".bed"),
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

p2g_res1 <- getPeak2GeneLinks(
  ArchRProj = proj3,
  corCutOff = 0.3,
  resolution = 1,
  FDRCutOff = 1e-04,
  returnLoops = FALSE
)

saveRDS(p2g_res1, "cCRE_gene_pairs.res1.rds")
export_p2g(p2g_res1, "cCRE_gene_pairs.res1.cor03")

p2g_res1000 <- getPeak2GeneLinks(
  ArchRProj = proj3,
  corCutOff = 0.3,
  resolution = 1000,
  returnLoops = FALSE
)

saveRDS(p2g_res1000, "cCRE_gene_pairs.res1000.rds")
export_p2g(p2g_res1000, "cCRE_gene_pairs.res1000.cor03")

saveArchRProject(
  ArchRProj = proj3,
  outputDirectory = workdir,
  load = TRUE
)

# -----------------------------------------------------------------------------
# Export grouped gene score and peak matrices
# -----------------------------------------------------------------------------
message("Exporting grouped gene score matrix by ATAC cell type...")
geneMat <- getGroupSE(
  ArchRProj = proj3,
  useMatrix = "GeneScoreMatrix",
  groupBy = "ATAC_CellTypes",
  divideN = TRUE,
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)

gene_mat <- assay(geneMat, "GeneScoreMatrix")
rownames(gene_mat) <- getFeatures(proj3)

write.table(
  data.frame(gene = rownames(gene_mat), gene_mat),
  file = "GeneMat.CellType.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE
)

proj3@cellColData$CellTypes_TimePoint <- paste0(
  proj3@cellColData$ATAC_CellTypes, ".", proj3@cellColData$TimePoint
)

message("Exporting grouped gene score matrix by cell type and time point...")
geneMat_tp <- getGroupSE(
  ArchRProj = proj3,
  useMatrix = "GeneScoreMatrix",
  groupBy = "CellTypes_TimePoint",
  divideN = TRUE,
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)

gene_mat_tp <- assay(geneMat_tp, "GeneScoreMatrix")
rownames(gene_mat_tp) <- getFeatures(proj3)

write.table(
  data.frame(gene = rownames(gene_mat_tp), gene_mat_tp),
  file = "GeneMat.CellTypes_TimePoint.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE
)

message("Exporting grouped peak matrix by ATAC cell type...")
PeakMatrix <- getGroupSE(
  ArchRProj = proj3,
  useMatrix = "PeakMatrix",
  groupBy = "ATAC_CellTypes",
  divideN = TRUE,
  scaleTo = 10000000,
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)

peak_mat <- assays(PeakMatrix)$PeakMatrix
peak_anno <- rowData(PeakMatrix)
rownames(peak_mat) <- paste0(peak_anno$seqnames, ":", peak_anno$start, "-", peak_anno$end)

write.table(
  peak_mat,
  "PeakScoreMatrix.CellType.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE
)

message("Exporting grouped peak matrix by cell type and time point...")
PeakMatrix_tp <- getGroupSE(
  ArchRProj = proj3,
  useMatrix = "PeakMatrix",
  groupBy = "CellTypes_TimePoint",
  divideN = TRUE,
  scaleTo = 10000000,
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)

peak_mat_tp <- assays(PeakMatrix_tp)$PeakMatrix
peak_anno_tp <- rowData(PeakMatrix_tp)
rownames(peak_mat_tp) <- paste0(peak_anno_tp$seqnames, ":", peak_anno_tp$start, "-", peak_anno_tp$end)

write.table(
  peak_mat_tp,
  "PeakScoreMatrix.CellTypes_TimePoint.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE
)

# -----------------------------------------------------------------------------
# Update annotations and export group bigWig files
# -----------------------------------------------------------------------------
message("Reading updated annotations for bigWig export...")
cell <- read.csv(
  file.path(rna_dir, "rna5.Cellbender.ATAC_UMAP.Overlap.June18_2025.obs.tsv.gz"),
  sep = "\t"
)

cell_index <- match(rownames(proj3@cellColData), cell$X)

proj3@cellColData$CellType.Time_Points <- paste0(
  proj3@cellColData$ATAC_CellType, ".", proj3@cellColData$TimePoint
)

message("Exporting group bigWig files...")
getGroupBW(
  ArchRProj = proj3,
  groupBy = "CellTypes_TimePoint",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 50000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)

message("Done.")