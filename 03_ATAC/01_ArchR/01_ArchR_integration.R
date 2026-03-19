#!/usr/bin/env Rscript

# =============================================================================
# ArchR preprocessing and integration pipeline for snATAC-seq
# =============================================================================

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

# -----------------------------------------------------------------------------
# Global configuration
# -----------------------------------------------------------------------------
addArchRThreads(threads = 48)
addArchRGenome("mm10")

# Project-relative paths
project_dir <- "."
raw_dir <- file.path(project_dir, "02_snATAC", "01_raw")
archr_dir <- file.path(project_dir, "02_snATAC", "02_ArchR")
h5ad_dir <- file.path(
  project_dir,
  "02_snATAC", "02_h5ad_basedOnClusterPeaks", "02_overall_cellbender"
)

arrow_project_dir <- file.path(archr_dir, "01_TSS2")
filtered_project_dir <- file.path(archr_dir, "01_TSS2_filterDoublets")
annotated_project_dir <- file.path(archr_dir, "03_proj3")

dir.create(archr_dir, recursive = TRUE, showWarnings = FALSE)

# Optional: set MACS2 path externally if needed
pathToMacs2 <- "macs2"

# -----------------------------------------------------------------------------
# Step 1. Create Arrow files
# -----------------------------------------------------------------------------
message("Listing input fragment files...")
inputFiles <- list.files(raw_dir, pattern = ".gz$", full.names = TRUE)

sample_names <- basename(inputFiles) %>%
  str_remove(".atac_fragments.tsv.gz")

names(inputFiles) <- sample_names

message("Creating Arrow files...")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 2,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  QCDir = file.path(raw_dir, "QualityControl"),
  threads = getArchRThreads()
)

ArrowFiles <- list.files(raw_dir, pattern = ".arrow$", full.names = TRUE)

# -----------------------------------------------------------------------------
# Step 2. Infer doublets
# -----------------------------------------------------------------------------
message("Computing doublet scores...")
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

# -----------------------------------------------------------------------------
# Step 3. Create ArchR project
# -----------------------------------------------------------------------------
message("Creating ArchR project...")
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = arrow_project_dir,
  copyArrows = TRUE
)

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = arrow_project_dir,
  load = TRUE
)

write.table(
  proj@cellColData,
  file.path(arrow_project_dir, "CellColData.TSS2.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

# -----------------------------------------------------------------------------
# Step 4. Filter doublets
# -----------------------------------------------------------------------------
message("Filtering doublets...")
proj2 <- filterDoublets(proj)

saveArchRProject(
  ArchRProj = proj2,
  outputDirectory = filtered_project_dir,
  load = TRUE
)

write.table(
  proj2@cellColData,
  file.path(filtered_project_dir, "CellColData.TSS2_filterDoublets.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

# -----------------------------------------------------------------------------
# Step 5. Dimensionality reduction
# -----------------------------------------------------------------------------
message("Loading filtered ArchR project...")
proj2 <- loadArchRProject(filtered_project_dir)

message("Running iterative LSI...")
proj2 <- addIterativeLSI(
  ArchRProj = proj2,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 20000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  force = TRUE
)

# -----------------------------------------------------------------------------
# Step 6. Clustering and UMAP
# -----------------------------------------------------------------------------
message("Adding clusters...")
proj2 <- addClusters(
  input = proj2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

message("Computing UMAP...")
proj2 <- addUMAP(
  ArchRProj = proj2,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine"
)

# -----------------------------------------------------------------------------
# Step 7. Visualization
# -----------------------------------------------------------------------------
message("Plotting UMAPs...")
p1 <- plotEmbedding(
  ArchRProj = proj2,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP"
)

p2 <- plotEmbedding(
  ArchRProj = proj2,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP"
)

ggAlignPlots(p1, p2, type = "h")

plotPDF(
  p1, p2,
  name = "Plot-UMAP-Sample-Clusters.pdf",
  ArchRProj = proj2,
  addDOC = FALSE,
  width = 5,
  height = 5
)

saveArchRProject(
  ArchRProj = proj2,
  outputDirectory = filtered_project_dir,
  load = TRUE
)

# -----------------------------------------------------------------------------
# Step 8. Peak calling
# -----------------------------------------------------------------------------
message("Adding group coverages...")
# NOTE:
# Peak calling was performed in two stages.
# First, we used default ArchR parameters (addGroupCoverages and addReproduciblePeakSet)
# to generate an initial peak set for each cluster. This peak set was then used as the
# genomic feature space for downstream analysis in Muon, where cell clustering and
# annotation were refined using Muon.
#
# After finalizing cell type annotations, we repeated the peak calling procedure
# using the updated cluster labels to obtain the final, high-confidence peak set.

proj2 <- addGroupCoverages(
  ArchRProj = proj2,
  groupBy = "Clusters",
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
proj2 <- addReproduciblePeakSet(
  ArchRProj = proj2,
  maxPeaks = 500000,
  cutOff = 0.01,
  force = TRUE,
  groupBy = "Clusters",
  pathToMacs2 = pathToMacs2
)

message("Adding peak matrix...")
proj2 <- addPeakMatrix(proj2)

peakSet_df <- data.frame(proj2@peakSet)
peakSet_df$Cell_type <- names(proj2@peakSet)

write.table(
  peakSet_df,
  file.path(filtered_project_dir, "proj2_cluster.peak.mtx.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

saveArchRProject(
  ArchRProj = proj2,
  outputDirectory = filtered_project_dir,
  load = TRUE
)

# -----------------------------------------------------------------------------
# Step 9. Add imputation weights
# -----------------------------------------------------------------------------
message("Adding imputation weights...")
proj2 <- addImputeWeights(proj2)

# -----------------------------------------------------------------------------
# Step 10. Integrate external UMAP and metadata
# -----------------------------------------------------------------------------
message("Reading external UMAP coordinates...")
X_umap <- read.csv(
  file.path(h5ad_dir, "snATAC.X_umap.tsv.gz"), # This is the umap coordinates from the h5ad file computed from muon
  sep = "\t"
)

X_umap_df <- X_umap[, -1, drop = FALSE]
rownames(X_umap_df) <- X_umap[, 1]

# Subset project to cells present in external UMAP
common_cells <- rownames(proj2@cellColData) %in% rownames(X_umap_df)
proj3 <- proj2[common_cells, ]

X_umap_matched <- X_umap_df[rownames(proj3@cellColData), , drop = FALSE]

proj3@embeddings$X_umap$df <- X_umap_matched
colnames(proj3@embeddings$X_umap$df) <- c(
  "IterativeLSI#UMAP_Dimension_1",
  "IterativeLSI#UMAP_Dimension_2"
)

# -----------------------------------------------------------------------------
# Step 11. Add cell-level annotations from external table
# -----------------------------------------------------------------------------
message("Reading cell metadata...")
cells <- read.csv(
  file.path(h5ad_dir, "snATAC.obs.tsv.gz"), # This is the h5ad.obs  from the h5ad file computed from muon
  sep = "\t"
)

cell_index <- match(rownames(proj3@cellColData), cells$X)

proj3@cellColData$batch <- cells[cell_index, "batch"]
proj3@cellColData$TimePoint <- cells[cell_index, "TimePoint"]
proj3@cellColData$CD45 <- cells[cell_index, "CD45"]
proj3@cellColData$Replicate <- cells[cell_index, "Replicate"]

proj3@cellColData$ATAC_CellTypes <- cells[cell_index, "ATAC_CellTypes"]
proj3@cellColData$ATAC_CellSubTypes <- cells[cell_index, "ATAC_CellSubTypes"]

# -----------------------------------------------------------------------------
# Step 12. Save annotated project
# -----------------------------------------------------------------------------
message("Saving final annotated ArchR project...")
saveArchRProject(
  ArchRProj = proj3,
  outputDirectory = annotated_project_dir,
  load = TRUE
)

message("Done.")