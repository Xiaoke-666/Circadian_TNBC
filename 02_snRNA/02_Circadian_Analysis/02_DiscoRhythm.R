#!/usr/bin/env Rscript

# =============================================================================
# DiscoRhythm pipeline for pseudobulk RNA-seq circadian analysis
# =============================================================================

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(DiscoRhythm))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------
parser <- ArgumentParser(description = "Run DiscoRhythm on pseudobulk data")

parser$add_argument("-i", "--input", required=TRUE, help="Input directory")
parser$add_argument("-d", "--data", required=TRUE, help="Expression matrix (TSV)")
parser$add_argument("-m", "--metadata", required=TRUE, help="Metadata (TSV)")
parser$add_argument("-p", "--prefix", required=TRUE, help="Output prefix")
parser$add_argument("-o", "--output", required=TRUE, help="Output directory")

args <- parser$parse_args()

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
setwd(args$input)

cat("Reading input data...\n")

expr_df <- read.csv(args$data, sep = "\t")
metadata <- read.csv(args$metadata, sep = "\t")

prefix <- args$prefix

setwd(args$output)

# -----------------------------------------------------------------------------
# Metadata processing
# -----------------------------------------------------------------------------
metadata$TimePoint <- factor(
  metadata$TimePoint,
  levels = c("ZT4","ZT10","ZT16","ZT22")
)

# Order by time then replicate
metadata <- metadata[order(metadata$TimePoint, metadata$Replicate),]

# Create DiscoRhythm sample names
metadata$tmp_Name <- paste0("Sample", seq_len(nrow(metadata)))
metadata$DiscoRhythm_Name <- paste0(
  metadata$TimePoint, "_",
  metadata$tmp_Name, "_rep", metadata$Replicate
)

# -----------------------------------------------------------------------------
# Expression matrix preparation
# -----------------------------------------------------------------------------
cat("Preparing expression matrix...\n")

indata <- expr_df[, c("GeneID", metadata$SampleName)]
colnames(indata) <- c("IDs", metadata$DiscoRhythm_Name)

# Filter low-expression genes
indata <- indata[rowSums(indata[,-1]) > 0.2, ]
indata <- na.omit(indata)

# -----------------------------------------------------------------------------
# Metadata for DiscoRhythm
# -----------------------------------------------------------------------------
tmp_metadata <- metadata[, c("DiscoRhythm_Name","TimePoint","Replicate")]
rownames(tmp_metadata) <- tmp_metadata$DiscoRhythm_Name

tmp_metadata$TimePoint <- as.numeric(str_replace(tmp_metadata$TimePoint, "ZT", ""))
tmp_metadata$Replicate <- paste0("rep", tmp_metadata$Replicate)

colnames(tmp_metadata) <- c("ID","Time","ReplicateID")

# -----------------------------------------------------------------------------
# DiscoRhythm pipeline
# -----------------------------------------------------------------------------
cat("Running DiscoRhythm pipeline...\n")

# 1. Convert to SE
se <- discoDFtoSE(indata)

# 2. Input check
selectDataSE <- discoCheckInput(se)

# 3–5. Outlier detection
CorRes <- discoInterCorOutliers(
  selectDataSE,
  cor_method="pearson",
  threshold=3,
  thresh_type="sd"
)

PCAres <- discoPCAoutliers(
  selectDataSE,
  threshold=3,
  scale=TRUE,
  pcToCut = c("PC1","PC2","PC3","PC4")
)

# Remove outliers
FilteredSE <- selectDataSE[, !PCAres$outliers & !CorRes$outliers]

# Save outliers
outliers <- metadata[
  metadata$DiscoRhythm_Name %in%
    colnames(selectDataSE)[PCAres$outliers | CorRes$outliers],
]

write.table(
  outliers,
  file = paste0(prefix, ".outliers.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Replicate consistency filtering
# -----------------------------------------------------------------------------
ANOVAres <- discoRepAnalysis(
  FilteredSE,
  aov_method="Equal Variance",
  aov_pcut=0.05,
  aov_Fcut=1,
  avg_method="Median"
)

FinalSE <- ANOVAres$se

# -----------------------------------------------------------------------------
# Period detection
# -----------------------------------------------------------------------------
PeriodRes <- discoPeriodDetection(
  FinalSE,
  timeType="linear",
  main_per=24
)

# -----------------------------------------------------------------------------
# Oscillation detection (core result)
# -----------------------------------------------------------------------------
FinalSE2 <- FinalSE[rowSums(assay(FinalSE)) > 0, ]

discoODAres <- discoODAs(
  FinalSE2,
  period=24,
  method="CS",
  ncores=1,
  circular_t=FALSE
)

results <- discoODAres$CS
results$GeneID <- rownames(results)

# -----------------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------------
cat("Saving results...\n")

write.table(
  results,
  paste0(prefix, ".discoODAres.tsv"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# Merge with expression matrix
results_sorted <- results[order(results$acrophase), ]
indata2 <- indata[indata$IDs %in% results_sorted$GeneID, ]
colnames(indata2)[1] <- "GeneID"

results_merged <- merge(results_sorted, indata2, by="GeneID")

write.table(
  results_merged,
  paste0(prefix, ".discoODAres_with_expression.tsv"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

cat("Done.\n")