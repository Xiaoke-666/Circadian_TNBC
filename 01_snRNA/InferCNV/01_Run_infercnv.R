# Install infercnv
# devtools::install_github("broadinstitute/infercnv")
setwd("~/cirdian")
# Run Infercnv
library(Seurat)
library(infercnv)
seurat_obj <- readRDS("./raw_data.rds")
# Subset for CD45- proportion
seurat_obj <- subset(seurat_obj, subset = CD45 == "CD45-")
# Extract raw counts
raw_counts = seurat_obj[["RNA"]]$counts
Idents(seurat_obj) <- seurat_obj$ATAC.Subclass.June18_2025
anno <- as.matrix(seurat_obj@active.ident)
# Mouse gene code
gene_order_file <- "./mouse_gencode.txt"
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = raw_counts,
                                     annotations_file = anno,
                                     gene_order_file = gene_order_file,
                                     delim = "\t",
                                     ref_group_names = unique(seurat_obj$ATAC.Subclass.June18_2025)[-1])
infercnv_obj_run <- run(infercnv_obj = infercnv_obj,
                        out_dir = "infercnv_output",
                        cutoff = 0.1,
                        output_format = "pdf",
                        denoise = T,
                        cluster_by_groups = T)