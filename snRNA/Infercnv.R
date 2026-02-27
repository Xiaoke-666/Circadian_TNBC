library(Seurat)
library(schard)
# Install infercnv
# devtools::install_github("broadinstitute/infercnv")
# Convert to rds
setwd("~/cirdian")
seurat_obj <- h5ad2seurat("raw_data.h5ad")
saveRDS(seurat_obj, file = "raw_data.rds")
seurat_obj <- readRDS("./raw_data.rds")
# Run Infercnv
library(Seurat)
library(infercnv)
seurat_obj <- readRDS("./raw_data.rds")
seurat_obj <- subset(seurat_obj, subset = CD45 == "CD45-")
raw_counts = seurat_obj[["RNA"]]$counts
Idents(seurat_obj) <- seurat_obj$ATAC.Subclass.June18_2025
anno <- as.matrix(seurat_obj@active.ident)
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
# read in data
library(tidyr)
library(dplyr)
library(tibble)
setwd("~/cirdian/infercnv_output")
infercnv_obj <- readRDS("run.final.infercnv_obj")
# Rescale from -1 to 1
exp <- infercnv_obj@expr.data
exp <- scale(t(exp))
tmp1 <- sweep(exp,2,apply(exp,2,min),"-")
tmp2 <- apply(exp,2,max)-apply(exp,2,min)
exp_1 <- t(2*sweep(tmp1,2,tmp2,"/")-1)
cnv_score <- as.data.frame(colSums(exp_1*exp_1))
colnames(cnv_score) <- "score"
# Add to seurat
seurat_obj <- readRDS("../raw_data.rds")
Idents(seurat_obj) <- seurat_obj$ATAC.Subclass.June18_2025
seurat_obj <- subset(seurat_obj, subset = CD45 == "CD45-")
metadata <- seurat_obj@meta.data
cnv_score <- cnv_score[match(rownames(metadata),rownames(cnv_score)),,drop = F]
metadata$infercnv_score <- cnv_score$score
seurat_obj@meta.data <- metadata
saveRDS(seurat_obj, file = "../CD45-subset.rds")
# Plot
library(Seurat)
library(scCustomize)
library(RColorBrewer)
library(ggplot2)
setwd("~/cirdian")
seurat_obj <- readRDS("./CD45-subset.rds")
p1 <- Plot_Density_Custom(seurat_object = seurat_obj, reduction = "Xumap_",
                          pt.size = 0.1,
                          features = c("infercnv_score"))+
  ggtitle("Infercnv_score")+
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5))
Idents(seurat_obj) <- seurat_obj$ATAC.Subclass.June18_2025
colors <- c("#a58874","#f4a85e","#e0a882","#f2c0d5","#ff9997")
names(colors) <- levels(seurat_obj)
p2 <- DimPlot(seurat_obj, reduction = "Xumap_", label = T, cols = colors)+
  theme_void()+
  theme(legend.position = "none")
p3 <- cowplot::plot_grid(p1,p2)
ggsave(p2, filename = "./umap.pdf", width = 4, height = 3, dpi = 300)
ggsave(p3, filename = "./figures/infercnv_score.pdf", width = 7, height = 3, dpi = 300)

df <- seurat_obj@meta.data[,c("infercnv_score","ATAC.Subclass.June18_2025"),drop = F]
write.csv(df, file = "./figures/infercnv_score.csv")