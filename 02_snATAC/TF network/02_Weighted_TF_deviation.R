library(Seurat)
library(schard)
library(tidyverse)

#### For cell type
setwd("~/cirdian/TF_connectivity/")
tf_score <- h5ad2seurat("./ChromVar_Deviations.Clean.July17_2025.h5ad")
id <- rownames(tf_score@meta.data)[-which(tf_score$ATAC.Subclass.June18_2025 %in% c("Normal_Epi","Endo","PVLs","Plasmablasts","DC","NK",NA))]
sub <- subset(tf_score, cell = id)
mat <- as.data.frame(t(as.matrix(sub[["RNA"]]$data)))
mat$cell_type <- sub@meta.data[match(rownames(mat),rownames(sub@meta.data)),"ATAC.Subclass.June18_2025"]
mat$cell_type <- as.character(mat$cell_type)

# pseudobulk
num_cols <- colnames(mat)[1:883]
pseudobulk <- mat %>%
  group_by(cell_type) %>%
  summarise(across(all_of(num_cols), mean)) %>%
  column_to_rownames("cell_type")

# Filter the TF for just circadian TF
TF_net <- read.csv("TF_connectivity.csv")
TFs <- unique(TF_net$source)
pseudobulk <- as.data.frame(t(as.matrix(pseudobulk[,TFs])))

weight <- data.frame(TF = TFs)
for (i in seq(length(files))) {
  name <- sub("\\..*", "", files[i])
  df <- read.table(paste0("./Integration/",files[i]), sep = "\t", header = T)
  df[,"log1/p"] <- log(1/df$pvalue)
  weight[,i+1] <- df[match(weight$TF,df$GeneID),"log1/p"]
  colnames(weight)[i+1] <- name
}
weight[is.na(weight)] <- 0
rownames(weight) <- weight$TF
weight <- weight[,-1]
weight <- weight[match(rownames(pseudobulk),rownames(weight)),]
pseudobulk <- pseudobulk[,match(colnames(weight),colnames(pseudobulk))]
pseudobulk_weight <- pseudobulk*weight
colnames(pseudobulk_weight) <- paste(colnames(pseudobulk_weight),"weighted",sep = "_")
pseudobulk_weight <- t(pseudobulk_weight)
pseudobulk <- t(pseudobulk)
pseudo_all <- rbind(pseudobulk,pseudobulk_weight)

TF_net <- read.csv("TF_connectivity.csv")
for (ct in rownames(pseudo_all)) {
  tf_scores <- setNames(as.numeric(pseudo_all[ct, ]), colnames(pseudo_all))
  TF_net[[ct]] <- tf_scores[TF_net$source]
}

# Remove ENSxxxx TF
TF_net <- TF_net[-c(grep("ENS",TF_net$source),grep("ENS",TF_net$target)),]
write.csv(TF_net, file = "./network_weighted_zero_anno_celltype.csv",quote = F)