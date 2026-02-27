library(data.table)
library(tidyverse)
setwd("~/cirdian")

#motif <- as.data.frame(fread(file = "./pseudobulk_motif_mtx.csv.gz",header = T))
#motif <- unique(motif$V1)
#motif_df <- unique(data.frame(motif_name = gsub("\\_.*","",motif)))
#write.csv(motif_df, file = "motif_mapped.csv",row.names = F)

# Read in mapped motif
motif_df <- read.csv("./motif_mapped.csv")
# Use motif name for those do not have gene_symbol
motif_df$gene_symbol[which(motif_df$gene_symbol %in% "")] <- motif_df$motif_name[which(motif_df$gene_symbol %in% "")]
# Read in gene expression
geneExp <- as.data.frame(fread(file = "./pseudobulk_GE_mtx.csv.gz",header = T))
# For those cannot be mapped, use motif name
motif_df$gene_symbol[-which(motif_df$gene_symbol %in% geneExp$V1)] <- motif_df$motif_name[-which(motif_df$gene_symbol %in% geneExp$V1)]
# Filter gene expression
geneExp <- geneExp[which(geneExp$V1 %in% motif_df$gene_symbol),]
colnames(geneExp)[1] <- "gene_symbol"


motif_df <- motif_df %>%
  left_join(geneExp, by = "gene_symbol")
motif_df <- na.omit(motif_df)
df_TF_exp <- data.frame()
for (i in 3:10) {
  celltype <- colnames(motif_df)[i]
  df_tmp <- motif_df[,c(1,2,i)]
  df_tmp <- df_tmp[order(df_tmp[[3]], decreasing = T),]
  df_tmp <- df_tmp[which(df_tmp[[3]] != 0),]
  df_tmp <- df_tmp[1:round(0.9*nrow(df_tmp),0),]
  colnames(df_tmp)[3] <- "Expression"
  df_tmp$celltype <- rep(celltype,nrow(df_tmp))
  df_TF_exp <- rbind(df_TF_exp,df_tmp)
}
#write.csv(df_TF_exp, file = "Left_TF_rank_GEX.csv", row.names = F)

# TF_network
TF_network <- read.csv("./network_weighted_zero_anno_celltype.csv")
TF_network <- unique(TF_network[,c(1,grep("weighted",colnames(TF_network)))])
#write.csv(TF_network, file = "TF_network_mapp.csv", row.names = F)

TF_network <- TF_network[which(TF_network$source %in% unique(df_TF_exp$motif_name)),]
rownames(TF_network) <- TF_network$source
TF_network <- TF_network[,-1]
colnames(TF_network) <- gsub("_weighted","",colnames(TF_network))
df_TF_exp <- df_TF_exp[which(df_TF_exp$motif_name %in% rownames(TF_network)),]

df_TF_exp_wide <- df_TF_exp %>%
  pivot_wider(names_from = celltype, values_from = Expression, values_fill = 0, values_fn = mean) %>%
  column_to_rownames(var = "motif_name")

TF_network <- TF_network[rownames(df_TF_exp_wide),colnames(df_TF_exp_wide)[-1]]
df_TF_exp_wide <- df_TF_exp_wide[,-1]
df_TF_exp_wide[df_TF_exp_wide > 0] <- 1 # Set the gene expression >0 to 1

# Multiple the enrichment score and the binarized expression
g.w.TF_network <- TF_network*df_TF_exp_wide
# Remove the TFs having 0 for all cell types
g.w.TF_network <- g.w.TF_network[rowSums(g.w.TF_network) >0,] #  TFs left 303

# Match the distance and connectivity
TF_network <- read.csv("./network_weighted_zero_anno_celltype.csv")
TF_network <- TF_network[,which(colnames(TF_network) %in% c("source","target","connectivity","distance"))]
TF_network <- TF_network[which(TF_network$source %in% rownames(g.w.TF_network)),]
TF_network <- TF_network[which(TF_network$target %in% rownames(g.w.TF_network)),]
g.w.TF_network <- t(g.w.TF_network)
rownames(g.w.TF_network) <-  paste(rownames(g.w.TF_network), "weighted",sep = "_")

# Map the weighted enrichment score to source
for (ct in rownames(g.w.TF_network)) {
  tf_scores <- setNames(as.numeric(g.w.TF_network[ct, ]), colnames(g.w.TF_network))
  TF_network[[ct]] <- tf_scores[TF_network$source]
}
write.csv(TF_network, file = "./network_motif_weighted_zero_celltype_filtered.csv",quote = F)
