import scanpy as sc
import pandas as pd
import os
import numpy as np
import glob

#########
## For all cell types
#########
# Read in TF matrix
tf_mat = sc.read_h5ad("ChromVar_Deviations.Clean.T.July17_2025.h5ad")
# Read in circadian genes
files = glob.glob("*.tsv")
df_list = [pd.read_table(f, sep="\t") for f in files]
df_comb = pd.concat(df_list, ignore_index=True)
gene = df_comb["GeneID"].unique()

# subset tf-matrix with circadian TF
sub_tf = tf_mat[tf_mat.obs_names.isin(gene)]
sub_tf = sub_tf[~sub_tf.obs_names.duplicated()].copy()

# Calculate neighbors
sc.pp.highly_variable_genes(sub_tf, min_mean=0.02, max_mean=4, min_disp=0.5)
sc.pp.pca(sub_tf)
sc.pp.neighbors(sub_tf, use_rep="X_pca")

conn = sub_tf.obsp["connectivities"]
dis = sub_tf.obsp["distances"]

#Convert format
coo_conn = conn.tocoo()
coo_dist = dis.tocoo()

edges = pd.DataFrame({
    "source": [sub_tf.obs_names[i] for i in coo_conn.row],
    "target": [sub_tf.obs_names[j] for j in coo_conn.col],
    "connectivity": coo_conn.data
})

dist_dict = {(i, j): d for i, j, d in zip(coo_dist.row, coo_dist.col, coo_dist.data)}

edges["distance"] = [
    dist_dict.get((i, j), np.nan) 
    for i, j in zip(coo_conn.row, coo_conn.col)
]

# Remove duplicates and self-loops
edges = edges[edges["source"] != edges["target"]]
edges = edges.drop_duplicates(subset=["source", "target"])

edges.to_csv("TF_connectivity.csv", index=False)

