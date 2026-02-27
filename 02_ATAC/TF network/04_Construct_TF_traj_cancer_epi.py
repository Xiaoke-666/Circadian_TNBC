import scanpy as sc
import pandas as pd
import os
import numpy as np
import glob

#########
## For cancer epi
#########

f_mat = sc.read_h5ad("./cirdian/ChromVar_Deviations.Clean.T.July17_2025.h5ad")
# Read in data and find TFs
df_cancer = pd.read_table("./cirdian/Cancer_Epi.discoODAres_CS_ANNOATION.outer.pvalue_0.02.tsv")
# Read in GEX filtered TF
tf_filtered = pd.read_csv("./cirdian/Left_TF_rank_GEX.csv")
tf_filtered = tf_filtered.loc[tf_filtered['celltype']=='Cancer_Epi',:]
# Filter the cancer circadian
df_cancer = df_cancer.loc[df_cancer["GeneID"].isin(tf_filtered["motif_name"]),:] # 109 TF lefts

TFs = df_cancer["GeneID"]

# Remove duplicated names in tf_mat.obs
tf_mat = tf_mat[~tf_mat.obs_names.duplicated()].copy()

# subset tf-matrix with circadian genes
sub_cancer = tf_mat[tf_mat.obs['TF'].isin(TFs),:]
# Subset variables to just cancer epi
sub_cancer = sub_cancer[:,sub_cancer.var["ATAC.Subclass.June18_2025"].isin(["Cancer_Epi"])]

# Calculate neighbors
sc.pp.highly_variable_genes(sub_cancer, min_mean=0.02, max_mean=4, min_disp=0.5)
sc.pp.pca(sub_cancer)
sc.pp.neighbors(sub_cancer, use_rep="X_pca")

# Export connectivity
conn = sub_cancer.obsp["connectivities"]
dis = sub_cancer.obsp["distances"]
#Convert format
coo_conn = conn.tocoo()
coo_dist = dis.tocoo()
edges = pd.DataFrame({
    "source": [sub_cancer.obs_names[i] for i in coo_conn.row],
    "target": [sub_cancer.obs_names[j] for j in coo_conn.col],
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
edges_anno  = edges.merge(df_cancer, left_on ="source", right_on = "GeneID", how = "left")
edges_anno = edges_anno[["source","target","connectivity","distance","acrophase","pvalue"]]

# Make the acrophase discrete
bins = [0, 3, 6, 9, 12, 15, 18, 21, 24]
labels = ["P1","P2","P3","P4","P5","P6","P7","P8"]
edges_anno["phase"] = pd.cut(edges_anno["acrophase"], bins = bins, labels=labels, right = False)

# Add acrphase information to sub_cancer
tf_phase = edges_anno[["source","phase"]].drop_duplicates()
tf_phase = tf_phase[tf_phase["source"].isin(sub_cancer.obs["TF"])]
phase_map = dict(zip(tf_phase["source"], tf_phase["phase"]))
sub_cancer.obs["phase"] = sub_cancer.obs["TF"].map(phase_map)

# Rerun paga
sc.tl.paga(sub_cancer, groups = "phase")
sc.pl.paga(sub_cancer, plot=False) 
sc.tl.umap(sub_cancer, init_pos = "paga", spread = 1.2)

sub_cancer.write_h5ad("./cirdian/Cancer_circadian.h5ad", compression = "gzip")

# Plot
from matplotlib import rc_context
import matplotlib.pyplot as plt
from adjustText import adjust_text
sub_cancer.uns["phase_colors"] = ["#FF4545","#FF9C73","#FBD288","#FCF596","#DFF2EB",
                                  "#B9E5E8","#7AB2D3","#4A628A"]
TFs_label = ["Jun","Junb","Jund","Fos","Fosb","Fosl1","Tead3", "Tead4","Smad2", "Smad3", 
  "Ppara","Ppard","Pparg","Thra","Thrb","Rarg","Nfkb1", "Rel", "Rela", "Relb","Bcl6", 
  "Bach1", "Bach2","Nfe2l2","Nr1d1", "Nr1d2","Rora", "Rorb"]
with rc_context({"figure.figsize": (5, 4)}):
    sc.pl.paga_compare(sub_cancer,edges=True,title = "", color = "phase", right_margin=0,
                       frameon=False, legend_loc = "none",edge_width_scale=0.3, show = False,
                       node_size_scale=3.0,fontsize=15,legend_fontsize=15, size = 150.0)
    

fig = plt.gcf()
axes_list = fig.axes
ax_umap = axes_list[-2]

# Modify scatter point edge width
for coll in ax_umap.collections:
    coll.set_linewidth(0.1)      # edge width
    coll.set_edgecolor("#C0C0C0")

ax_paga = axes_list[-1]

from matplotlib.collections import LineCollection

for artist in ax_paga.collections:
    if isinstance(artist, LineCollection):
        artist.set_color("gray")     # edge color
        artist.set_alpha(0.8)        # optional


coords = sub_cancer.obsm["X_umap"]
obs_names = np.array(sub_cancer.obs_names)
texts = []
# Add labels for specific TFs
for tf in TFs_label:
    if tf in obs_names:
        i = np.where(obs_names == tf)[0][0]
        x, y = coords[i, 0], coords[i, 1]
        txt = ax_umap.text(x, y, tf, fontsize=10, weight="bold", color="black")
        texts.append(txt)

adjust_text(
    texts,
    ax=ax_umap,
    arrowprops=dict(arrowstyle="-", color='gray', lw=0.5),
    force_points=0.5,
    force_text=0.5,
    expand_text=(1.05, 1.2),
    expand_points=(1.05, 1.2)
)

plt.savefig("./cirdian/umap.pdf", bbox_inches="tight")