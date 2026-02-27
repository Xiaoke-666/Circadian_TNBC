#!/usr/bin/env python3

import pandas as pd
import scanpy as sc
import numpy as np
import os
import sys
import muon as mu
from muon import atac as ac
from matplotlib.pyplot import rc_context
from pandas.api.types import CategoricalDtype
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import rc_context

# Read in the data
os.chdir("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP")
sc.figure_size=(8,8)
sc.set_figure_params(scanpy=True, dpi=300, dpi_save=300,fontsize=6, color_map="viridis")
pd.options.display.max_rows = 999

def pl_umap(adata,lsts,prefix):
    for lst in lsts:
        sc.pl.umap(adata, color=lst, legend_loc='on data',legend_fontsize=3,
                    save="."+prefix+"."+lst+".ondata.png")
        sc.pl.umap(adata, color=lst, save="."+prefix+"."+lst+".png")


def pl_umap2(adata,lsts,prefix):
    for lst in lsts:
        sc.pl.umap(adata, color=lst, 
                   vmin='p5',vmax='p95',
                   save="."+prefix+"."+lst+".png")


def pl_umap3(adata,lsts,prefix):
    for lst in lsts:
        print(lst)
        sc.pl.umap(adata, color=lst, 
                   frameon=False,title="",
                   legend_loc=None,colorbar_loc=None,
                   cmap="coolwarm",
                   vmin="p5",vmax="p95", 
                   show=False)
        plt.savefig("./figures/UMAP."+prefix+"."+lst+".png", 
                    dpi=300, bbox_inches='tight', transparent=True)
        


BuRd=["#2166ACFF", "#4393C3FF", "#92C5DEFF", "#D1E5F0FF", "#F7F7F7FF", "#FDDBC7FF", "#F4A582FF", "#D6604DFF", "#B2182BFF"]

#RdBu_r
def pl_umap3(adata,lsts,prefix):
    for lst in lsts:
        sc.pl.umap(adata, color=lst, 
                   vmin='p20',vmax='p95',frameon=False,
                   color_map="RdBu_r",
                   save="."+prefix+"."+lst+".png")

rna3=sc.read("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP/rna4.Cellbender.ATAC_UMAP.Overlap.Apr22_2025.h5ad")

rna3.X=rna3.layers["raw"].copy()
for x in ['GeneSetScore_CircadianModule.M3', 'GeneSetScore_CircadianModule.M12', 'GeneSetScore_CircadianModule.M14', 'GeneSetScore_CircadianModule.M4', 'GeneSetScore_CircadianModule.M7', 'GeneSetScore_CircadianModule.M1', 'GeneSetScore_CircadianModule.M9', 'GeneSetScore_CircadianModule.M5', 'GeneSetScore_CircadianModule.M10', 'GeneSetScore_CircadianModule.M6', 'GeneSetScore_CircadianModule.M2', 'GeneSetScore_CircadianModule.M11', 'GeneSetScore_CircadianModule.M15', 'GeneSetScore_CircadianModule.M8', 'GeneSetScore_CircadianModule.M16', 'GeneSetScore_CircadianModule.M13', 'GeneSetScore_CircadianModule.M17', 'GeneSetScore_CircadianModule.M18', 'GeneSetScore_CircadianModule.M19', 'GeneSetScore_CircadianModule.M20']:
    del rna3.obs[x]

del rna3.layers["raw"]


celltypes=['Cancer_Epi',
 'CAFs',
 'CD4_T',
 'B',
 'Macro_Mono',
 'CD8_T',
 'NKT',
 'Neu'
]

atac3=sc.read("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/02_snATAC/04_Genescore_Matrix/GeneScoreMatrix_combined.h5ad")
atac3.var_names= atac3.var["name"]
atac3.var_names_make_unique()
# index name as None
atac3.var.index.name = None


intersected_genes= set(atac3.var_names).intersection(set(rna3.var_names))

intersected_genes= list(intersected_genes)
# write intersected_genes to a tsv
intersected_genes_df= pd.DataFrame(intersected_genes, columns=["gene"])
intersected_genes_df.to_csv("ATAC_RNA.intersected_genes.tsv", sep="\t", index=False)
rna3 = rna3[:,intersected_genes]

for celltype in celltypes:
    print(celltype)
    tmp= rna3[rna3.obs["ATAC.Subclass.Nov23_2024"]==celltype].copy()
    circadian_genes= pd.read_csv(f"/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_circadian/03_snMultiome_Tumor_Jan19_2024/01_snRNA_cellbender/01_map_snATAC_UMAP/CircadianGenes_v2/Integration/{celltype}.discoODAres_CS_ANNOATION.outer.tsv",
                            sep="\t",index_col=0)
    sig_circadian_genes= circadian_genes[circadian_genes["pvalue"]<0.02]
    sig_circadian_genes=sig_circadian_genes.sort_values(by="acrophase", ascending=True)
    sig_circadian_genes= sig_circadian_genes[sig_circadian_genes.index.isin(intersected_genes)]
    tmp.var["CircadianGene"]= tmp.var_names.isin(sig_circadian_genes.index)
    sc.pp.calculate_qc_metrics(tmp, qc_vars=["CircadianGene"], inplace=True, percent_top=[], log1p=True)
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    lst= []
    nCells=tmp.obs["TimePoint"].value_counts().min()
    print(f"Number of cells : {nCells}")
    print(f"Number of circadian genes : {len(sig_circadian_genes.index)}")

    for tp in ['ZT4', 'ZT10','ZT16' , 'ZT22']:
        print(tp)
        adata=tmp[tmp.obs["TimePoint"]==tp].copy()
        #adata=adata[adata.obs["total_counts_CircadianGene"]>0, :]
        adata=adata[np.argsort(adata.obs["total_counts_CircadianGene"])[-nCells:], :]
        lst=lst + adata.obs.index.tolist()
    tmp2= tmp[lst,sig_circadian_genes.index].copy()
    print(f"writing {celltype} data with {len(tmp2.obs.index)} cells and {len(tmp2.var_names)} genes")
    tmp2.write(f"rna4.Cellbender.{celltype}.CircadianGenes.ZT4_10_16_22.top{nCells}.h5ad")
    df= pd.DataFrame(tmp2.X.T.toarray(), index=tmp2.var_names.tolist(),
                        columns=tmp2.obs.index.tolist())
    df2= df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    #df2.to_csv(f"rna4.Cellbender.{celltype}.CircadianGenes.ZT4_10_16_22.top{nCells}.tsv.gz", sep="\t", compression="gzip")
    # rank cells
    time_order = ["ZT4", "ZT10", "ZT16", "ZT22"]
    tmp2.obs["TimePoint"] = pd.Categorical(tmp2.obs["TimePoint"], categories=time_order, ordered=True)
    top_genes=sig_circadian_genes[sig_circadian_genes["acrophase"]<6].index.tolist()
    ranks = []
    for tp in time_order:
        print(f"Ranking cells in {tp}")
        cell_ids = tmp2.obs[tmp2.obs["TimePoint"] == tp].index
        sub_df = df2.loc[top_genes, cell_ids]

        # 用top基因表达的总和排序细胞
        sorted_cells = sub_df.sum(axis=0).sort_values(ascending=False).index.tolist()
        
        # 添加到整体顺序中
        ranks.extend(sorted_cells)
    # 排序矩阵列
    df2_sorted = df2[ranks]
    print(f"Writing sorted data for {celltype} with {len(df2_sorted.columns)} cells and {len(df2_sorted.index)} genes")
    df2_sorted.to_csv(f"rna4.Cellbender.{celltype}.CircadianGenes.ZT4_10_16_22.top{nCells}.sorted.tsv.gz", sep="\t", compression="gzip")




'''


from scipy.ndimage import gaussian_filter1d

smoothed_df = CD4_T2_df_sorted.apply(lambda x: gaussian_filter1d(x, sigma=2), axis=1, result_type='broadcast')
smoothed_df.to_csv("CD4_T.CircadianGenes.ZT4_10_16_22.top200.cells.smoothed.tsv.gz", sep="\t", compression="gzip")


from statsmodels.nonparametric.smoothers_lowess import lowess
import pandas as pd
import numpy as np

def loess_smooth_row(row, frac=0.1):
    x = np.arange(len(row))
    y = row.values
    # frac 是相当于 span，越大越平滑
    smoothed = lowess(y, x, frac=frac, return_sorted=False)
    return pd.Series(smoothed, index=row.index)

# 应用到整个表达矩阵
smoothed_df = CD4_T2_df_sorted.apply(loess_smooth_row, axis=1)

smoothed_df = smoothed_df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)


smoothed_df.to_csv("CD4_T.CircadianGenes.ZT4_10_16_22.top200.cells.lowess_smoothed.tsv.gz", sep="\t", compression="gzip")

'''


############## 2000 Cells 

lst=[]
for tp in ['ZT4', 'ZT10','ZT16' , 'ZT22']:
    print(tp)
    adata=CD4_T[CD4_T.obs["TimePoint"]==tp].copy()
    adata=adata[adata.obs["total_counts_CircadianGene"]>0, :]
    adata=adata[np.argsort(adata.obs["total_counts_CircadianGene"])[-2000:], :]
    lst=lst + adata.obs.index.tolist()

# normalize
CD4_T2= CD4_T[lst,sig.index].copy()

# save the data
CD4_T2_df= pd.DataFrame(CD4_T2.X.T.toarray(), index=CD4_T2.var_names.tolist(),
                        columns=CD4_T2.obs.index.tolist())

# scale by row
CD4_T2_df= CD4_T2_df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

#CD4_T2_df.to_csv("CD4_T.CircadianGenes.ZT4_10_16_22.top2000.cells.tsv.gz", sep="\t", compression="gzip")

from scipy.ndimage import gaussian_filter1d
#rank columns for 1:200, 201:400, 401:600, 601:800, 801:1000, by top 100 row from highest to lowest
time_order = ["ZT4", "ZT10", "ZT16", "ZT22"]
CD4_T2.obs["TimePoint"] = pd.Categorical(CD4_T2.obs["TimePoint"], categories=time_order, ordered=True)

top100_genes = CD4_T2_df.index[:100]  # 前100行

ranks = []

for tp in time_order:
    print(f"Ranking cells in {tp}")
    cell_ids = CD4_T2.obs[CD4_T2.obs["TimePoint"] == tp].index
    sub_df = CD4_T2_df.loc[top100_genes, cell_ids]

    # 用前100个基因表达的总和排序细胞
    sorted_cells = sub_df.sum(axis=0).sort_values(ascending=False).index.tolist()
    
    # 添加到整体顺序中
    ranks.extend(sorted_cells)

# 排序矩阵列
CD4_T2_df_sorted = CD4_T2_df[ranks]
CD4_T2_df_sorted.to_csv("CD4_T.CircadianGenes.ZT4_10_16_22.top2000.cells.sorted.tsv.gz", sep="\t", compression="gzip")






############## 3000 Cells 

lst=[]
for tp in ['ZT4', 'ZT10','ZT16' , 'ZT22']:
    print(tp)
    adata=CD4_T[CD4_T.obs["TimePoint"]==tp].copy()
    adata=adata[adata.obs["total_counts_CircadianGene"]>0, :]
    adata=adata[np.argsort(adata.obs["total_counts_CircadianGene"])[-3000:], :]
    lst=lst + adata.obs.index.tolist()

# normalize
CD4_T2= CD4_T[lst,sig.index].copy()

# save the data
CD4_T2_df= pd.DataFrame(CD4_T2.X.T.toarray(), index=CD4_T2.var_names.tolist(),
                        columns=CD4_T2.obs.index.tolist())

# scale by row
CD4_T2_df= CD4_T2_df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

#CD4_T2_df.to_csv("CD4_T.CircadianGenes.ZT4_10_16_22.top2000.cells.tsv.gz", sep="\t", compression="gzip")

from scipy.ndimage import gaussian_filter1d
#rank columns for 1:200, 201:400, 401:600, 601:800, 801:1000, by top 100 row from highest to lowest
time_order = ["ZT4", "ZT10", "ZT16", "ZT22"]
CD4_T2.obs["TimePoint"] = pd.Categorical(CD4_T2.obs["TimePoint"], categories=time_order, ordered=True)

top100_genes = CD4_T2_df.index[:100]  # 前100行

ranks = []

for tp in time_order:
    print(f"Ranking cells in {tp}")
    cell_ids = CD4_T2.obs[CD4_T2.obs["TimePoint"] == tp].index
    sub_df = CD4_T2_df.loc[top100_genes, cell_ids]

    # 用前100个基因表达的总和排序细胞
    sorted_cells = sub_df.sum(axis=0).sort_values(ascending=False).index.tolist()
    
    # 添加到整体顺序中
    ranks.extend(sorted_cells)

# 排序矩阵列
CD4_T2_df_sorted = CD4_T2_df[ranks]
CD4_T2_df_sorted.to_csv("CD4_T.CircadianGenes.ZT4_10_16_22.top3000.cells.sorted.tsv.gz", sep="\t", compression="gzip")
