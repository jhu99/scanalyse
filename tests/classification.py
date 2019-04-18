import sys
sys.path.insert(0,'./lib/pylib/')
import preprocess_10x as prep
import zinb_AE as zae
import scanpy as sc
import matplotlib.pyplot as pl
import seaborn as sns
import pandas as pd
import os
sns.set(style='white', rc={'figure.figsize':(8,6), 'figure.dpi':150})
# Read projection after louvain and leiden clusters
method="t_test"
# path="./result/ica_bm_qc4/"
path=sys.argv[1]
cluster_type="louvain"

h5adfile_after_cluster=path+"ica_markers_"+method+".h5ad"
adata = prep.read_10x_data(h5adfile_after_cluster,"10x_h5ad")
num_clusters=len(adata.obs[cluster_type].astype('category').cat.categories)

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
df = pd.DataFrame(
	{group + '_' + key[:1]: result[key][group]
	for group in groups for key in ['names', 'pvals']}).head(10)
df.to_csv(path+"ica_markers_t_test_top10.csv")

top_ranked_genes=pd.DataFrame(adata.uns['rank_genes_groups']['names'][range(1)])
gene_index = pd.Index(top_ranked_genes.values.flatten()).drop_duplicates(keep='first')
gene_index = gene_index.union(pd.read_csv("./data/marker_genes.csv").x)
marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14', 'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1','FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
gene_index=pd.Index(marker_genes).union(gene_index)
gene_index=pd.Index(gene_index).intersection(adata.var_names)
ofile="_ica_louvain_marker_"+method+".png"
gs = sc.pl.matrixplot(adata, gene_index, groupby='louvain', dendrogram=True, show=False, save=ofile)
filename="matrixplot"+ofile
os.rename("./figures/"+filename,path+filename)

obs=adata.obs
cells_in_cluster=[]
cell_index=pd.Index([])
import random
for cluster in range(num_clusters):
	cells = obs.index[obs[cluster_type].astype(int)==cluster]
	cell_index = cell_index.union(pd.Index(random.sample(list(cells),200)))
	cells_in_cluster = cells_in_cluster + [len(cells)]
	print([cluster,len(cells)])
# cell_index=random.sample(list(adata.obs_names),int(0.01*len(adata.obs_names)))
adata=adata[cell_index,:]

pd.DataFrame(cells_in_cluster).to_csv(path+"ica_cluster_size.csv",index=False,header=False)
ofile="_plot_ica_louvain_marker_"+method+".png"
pl.subplot()
sc.pl.heatmap(adata,gene_index, use_raw=False, swap_axes=True, show=False, show_gene_labels=True, save=ofile, groupby='louvain', dendrogram=True)
filename="heatmap"+ofile
os.rename("./figures/"+filename,path+filename)

# 
# # Differential expression by comparing to the rest gronps
# sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)
# sc.pl.violin(adata, ['CST3', 'NKG7', 'PPBP'], groupby='louvain')
# #
# new_cluster_names = [
# 'CD4 T', 'CD14+ Monocytes',
# 'B', 'CD8 T',
# 'NK', 'FCGR3A+ Monocytes',
# 'Dendritic', 'Megakaryocytes']
# adata.rename_categories('louvain', new_cluster_names)
# ax = sc.pl.dotplot(adata, marker_genes, groupby='louvain')
# ax = sc.pl.stacked_violin(adata, marker_genes, groupby='louvain', rotation=90)
