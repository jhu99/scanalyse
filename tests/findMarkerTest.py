import sys
sys.path.insert(0,'./lib/pylib/')
import preprocess_10x as prep
import zinb_AE as zae
import scanpy as sc
import scanpy.api
import matplotlib.pyplot as pl
import seaborn as sns
import random
sns.set(style='white', rc={'figure.figsize':(8,6), 'figure.dpi':150})

# h5adfile="./data/hca/ica_bm_after_filtered_recipe_zheng_qc3.h5ad"
# filepath="./result/ica_bm_qc3/"
# # # Read projection after louvain and leiden clusters
# h5adfile_after_cluster=filepath+"ica_clusters.h5ad"
h5adfile=sys.argv[1]
filepath=sys.argv[2]
h5adfile_after_cluster=sys.argv[3]

adata_cluster = prep.read_10x_data(h5adfile_after_cluster,"10x_h5ad")
adata= prep.read_10x_data(h5adfile,"10x_h5ad")
# Join the cluster annotation
adata.obs = adata.obs.join(adata_cluster.obs)
adata.uns = adata_cluster.uns
adata.obsm = adata_cluster.obsm
# Rank cluster specific genes using t-test
sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
#sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
adata.write(filepath+"ica_markers_t_test.h5ad")
# sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
# #sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
# adata.write(filepath+"ica_markers_wilcoxon.h5ad")
# sc.tl.rank_genes_groups(adata, 'louvain', method='logreg')
# #sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
# adata.write(filepath+"ica_markers_logreg.h5ad")



