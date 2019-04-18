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
# path="./result/ica_bm_qc3/"
path=sys.argv[1]
cluster_type="louvain"
h5adfile_after_cluster=path+"ica_markers_"+method+".h5ad"
adata = prep.read_10x_data(h5adfile_after_cluster,"10x_h5ad")

sc.pl.tsne(adata,color=['LYZ','CST3','CD14','MS4A7','FCGR3A','FCER1A'],show=False,ncols=3)
pl.savefig(path+"ica_tsne_marker1.png")
pl.close()
sc.pl.tsne(adata,color=['GNLY','NKG7','KLRB1','IL7R','CD8B','CD27'],show=False,ncols=3)
pl.savefig(path+"ica_tsne_marker2.png")
pl.close()
sc.pl.tsne(adata,color=['CD79A','MME','MS4A1','SEPP1','SDC1','MZB1'],show=False,ncols=3)
pl.savefig(path+"ica_tsne_marker3.png")
pl.close()
sc.pl.tsne(adata,color=['PPBP','IL3RA','HBB','CD34','PF4','RPL34'],show=False,ncols=3)
pl.savefig(path+"ica_tsne_marker4.png")
pl.close()

sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True,show=False)
pl.savefig(path+"violin_after_quality_control.png")
pl.close()

sc.pl.tsne(adata,color=['louvain'],show=False)
# 	#pl.title("Visualization of ~700K HCA immune cells via t-SNE")
pl.title("")
pl.legend(loc=3,fontsize=6,mode="expand",bbox_to_anchor=(0.0, 1.01, 1, 0.2),ncol=17)
pl.savefig(path+"ica_tsne_louvain.png")
pl.close()

# new_cluster_names = [
# 	'Mature B Cells 1','CD14+ Monocytes 1','DNT 1','CD8+ CCR7-','NK Cells 1',
# 	'DNT 2','NK Cells 2','cDC 1','CD8+ CM','MT-hi Cells',
# 	'CD8+ naïve', 'Immature B Cells 1', 'HBB','CD14+ Monocytes 2','NK Cells 3',
# 	'Non-classical Monocytes', 'Immature B Cells 2', 'HSPC', 'cDC 2', 'Plasma B Cell',
# 	'CD14+ Monocytes 3', 'pDC', 'Mature B Cells 2', 'SEPP1', 'Megakaryocytes', 'NK Cells 4']
# adata.rename_categories('louvain', new_cluster_names)
# sc.pl.tsne(adata, color='louvain', legend_loc='on data', title='', legend_fontsize=6, show=False)
# pl.savefig(path+"ica_tsne_louvain_annotated.png")
# pl.close()

marker_genes = ['CD14','CD79A','LYZ','NKG7','GNLY','CD8B','PPBP','SEPP1','HBB','CD34','SDC1','MME']
ax = sc.pl.stacked_violin(adata, marker_genes, groupby='louvain', dendrogram=True, swap_axes=True,show=False)
pl.savefig(path+"ica_violin_louvain_marker_genes.png")
pl.close()

# sc.pl.rank_genes_groups_heatmap(adata,groupby='louvain',show=False)
# pl.savefig(path+"ica_heatmap_louvain_genes_group.png")
# pl.close()


