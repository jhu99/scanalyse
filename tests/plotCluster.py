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
############## Read projection and plot clusters########
filepath="result/ica_bm_qc1/"
# adata = prep.read_10x_data(filepath+"projection.csv","10x_csv")
adata = prep.read_10x_data(filepath+"ica_clusters.h5ad","10x_h5ad")
zae.plotCluster(adata,filepath=filepath,dm_reduction=False)
############## Plot Clusters ##############################
# ifile="result/zheng/frozen_pbmc_donor_b/projection.csv"
# h5adfile="./result/zheng/frozen_pbmc_donor_b/frozen_pbmc_donor_b_clusters.h5ad"
# o1="result/figures/plot_umap_frozen_pbmc_donor_b.pdf"
# o2="result/figures/plot_tsne_frozen_pbmc_donor_b.pdf"
# adata = prep.read_10x_data(ifile,"10x_csv")
# zae.plotCluster(adata,h5adfile=h5adfile,umapfile=o1,tsnefile=o2)
# #
# ifile="result/zheng/fresh_68k_pbmc_donor_a/projection.csv"
# h5adfile="./result/zheng/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_clusters.h5ad"
# o1="result/figures/plot_umap_fresh_68k_pbmc_donor_a.pdf"
# o2="result/figures/plot_tsne_fresh_68k_pbmc_donor_a.pdf"
# adata = prep.read_10x_data(ifile,"10x_csv")
# zae.plotCluster(adata,h5adfile=h5adfile,umapfile=o1,tsnefile=o2)
############## Read Zheng's data in mtx#######################
# ifile = 'data/zheng/fresh_68k_pbmc_donor_a/hg19/'
# gfile = 'result/ica_all/gene_list.csv'
# ofile = 'result/zheng/fresh_68k_pbmc_donor_a/projection.csv'
#adata = prep.read_10x_data(ifile, format_type="10x_mtx")
#adata = prep.pre_process_input_data(adata, gfile)
