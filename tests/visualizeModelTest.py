import sys
sys.path.insert(0,'./lib/pylib/')
import preprocess_10x as prep
import zinb_AE as zae

############## Read Zheng's data in mtx#######################
ifile = 'data/zheng/fresh_68k_pbmc_donor_a/hg19/'
gfile = 'result/ica_all/gene_list.csv'
ofile = 'result/zheng/fresh_68k_pbmc_donor_a/projection.csv'
#adata = prep.read_10x_data(ifile, format_type="10x_mtx")
#adata = prep.pre_process_input_data(adata, gfile)
############## Read reference panel in h5ad###################
# ifile='./data/hca/ica_all_after_filtered_recipe_zheng.h5ad'
# ifile_raw='./data/hca/ica_all_after_filtered_raw_recipe_zheng.h5ad'
# adata = prep.read_10x_data(ifile,"10x_h5ad",'r')
# adataraw = prep.read_10x_data(ifile_raw,"10x_h5ad")
# adata.raw = adataraw.copy()
############## load pretrained weight #####################
# zae.prediction(adata,ofile)
############## Plot Clusters ##############################
ifile="result/zheng/frozen_pbmc_donor_b/projection.csv"
h5adfile="./result/zheng/frozen_pbmc_donor_b/frozen_pbmc_donor_b_clusters.h5ad"
o1="result/figures/plot_umap_frozen_pbmc_donor_b.pdf"
o2="result/figures/plot_tsne_frozen_pbmc_donor_b.pdf"
adata = prep.read_10x_data(ifile,"10x_csv")
zae.plotCluster(adata,h5adfile=h5adfile,umapfile=o1,tsnefile=o2)
#
ifile="result/zheng/fresh_68k_pbmc_donor_a/projection.csv"
h5adfile="./result/zheng/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_clusters.h5ad"
o1="result/figures/plot_umap_fresh_68k_pbmc_donor_a.pdf"
o2="result/figures/plot_tsne_fresh_68k_pbmc_donor_a.pdf"
adata = prep.read_10x_data(ifile,"10x_csv")
zae.plotCluster(adata,h5adfile=h5adfile,umapfile=o1,tsnefile=o2)
#
ifile="result/ica_all/latent.csv"
h5adfile="./result/ica_all/ica_clusters.h5ad"
o1="result/figures/plot_umap_ica.pdf"
o2="result/figures/plot_tsne_ica.pdf"
adata = prep.read_10x_data(ifile,"10x_csv")
zae.plotCluster(adata,h5adfile=h5adfile,umapfile=o1,tsnefile=o2)

