import sys
sys.path.insert(0,'./lib/pylib/')
import preprocess_10x as prep

ifile='./data/hca/ica_all.h5'
ofile='./data/hca/ica_all_after_filtered_recipe_zheng.h5ad'
ofilex='./data/hca/ica_all_after_filtered_raw_recipe_zheng.h5ad'
############## Read 10x datasets###################
adata = prep.read_10x_data(ifile,"10x_h5")
adataraw=adata.copy()
############## Filter and normalization###################
adata = prep.recipe_zheng(adata)
cell_idx=adata.obs_names
gene_idx=adata.var_names
adata_subset = adataraw[:,gene_idx]
adata_subset = adata_subset[cell_idx,:]
adata.write(ofile,compression='gzip',force_dense=False)
adata_subset.write(ofilex,compression='gzip',force_dense=False)

