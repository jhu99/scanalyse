import sys
import os
sys.path.insert(0,'./lib/pylib/')
import preprocess_10x as prep

# b_cells_filtered_matrices_mex, cd14_monocytes_filtered_matrix_mex,
# cd56_nk_filtered_matrices_mex, 
# ifile='./data/hca/ica_bone_marrow_h5.h5'
# # ifile='./data/hca/ica_bone_marrow_filtered_gene_bc_matrices_mex2/hg19/'
# ofile='./data/hca/ica_bm_after_filtered_recipe_zheng_qc2.h5ad'
# ofilex='./data/hca/ica_bm_after_filtered_raw_recipe_zheng_qc2.h5ad'
# f='data/hca/ica_bone_marrow_filtered_gene_bc_matrices_mex2/hg19/barcodes_ann_revised.tsv'
ifile=sys.argv[1]
ofile=sys.argv[2]
ofilex=sys.argv[3]
path=sys.argv[4]
############## Read 10x datasets###################
adata = prep.read_10x_data(ifile,"10x_h5")
adataraw=adata.copy()
############## Filter and normalization###################
adata = prep.recipe_zheng(adata, n_top_genes=2000)
os.rename("./figures/violin_quality_control.png",path+"violin_quality_control.png")
cell_idx=adata.obs_names
gene_idx=adata.var_names
adata_subset = adataraw[:,gene_idx]
adata_subset = adata_subset[cell_idx,:]
adata.write(ofile,compression='gzip',force_dense=False)
adata_subset.write(ofilex,compression='gzip',force_dense=False)

