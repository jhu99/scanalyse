import sys
sys.path.insert(0,'./lib/pylib/')
import preprocess_10x as prep
import zinb_AE as zae
############## train model from ./test/train_model_test.py
############## Read reference panel in h5ad###################
ifile='./data/hca/ica_bm_after_filtered_recipe_zheng_qc1.h5ad'
ifile_raw='./data/hca/ica_bm_after_filtered_raw_recipe_zheng_qc1.h5ad'
outputpath="./result/ica_bm_qc1/"
adata = prep.read_10x_data(ifile,"10x_h5ad",'r')
adataraw = prep.read_10x_data(ifile_raw,"10x_h5ad")
adata.raw = adataraw.copy()
############# Predict latent layer#####################
zae.prediction(adata,outputpath)
