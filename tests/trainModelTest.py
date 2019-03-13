import sys
sys.path.insert(0,'./lib/pylib/')
import preprocess_10x as prep
import zinb_AE as net

ifile='./data/hca/ica_all_after_filtered_recipe_zheng.h5ad'
ifile_raw='./data/hca/ica_all_after_filtered_raw_recipe_zheng.h5ad'
#ofile='./result/'
############## Read h5ad###################
adata = prep.read_10x_data(ifile,"10x_h5ad",'r')
adataraw = prep.read_10x_data(ifile_raw,"10x_h5ad")
############## Train model#####################
adata.raw=adataraw.copy()
net.train_zinb_model(adata)

