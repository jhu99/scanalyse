import sys
sys.path.insert(0,'./lib/pylib/')
import preprocess_10x as prep
import zinb_AE as net

# ifile='./data/hca/ica_bm_after_filtered_recipe_zheng_qc1.h5ad'
# ifile_raw='./data/hca/ica_bm_after_filtered_raw_recipe_zheng_qc1.h5ad'
# path='./result/ica_bm_qc1/'
ifile=sys.argv[1]
ifile_raw=sys.argv[2]
path=sys.argv[3]
############## Read h5ad###################
adata = prep.read_10x_data(ifile,"10x_h5ad",'r')
adataraw = prep.read_10x_data(ifile_raw,"10x_h5ad")
############## Train model#####################
adata.raw=adataraw.copy()
adata.var.to_csv(path+'genelist.csv')
net.train_zinb_model(adata, filepath=path)

