import sys
sys.path.insert(0,'./lib/pylib/')
import preprocess_10x as prep

ifile='./data/hca/ica_bone_marrow_h5.h5'
path='./data/hca/ica_bone_marrow_filtered_gene_bc_matrices_mex2/hg19/'
############## Read 10x datasets###################
adata = prep.read_10x_data(ifile,"10x_h5")
prep.filter_basic(adata)
prep.write2mtx(adata,path=path)

