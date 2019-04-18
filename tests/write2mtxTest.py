import sys
sys.path.insert(0,'./lib/pylib/')
import preprocess_10x as prep

ifile = sys.argv[1]
ifile2 = sys.argv[2]
path = sys.argv[3]
############## Read 10x datasets###################
adata = prep.read_10x_data(ifile,"10x_h5")
adata_cluster = prep.read_10x_data(ifile2,"10x_h5ad")
cellindex = adata_cluster.obs_names
adata=adata[cellindex,:]
# adata.obs = adata.obs.join(adata_cluster.obs)
adata.obs['louvain']=adata_cluster.obs['louvain']
adata.uns = adata_cluster.uns
adata.obsm = adata_cluster.obsm
prep.write2mtx(adata,path=path)

