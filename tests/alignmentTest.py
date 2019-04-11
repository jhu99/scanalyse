import sys, os
sys.path.insert(0,'./lib/pylib/')
import preprocess_10x as prep
import numpy as np
from anndata import AnnData

datapath=os.listdir("./data/zheng/")
datapath.remove("293t_filtered_gene_bc_matrices_mex")
datapath.remove("frozen_pbmc_donor_b")
datapath.remove("frozen_pbmc_donor_c")
datapath.remove("fresh_68k_pbmc_donor_a")

for i in range(len(datapath)):
	file= "./data/zheng/"+datapath[i]+"/hg19/"
	adata = prep.read_10x_data(file,"10x_mtx")
	adata.obs['cell_type']=datapath[i]
	if i==0:
		X=adata.X.toarray()
		obs=adata.obs
		var=adata.var
	else:
		X= np.vstack((X,adata.X.toarray()))
		obs = obs.append(adata.obs)

adata = AnnData(X,
			obs=obs,
			var=var,
			dtype=X.dtype.name,
			filemode=True)

adata.write("./data/zheng/ten_mixed_cell_types.h5ad",compression="gzip")
