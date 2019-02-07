import scipy as sci
from anndata import AnnData
import numpy as np
import pandas as pd
from getAnnData import getAnnData, getAnnData_10x_h5, getAnnData_10x_mtx

def pre_process_input_data(gene_file,input_file,format_type="10x_mtx"):
	if format_type == "10x_h5":
		adata = getAnnData_10x_h5(input_file)
	else:
		adata = getAnnData_10x_mtx(input_file)
	
	rownames = adata.obs_names.values
	colnames = adata.var_names.values
	X2=pd.DataFrame(adata.X.todense(), index=rownames, columns=colnames)
	
	gene_input = pd.read_csv(gene_file,index_col=1)
	
	
	common_ind = adata.var.index.intersection(gene_input.index.values)
	left_ind = gene_input.index.difference(common_ind)
	X3=X2[common_ind]
	for ind in left_ind:
		X3.insert(0,ind,0)
	X3=X3.reindex(columns=gene_input.index.values)
	X=sci.sparse.csr_matrix(X3,dtype='float32')
	adata = AnnData(X,
				obs=pd.DataFrame(index=rownames),
				var=pd.DataFrame(index=gene_input.index.values),
				dtype=X.dtype.name,
				filemode=True)
	return adata


