import scipy as sci
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
import scanpy.api as sc
import anndata
from anndata import h5py, logging
from anndata import AnnData, read_h5ad
import pandas as pd


def getAnnData(input_file):
    h5 = h5py.File(input_file,'r')
    data = h5['/GRCh38/data']
    indices = h5['/GRCh38/indices']
    barcodes = h5['/GRCh38/barcodes']
    indptr = h5['/GRCh38/indptr']
    genes = h5['/GRCh38/genes']
    gene_names = h5['/GRCh38/gene_names']
    shape = h5['/GRCh38/shape']

    X = csr_matrix((data,indices,indptr),shape=(shape.value[1],shape.value[0]),dtype='float32')
    adata = AnnData(X,
		   obs=pd.DataFrame(index=barcodes.value),
		   var=pd.DataFrame(index=genes.value),
		   dtype=X.dtype.name,
		   filemode=True)
    return adata
 
def getAnnData_10x_h5(input_file):
	adata = sc.read_10x_h5(input_file,"GRCh38")
	return adata

def getAnnData_10x_mtx(input_file):
	adata = sc.read_10x_mtx(input_file)
	return adata

def pre_process_input_data(gene_file,input_file,filtered=False,format_type="10x_mtx"):
	if format_type == "10x_h5":
		adata = getAnnData_10x_h5(input_file)
	else:
		adata = getAnnData_10x_mtx(input_file)
	if filtered:
		return adata
	rownames = adata.obs_names.values
	colnames = adata.var['gene_ids'].values
	X2=pd.DataFrame(adata.X.todense(), index=rownames, columns=colnames)
	
	gene_input = pd.read_csv(gene_file,index_col=1)
	
	common_ind = pd.Index(colnames).intersection(gene_input.index.values)
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
