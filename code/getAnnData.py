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
