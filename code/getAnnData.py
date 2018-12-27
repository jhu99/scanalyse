import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
import anndata
from anndata import h5py, logging
from anndata import AnnData, read_h5ad

def getAnnData(path):
    h5 = h5py.File(path,'r')
    data = h5['/GRCh38/data']
    indices = h5['/GRCh38/indices']
    barcodes = h5['/GRCh38/barcodes']
    indptr = h5['/GRCh38/indptr']
    gene_names = h5['/GRCh38/gene_names']
    shape = h5['/GRCh38/shape']

    X = csr_matrix((data,indices,indptr),shape=(len(barcodes),len(gene_names)))
    adata = AnnData(X,obs=np.array(barcodes),var=np.array(gene_names),dtype=X.dtype.name, filemode=True)
   
    return adata
