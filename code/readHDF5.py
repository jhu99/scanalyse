#!/usr/bin/python3.5
import h5py
import numpy as np
import anndata
from scipy.sparse import csr_matrix
import pandas as pd

def my_read_hdf(path):
    f = h5py.File(path, 'r')
    group_GRCh38 = f['/GRCh38']
    cell_names = group_GRCh38['barcodes']
    mat_shape = group_GRCh38['shape']
    indptr = np.array(group_GRCh38['indptr'])
    indices = np.array(group_GRCh38['indices'])
    data = group_GRCh38['data']
    genes = group_GRCh38['genes']
    X=csr_matrix((data, indices, indptr), shape=(mat_shape[1], mat_shape[0])).toarray()
    print("aaaa")
    adata = anndata.AnnData(X, pd.DataFrame(index=cell_names.value), pd.DataFrame(index=genes.value), dtype=X.dtype.name)
    print("bbb")
    return adata

a=my_read_hdf('../data/ica_cord_blood_h5.h5')

