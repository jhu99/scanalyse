from autoEncoder import AutoEncoder

import h5py
from math import fabs
import numpy as np
from scipy.sparse import csc_matrix,csr_matrix

f = h5py.File("../../data/mtx_rank01.h5","r")
print(f)
genes = f['/GRCh38/genes'][:]
barcodes = f['/GRCh38/barcodes'][:]
data = np.array(f['/GRCh38/data'][:])
indices=np.array(f['/GRCh38/indices'][:])
indptr=np.array(f['/GRCh38/indptr'][:])
rank_zero=np.array(f['/GRCh38/zero_value'][:])
print(len(barcodes))
print(len(genes))
print(len(data))
sparse_matrix=csr_matrix((data, indices, indptr), shape=(len(barcodes), len(genes))).toarray()
for i in range(len(barcodes)):
    for j in range(len(genes)):
        #print(sparse_matrix[i]," ")
        if not (sparse_matrix[i][j] > 0):
            sparse_matrix[i][j]=rank_zero[i]

autoEncoder = AutoEncoder(len(genes))
autoEncoder.build(sparse_matrix)
