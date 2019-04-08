from scipy.sparse import csr_matrix
from scipy.io import mmwrite
from anndata import h5py
import numpy as np


class MatrixCSR():
    def __init__(self):
        self.shape = None
        self.indptr = None
        self.indices = None
        self.barcodes = None
        self.genes = None
        self.data = None
        self.coo = None
        self.gene_names = None
        self.trans_m=None

    def read_h5(self, input_file):
        h5 = h5py.File(input_file, 'r')
        self.data = h5['/GRCh38/data']
        self.indices = h5['/GRCh38/indices']
        self.indptr = h5['/GRCh38/indptr']
        self.genes = h5['/GRCh38/genes']
        self.barcodes = h5['/GRCh38/barcodes']
        self.shape = h5['/GRCh38/shape']
        self.gene_names = h5['/GRCh38/gene_names']
        csr_m = csr_matrix((self.data, self.indices, self.indptr), shape=(self.shape[1], self.shape[0]))
        self.trans_m=np.transpose(csr_m)

    def write2mtx(self, write_path):
        barcodes_path = write_path+'barcodes.tsv'
        file = open(barcodes_path, 'w')
        print(str(self.barcodes[0],'utf-8'))
        for i in range(self.shape[1]):
            file.write(str(self.barcodes[i],'utf-8'))
            file.write('\n')
        file.close()
        print("barcodes writed")

        genes_path=write_path+'genes.tsv'
        file = open(genes_path, 'w')
        print(self.gene_names[0])
        print(str(self.gene_names[0],encoding ="utf8"))
        for i in range(self.shape[0]):
            file.write(str(self.genes[i],'utf-8'))
            file.write('\t')
            file.write(str(self.gene_names[i],'utf-8'))
            file.write('\n')
        file.close()
        print("genes writed")
        print(self.shape[0])
        data_path = write_path + 'matrix.mtx'
        mmwrite(data_path,self.trans_m)
        print("data writed")

def convert_h52mtx(input_file,write_path):
    mat=MatrixCSR()
    mat.read_h5(input_file)
    mat.write2mtx(write_path)