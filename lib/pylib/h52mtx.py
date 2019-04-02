from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix
from anndata import h5py

class MatrixCSR():
    def __init__(self):
        self.shape=None
        self.indptr=None
        self.indices=None
        self.barcodes=None
        self.genes=None
        self.data=None
        self.coo=None
        self.gene_names=None

    def read_h5(self,input_file):
        h5 = h5py.File(input_file, 'r')
        self.data = h5['/GRCh38/data']
        self.indices = h5['/GRCh38/indices']
        self.indptr = h5['/GRCh38/indptr']
        self.genes = h5['/GRCh38/genes']
        self.barcodes = h5['/GRCh38/barcodes']
        self.shape = h5['/GRCh38/shape']
        self.gene_names=h5['/GRCh38/gene_names']
        csr_m = csr_matrix((self.data, self.indices, self.indptr), shape=(self.shape[1],self.shape[0]))
        self.coo = csr_m.tocoo()
        print(self.coo.col[0])
        print(self.barcodes[0])
        print(self.genes[0])

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
        data_path=write_path+'matrix.mtx'
        file = open(data_path, 'w')
        file.write('%%MatrixMarket matrix coordinate real general\n')
        file.write('%\n')
        file.write(str(self.shape[0]))
        file.write('\t')
        file.write(str(self.shape[1]))
        file.write('\t')
        file.write(str(len(self.data)))
        print(len(self.data))
        print(len(self.coo.col))
        file.write(' \n')
        for i in range(len(self.coo.col)):
            file.write(str(self.coo.col[i]+1))
            file.write(' ')
            file.write(str(self.coo.row[i]+1))
            file.write(' ')
            file.write(str(self.coo.data[i]))
            file.write(' ')
            file.write('\n')
        file.close()
        print("data writed")

def convert_h52mtx(input_file,write_path):
    mat=MatrixCSR()
    mat.read_h5(input_file)
    mat.write2mtx(write_path)