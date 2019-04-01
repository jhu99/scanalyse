library(garnett)
barcodes<-read.table("E:\\graduated\\data\\mtx_test02\\barcodes.tsv",header = F,sep="\t")
genes<-read.table("E:\\graduated\\data\\mtx_test02\\genes.tsv",header = F,sep="\t")
mat <- readMM("E:\\graduated\\data\\mtx_test02\\matrix.mtx")
row.names(mat) <- row.names(genes)
colnames(mat) <- row.names(barcodes)
pd <- new("AnnotatedDataFrame", data = barcodes)
fd <- new("AnnotatedDataFrame", data = genes)
pbmc_cds <- newCellDataSet(as(mat, "dgCMatrix"),
                             phenoData = pd,
                             featureData = fd)