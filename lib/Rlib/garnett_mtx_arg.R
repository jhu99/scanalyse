Args <- commandArgs()
library(garnett)
barcodes<-read.table(paste(Argvs[1],"barcodes.tsv",sep=""),header = F,sep="\t")
genes<-read.table(paste(Argvs[1],"genes.tsv",sep=""),header = F,sep="\t")
mat <- readMM(paste(Argvs[1],"matrix.mtx",sep=""))
row.names(mat) <- row.names(genes)
colnames(mat) <- row.names(barcodes)
pd <- new("AnnotatedDataFrame", data = barcodes)
fd <- new("AnnotatedDataFrame", data = genes)
pbmc_cds <- newCellDataSet(as(mat, "dgCMatrix"),
                             phenoData = pd,
                             featureData = fd)