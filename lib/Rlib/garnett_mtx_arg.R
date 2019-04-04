Args <- commandArgs()
library(garnett)
barcodes<-read.table(paste(Argvs[1],"barcodes.tsv",sep=""),header = F,sep="\t",row.names=1)
genes<-read.table(paste(Argvs[1],"genes.tsv",sep=""),header = F,sep="\t",row.names=1)
mat <- readMM(paste(Argvs[1],"matrix.mtx",sep=""))
row.names(mat) <- row.names(genes)
colnames(mat) <- row.names(barcodes)
pd <- new("AnnotatedDataFrame", data = barcodes)
fd <- new("AnnotatedDataFrame", data = genes)
pbmc_cds <- newCellDataSet(as(mat, "dgCMatrix"),
                             phenoData = pd,
                             featureData = fd)
pbmc_cds <- estimateSizeFactors(pbmc_cds)
library(org.Hs.eg.db)
marker_file_path <-"E:\\graduated\\data\\bones.txt"

marker_check <- check_markers(pbmc_cds, marker_file_path,
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")
head(pData(pbmc_cds))
pbmc_classifier <- train_cell_classifier(cds = pbmc_cds,
                                         marker_file = marker_file_path,
                                         db=org.Hs.eg.db,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL")
feature_genes <- get_feature_genes(pbmc_classifier,
                                   node = "root",
                                   db = org.Hs.eg.db)
head(feature_genes)
