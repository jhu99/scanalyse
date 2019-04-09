Argvs<- commandArgs()
library(garnett)
print(paste(Argvs[6],"barcodes.tsv",sep=""))
barcodes<-read.table(paste(Argvs[6],"barcodes.tsv",sep=""),header = F,sep="\t",row.names=1)
genes<-read.table(paste(Argvs[6],"genes.tsv",sep=""),header = F,sep="\t",row.names=1)
mat <- readMM(paste(Argvs[6],"matrix.mtx",sep=""))
row.names(mat) <- row.names(genes)
colnames(mat) <- row.names(barcodes)
pd <- new("AnnotatedDataFrame", data = barcodes)
fd <- new("AnnotatedDataFrame", data = genes)
pbmc_cds <- newCellDataSet(as(mat, "dgCMatrix"),
                             phenoData = pd,
                             featureData = fd)
pbmc_cds <- estimateSizeFactors(pbmc_cds)
library(org.Hs.eg.db)
marker_file_path <-Argvs[7]

marker_check <- check_markers(pbmc_cds, marker_file_path,
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "ENSEMBL",
                              marker_file_gene_id_type = "SYMBOL")
head(pData(pbmc_cds))
pbmc_classifier <- train_cell_classifier(cds = pbmc_cds,
                                         marker_file = marker_file_path,
                                         db=org.Hs.eg.db,
                                         cds_gene_id_type = "ENSEMBL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL")
feature_genes <- get_feature_genes(pbmc_classifier,
                                   node = "root",
                                   db = org.Hs.eg.db)
head(feature_genes)
pbmc_cds <- classify_cells(pbmc_cds, pbmc_classifier,
                           db = org.Hs.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "ENSEMBL")
head(pData(pbmc_cds))
type<-pData(pbmc_cds)$cell_type
write.table(type,Argvs[8],row.names = FALSE,col.names = FALSE)


