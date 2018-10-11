
	library(biomaRt)
	ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
	data_truseq<-read.csv(file="E:\\cell data\\GSE52529\\GSE52529_fpkm_matrix.csv")
	v1=data_truseq[,1]
	hgnc_swissprot <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),filters = 'ensembl_gene_id', values = v1, mart = ensembl)
	write.csv(hgnc_swissprot,file="E:\\cell data\\GSE52529\\GSE52529_fpkm_mapping.csv")
