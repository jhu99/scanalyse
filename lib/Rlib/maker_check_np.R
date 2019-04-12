library(garnett)
path="./data/hca/ica_bone_marrow_filtered_gene_bc_matrices_mex2/hg19/"
X <- readMM(paste0(path,"matrix.mtx"))
fData <- read.table(paste0(path,"genes.tsv"),sep="\t",stringsAsFactors = FALSE, row.names = 1, col.names = c("gene symbol","ensembl","counts"))
pData <- read.table(paste0(path,"barcodes.tsv"),stringsAsFactors = FALSE, row.names = 1)
rownames(X) <- row.names(fData)
colnames(X) <- row.names(pData)
fData$gene_short_name = row.names(fData)

# create a new CDS object
pd <- new("AnnotatedDataFrame", data = pData)
fd <- new("AnnotatedDataFrame", data = fData)
ica_bm_cds <- newCellDataSet(as(X, "dgCMatrix"),
													 phenoData = pd,
													 featureData = fd)

# generate size factors for normalization later
ica_bm_cds <- estimateSizeFactors(ica_bm_cds)

# check marker file
library(org.Hs.eg.db)
marker_file_path="./data/hca/ica_bones_cell_type_markers.txt"
marker_check <- check_markers(ica_bm_cds, marker_file_path,
															db=org.Hs.eg.db,
															cds_gene_id_type = "SYMBOL",
															marker_file_gene_id_type = "SYMBOL")
save(ica_bm_cds,file = paste0(path,"garnett_cds.Rdata"))
write.table(marker_check,file = paste0(path,"ica_bone_marker_check.txt"))
df=marker_check[which(marker_check$ambiguity>0.5),]
write.table(df,file = paste0(path,"ica_bone_marker_check_remove.tsv"))

for(i in 0:10)
{
	print(i)
	t = i*30
	if(i<10){
		t=t+1:30
	}else{
		t=t+1:33
	}
	print(t)
	marker_check_sub = marker_check[t,]
	browser()
	# pdf(paste0(path,"marker-",i,".pdf"))
	# plot_markers(marker_check_sub)
	# dev.off()
}
# pdf(paste0(path,"marker-",i,".pdf"))
# plot_markers(marker_check_sub)
# dev.off()


# revise marker files after plot_markers
##################################
# train classifier based on the revised marker file in the next step.


