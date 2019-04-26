library('SingleCellExperiment')
library('scmap')
library('Matrix')
library('readr')
args = commandArgs(trailingOnly=TRUE)

path=args[1]
respath=args[2]
print(c(path,respath))
X <- as.matrix(readMM(paste0(path,"matrix.mtx")))
Xm <- scale(X,center = FALSE,scale = colSums(X)/1e6 )
fData <- read.table(paste0(path,"genes.tsv"),sep="\t",stringsAsFactors = FALSE, row.names = 1, col.names = c("gene symbol","ensembl"))
pData <- read.table(paste0(path,"barcodes.tsv"),sep="\t",stringsAsFactors = FALSE, row.names = 1, col.names = c("barcode","louvain","test"))
pData["test"] <- sample(c(TRUE,FALSE),dim(pData)[1],replace = TRUE,prob = c(0.1,0.9))
# Index scmap-cluster
bm <- SingleCellExperiment(assays = list(normcounts=Xm), colData=pData, rowData=fData)
logcounts(bm) <- log2(normcounts(bm) + 1)
rowData(bm)$feature_symbol <- rownames(bm)
bmref <- bm[,colData(bm)$test==FALSE]
bmtest <- bm[,colData(bm)$test==TRUE]
pdf(file=paste0(respath,"feture_scmap_cluster.pdf"))
bmref <- selectFeatures(bmref, suppress_plot = FALSE)
dev.off()
bmref <- indexCluster(bmref,cluster_col = "louvain")
pdf(file=paste0(respath,"heatmap_scmap_cluster.pdf"))
heatmap(as.matrix(metadata(bmref)$scmap_cluster_index))
dev.off()

scmapCluster_results <- scmapCluster(
  projection = bmtest,
  index_list = list(
    m3drop = metadata(bmref)$scmap_cluster_index
  )
)
# plot(
#   getSankey(
#     colData(bm)$louvain,
#     scmapCluster_results$scmap_cluster_labs[,'m3drop']
#     plot_height = 400
#   )
# )
res <- data.frame(louvain=colData(bmtest)$louvain,
                  scmapcluster=scmapCluster_results$scmap_cluster_labs[,'m3drop'])


# Index Cell scmap-cell
set.seed(1)
bmref <- indexCell(bmref)
scmapCell_results <- scmapCell(
  bmtest,
  list(
    m3drop = metadata(bmref)$scmap_cell_index
  )
)
scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results,
  list(
    as.character(colData(bm)$louvain)
  )
)
# plot(
#   getSankey(
#     colData(bm)$louvain,
#     scmapCell_clusters$scmap_cluster_labs[,"m3drop"],
#     plot_height = 400
#   )
# )
res$scmapcell=scmapCell_clusters$scmap_cluster_labs[,"m3drop"]
write_csv(res,path = paste0(respath,"scmap_assignment.csv"))
