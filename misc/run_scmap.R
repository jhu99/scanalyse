#!/usr/bin/env Rscript
library('SingleCellExperiment')
library('scmap')

path=sys.argv[1]
respath=sys.argv[2]

X <- as.matrix(readMM(paste0(path,"matrix.mtx")))
Xm <- scale(X,center = FALSE,scale = colSums(X)/1e6 )
fData <- read.table(paste0(path,"genes.tsv"),sep="\t",stringsAsFactors = FALSE, row.names = 1, col.names = c("gene symbol","ensembl"))
pData <- read.table(paste0(path,"barcodes.tsv"),sep="\t",stringsAsFactors = FALSE, row.names = 1, col.names = c("barcode","louvain"))

# Index scmap-cluster
# bm <- SingleCellExperiment(assays = list(counts=as(X,"matrix")), colData=pData, rowData=fData)
bm <- SingleCellExperiment(assays = list(normcounts=Xm), colData=pData, rowData=fData)
logcounts(bm) <- log2(normcounts(bm) + 1)
rowData(bm)$feature_symbol <- rownames(bm)
bm <- selectFeatures(bm, suppress_plot = FALSE)
bm <- indexCluster(bm,cluster_col = "louvain")
heatmap(as.matrix(metadata(bm)$scmap_cluster_index))

scmapCluster_results <- scmapCluster(
  projection = bm, 
  index_list = list(
    m3drop = metadata(bm)$scmap_cluster_index
  )
)
plot(
  getSankey(
    colData(bm)$louvain, 
    scmapCluster_results$scmap_cluster_labs[,'m3drop'],
    plot_height = 400
  )
)
# Index Cell scmap-cell
set.seed(1)
bm <- indexCell(bm)
names(metadata(bm)$scmap_cell_index)
scmapCell_results <- scmapCell(
  bm, 
  list(
    md3drop = metadata(bm)$scmap_cell_index
  )
)
scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results, 
  list(
    as.character(colData(bm)$louvain)
  )
)
plot(
  getSankey(
    colData(bm)$louvain, 
    scmapCell_clusters$scmap_cluster_labs[,"md3drop"],
    plot_height = 400
  )
)

