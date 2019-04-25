library(SingleCellExperiment)
library(scmap)

Argvs<- commandArgs()
print(Argvs[6])
print(Argvs[7])
#input
test<-read.csv(Argvs[6],header = T,sep=",",row.names=1)
types<-read.table(Argvs[7],header = T,sep=',')
sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(test)), colData = types)
logcounts(sce) <- log2(normcounts(sce) + 1)
rowData(sce)$feature_symbol <- rownames(sce)
isSpike(sce, "ERCC") <- grepl("^ERCC-", rownames(sce))
sce <- sce[!duplicated(rownames(sce)), ]

#Feature selection
sce <- selectFeatures(sce, suppress_plot = FALSE)
table(rowData(sce)$scmap_features)

#scmap-cluster
sce <- indexCluster(sce)
head(metadata(sce)$scmap_cluster_index)
heatmap(as.matrix(metadata(sce)$scmap_cluster_index))

scmapCluster_results <- scmapCluster(
  projection = sce, 
  index_list = list(
    yan = metadata(sce)$scmap_cluster_index
  )
)
print("finish cluster")
#Results
plot(
  getSankey(
    colData(sce)$cell_type1, 
    scmapCluster_results$scmap_cluster_labs[,'yan'],
    plot_height = 400
  )
)
print("finish plot")
