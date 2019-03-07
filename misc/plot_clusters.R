# Plot clusters
library(Rtsne)
ifile <- "./result/ica_all/latent.csv"
proj <- read.csv(ifile, row.names=1, header=FALSE)
proj <- proj[!duplicated(proj),]
tsne <- Rtsne(as.matrix(proj), checkduplicates=FALSE, pca= FALSE, dims=2, theta =0.5, verbose= TRUE, num_threads=12)
write.csv(tsne$Y, file="./result/ica_all/projection_tsne.csv")
