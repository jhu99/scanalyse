# R script for performance comparison

filenames_saverX=list.files("./result",recursive=T,pattern="rds",full.names=T)
filenames_myres=list.files("./result",recursive=T,pattern="mean",full.names=T)
filenames_raw=list.files("./data/geneFilterResult",pattern="csv",full.names=T)
filenames_log=list.files("./data", recursive=T, pattern="293.*2_[0-9]{1,2}_log.csv", full.names=T)
raw_x <- as.matrix(read.csv(filenames_raw,row.names=1))
maskInd <- rep(0L,2000*3000)

saverx_mse_samples<-rep(0,10)
myres_mse_samples<-rep(0,10)
for(i in 1:10)
{
	fsaverx <- filenames_saverX[i]
	fmyres <- filenames_myres[i]
	flog <- filenames_log[i]
	x1 <- readRDS(fsaverx)
	x2 <- as.matrix(read.csv(fmyres,row.names=1))
	ind <- as.matrix(read.csv(flog,header=F,col.names=paste('V',1:500,sep='')))
	dimInd <- dim(ind)
	rownum <- dimInd[1]
	colnum <- dimInd[2]
	try(if(rownum!=2885) stop("Invalid log files!"))
	
	locationInd=1
	for(m in 1:rownum)
	{
		cellid = m-1	
		for(n in 2:colnum){
			if(is.na(ind[m,n]))next
			geneid=ind[m,n]+1 
			location = cellid*colnum+geneid
			maskInd[locationInd] = location
			locationInd=locationInd+1
		}
	}
	
	try(if(any(raw_x[maskIndValid]==0)) stop("Invalid log files!"))

	maskIndValid = maskInd[maskInd>0]
	saverx_mse_samples[i] <- mean((x1[maskIndValid]-raw_x[maskIndValid])^2)
	myres_mse_samples[i] <- mean((x2[maskIndValid]-raw_x[maskIndValid])^2)
	
	maskInd[] <- 0
	
}































