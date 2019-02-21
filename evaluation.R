# R script for performance comparison
library(readRDS)

filenames_saverX=list.files("./result",recursive=T,pattern="rds",full.names=T)
filenames_myres=list.files("./result",recursive=T,pattern="mean",full.names=T)
filenames_raw=list.files("./data/geneFilterResult",pattern="csv",full.names=T)
filenames_log=list.files("./data", recursive=T, pattern="293.*2_[0-9]{1,2}_log.csv", full.names=T)
raw_x <- read.csv(filenames_raw,row.names=1)
maskInd <- rep(0L,2000*3000)

for(i in 1:10)
{
	
	fsaverx <- filenames_saverX[i]
	fmyres <- filenames_myres[i]
	flog <- filenames_log[i]
	x1 <- readRDS(fsaverx)
	x2 <- read.csv(fmyres,row.names=1)
	ind <- read.csv(flog,header=F,col.names=paste('V',1:2000,sep=''))
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
	
	maskIndValid = maskInd[maskInd>0]
	mse1 <- mse(x1[maskIndValid]-rawx[maskIndValid])
	mse2 <- mse(x2[maskIndValid]-rawx[maskIndValid])
	
	maskInd[] <- 0
	
}






























