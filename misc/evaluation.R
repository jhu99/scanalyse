# R script for performance comparison

filenames_saverX=list.files("./result",recursive=T,pattern="rds",full.names=T)
filenames_myres=list.files("./result",recursive=T,pattern="mean",full.names=T)
filenames_raw=list.files("./data/geneFilterResult",pattern="csv",full.names=T)
filenames_log=list.files("./data", recursive=T, pattern="293.*2_[0-9]{1,2}_log.csv", full.names=T)

raw_x <- as.matrix(read.csv(filenames_raw,row.names=1))
maskInd <- rep(0L,2000*3000)
msev<-c()
corv<-c()
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
	genenum <- nrow(x1)
	try(if(rownum!=2885) stop("Invalid log files!"))
	
	locationInd=1
	maskInd[] <- 0
	for(m in 1:rownum)
	{
		cellid = m-1	
		for(n in 2:colnum){
			if(is.na(ind[m,n])){
				next
			}
			geneid=ind[m,n]+1 
			location = cellid*genenum+geneid
			maskInd[locationInd] = location
			locationInd=locationInd+1
		}
	}
	
	maskIndValid = maskInd[maskInd>0]
	try(if(any(raw_x[maskIndValid]==0)) stop("Invalid log files!"))

	msev <- c(msev,mean((x1[maskIndValid]-raw_x[maskIndValid])^2))
	msev <- c(msev,mean((x2[maskIndValid]-raw_x[maskIndValid])^2))
	corv <- c(corv,cor.test(x1[maskIndValid],raw_x[maskIndValid])$estimate)
	corv <- c(corv,cor.test(x2[maskIndValid],raw_x[maskIndValid])$estimate)
}
data <- data.frame(mse=msev,cor=corv,method=rep(c("saver-X","myres"),times=10))
save(data, file="./result/cor_mse_comparison.Rdata")
pdf("./result/figures/f1.pdf")
boxplot(mse~method,data=data)
boxplot(cor~method,data=data)
dev.off()






