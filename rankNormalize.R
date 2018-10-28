Rank.Normalize <- function(m, scale = NULL, threshold = NULL, ceiling = NULL, shift = NULL) {
   if (is.null(scale)) { scale = length(m[,1]) }
   if (!is.null(threshold)) { m[m < threshold] <- threshold }
   if (!is.null(ceiling)) { m[m > ceiling] <- ceiling }
   if (!is.null(shift)) { m <- m + shift }

   # column rank normalization from 1/N to 1 with "average" used to break ties.
   cols <- length(m[1,])
   for (j in 1:cols) {
      m[,j] <- rank(m[,j], ties.method = "average")
   }
   m <- m*(scale/length(m[,1]))
   return(m)
}
m<-read.csv(file="E:\\testcell\\GSE52529_truseq_fpkm_matrix.csv")
m<-Rank.Normalize(m,NULL,NULL,NULL,NULL)
write.csv(m,file="E:\\testresult03.csv")
