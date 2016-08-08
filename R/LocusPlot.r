########### function for plotting individual allele distances, colored by ploidy ########
LocusPlot <- function(x){
  par(mfrow=c(1,1))
  locus <- x
  if (locus %in% allMat[,1]==TRUE){
    locusMat <- allMat[allMat[,1]==locus,]
    plot(1:length(locusMat[,1]),rep(0,length(locusMat[,1])), ylim=c(-.5,.5),cex=.01, main=toString(locus), ylab="Allele Distance", xlab = "Individual")
    for (i in 1:length(locusMat[,1])){
      points(i,locusMat[i,4], col = if(locusMat[i,3]==0){"green"} else if(locusMat[i,3]==2){"blue"} 						else{"red"},cex=1)
    }
    plot(1:length(locusMat[,1]),rep(0,length(locusMat[,1])), ylim=c(-.1,.1),cex=.01,main=toString(locus), ylab="Allele Distance", xlab = "Individual")
    for (i in 1:length(locusMat[,1])){
      points(i,locusMat[i,4], col = if(locusMat[i,3]==0){"green"} else if(locusMat[i,3]==2){"blue"} 						else{"red"},cex=1)
    }
  }
  
}