#' @export
DeterminePloidy <- function(meanThresh, sdThresh, inds){
  ########### Ind x Locus agreement with informative ######################################
  if(missing(inds)){} else {cutMat <- rankMat[1:inds,]}
  if(missing(meanThresh)){indThresh <- .75} else {indThresh <- meanThresh}
  if(missing(sdThresh)){sThresh <- .8} else {sThresh <- sdThresh}
  
  meanLocusDif <- (alleleCheckMat[,2]+alleleCheckMat[,1])/2
  individuals <- unique(allMat[,2])
  infTestMat <- matrix(nrow=length(individuals)+1, ncol=length(cutMat[,1])+1)
  infTestMat[,1] <- c(0,individuals)
  infTestMat[1,] <- c(0,cutMat[,1])
  ratioMat <- matrix (nrow=length(individuals), ncol = 3)
  totalsMat <- matrix(nrow=length(individuals),ncol=3)
  for (i in 1:length(individuals)){
    tempMeanSum <- 0
    
    rowRef <- which(grepl(paste(c("^I",individuals[i],"$"), collapse = ""), info$Individual))
    
    tSum <- 0
    fSum <- 0
    iCount = i+1
    indMat <- allMat[allMat[,2]==individuals[i],]
    for (j in 1:length(cutMat[,1])){
      jCount = j+1
      locusMean <- meanLocusDif[which(testLoci==cutMat[j,1])]
      tetSd <- sdOutsCheckMat[which(testLoci==cutMat[j,1]),2]
      indLocMean <- mean(indMat[indMat[,1]==cutMat[j,1],4])
      indSdofOuts <- mean(indMat[indMat[,1]==cutMat[j,1],5])
      tempMeanSum <- tempMeanSum + indLocMean
      if (indLocMean != "NaN"){
        if (indLocMean >= locusMean+-(locusMean*indThresh) & indSdofOuts >=tetSd+-(tetSd*sThresh)){
          infTestMat[iCount,jCount] <- 1
          tSum <- tSum +1
        }
        else{
          infTestMat[iCount,jCount] <- 0
          fSum <- fSum +1
        }
      }	
      
      else{}
    }
    
    ratioMat[i,] <- c(individuals[i],tSum,info$PrePloidy[rowRef])
    
    
  }
  
  plot(ratioMat[,3]+runif(length(ratioMat[,3]),-0.1,0.1),ratioMat[,2], xlab="Ploidy", ylab="Number of Agreeing Informative Loci")
  
  tetMin <- 1
  dipMax <-1
  tetMean <-1
  dipMean <-1
  
  tets <- ratioMat[ratioMat[,3]==4,]
  avgAgree <- mean(tets[,2])
  sdAgree <- sd(tets[,2])
  badTet <- tets[tets[,2]<(avgAgree-(1.96*sdAgree)),]
  
  dips <- ratioMat[ratioMat[,3]==2,]
  badDip <- dips[dips[,2]>(avgAgree-(1.96*sdAgree)),]
  
  unknowns <- ratioMat[ratioMat[,3]==0,]
  
  assign("badTet", badTet, globalenv())
  assign("badDip", badDip, globalenv())
  assign("unknowns", unknowns, globalenv())
  assign("infMat", ratioMat, globalenv())
}
