FindInformativeLoci <- function(outlierThresh, numOutliers, informativeThresh) {
  ######## Function to plot and calculate informative loci based on outlier threshold and difference threshold #####
  
  ######## Defaults ##########
  #  outlierThresh <- 0.1  ### level of divergence determining that individual as outlier
  #  numOutliers <- 4      ### number of outliers allowed per locus
  #  informativeThresh <- 0.0055  
  
  if(missing(outlierThresh)){outThresh <- 0.1} else {outThresh <- outlierThresh}
  if(missing(numOutliers)){numOutThresh <- 4} else {numOutThresh <- numOutliers}
  if(missing(informativeThresh)){infoThresh <- 0.0055} else {infoThresh <- informativeThresh}
  alleleCheckMat <- matrix(nrow=numLoci,ncol=2)
  sdOutsCheckMat <- matrix(nrow=numLoci, ncol=2)
  sdDifs <- matrix(nrow=numLoci,ncol=1)  
  outMat <- matrix(nrow=numLoci,ncol=2)
  dipTemp <- c()
  tetTemp <- c()
  for (i in 1:numLoci){
    locusMat <- allMat[allMat[,1]==testLoci[i],]
    sdDifs[i] <- sd(locusMat[,4])
    dipTemp <- locusMat[locusMat[,3]==2,4]
    dipSdTemp <- locusMat[locusMat[,3]==2,5]
    outMat[i,1] <- sum(dipTemp > outThresh)
    dipTemp <- dipTemp[dipTemp <= outThresh]	 
    tetTemp <- locusMat[locusMat[,3]==4,4]
    tetSdTemp <- locusMat[locusMat[,3]==4,5]
    outMat[i,2] <- sum(tetTemp > outThresh)
    tetTemp <- tetTemp[tetTemp <= outThresh]
    alleleCheckMat[i,1] <- mean(dipTemp)
    alleleCheckMat[i,2] <- mean(tetTemp)
    sdOutsCheckMat[i,1] <- mean(dipSdTemp)
    sdOutsCheckMat[i,2] <- mean(tetSdTemp)
  }
  alleleDifMat <- alleleCheckMat[,2] - alleleCheckMat[,1]
  iLoci <- c()	
  InfsPDFPath = paste (c(str1,"/PloidyRFiles/InformativeLocigraph.pdf"), collapse = "")
  pdf(file = InfsPDFPath )
  plot(1:testLoci[length(testLoci)],rep(0,testLoci[length(testLoci)]), ylim=c(-.025,.03),cex=.01, xlab = "Locus", ylab = "Level of Divergence")
  for (i in 1:numLoci){
    if(sum(outMat[i,]) <= numOutThresh){
      if(alleleDifMat[i] != "NaN"){
        if(alleleDifMat[i] >= infoThresh){
          color <- "blue"
          iLoci<- rbind(iLoci, as.numeric(testLoci[i]))
          text(testLoci[i]+9,alleleDifMat[i],testLoci[i])
        }
        else{
          color <- "black"
        }
      }
      else{}
      
    } 
    else{
      color = "red"
    }
    points(testLoci[i], alleleDifMat[i], col = color)
  }
  
  ##################################   RANKING   ##############################################################		
  rankMat	<- matrix(nrow=length(iLoci),ncol=2)
  for (i in 1:length(iLoci)){
    rankMat[i,1] <- iLoci[i,]
    rankMat[i,2] <- alleleDifMat[which(testLoci==iLoci[i,])]	
  }
  rankMat <- rankMat[sort.list(rankMat[,2],decreasing=TRUE),]
  cutMat <- rankMat	
  assign("iLoci", iLoci, globalenv())
  assign("alleleCheckMat",alleleCheckMat,globalenv())	
  assign("alleleDifMat",alleleDifMat,globalenv())	
  assign("cutMat", cutMat, globalenv())
  assign("rankMat", cutMat, globalenv())
  assign("sdOutsCheckMat", sdOutsCheckMat, globalenv())
  write.table(sdOutsCheckMat, file = paste(c(str1,"/PloidyRFiles/sdOutsCheckMat.csv"), collapse = ""), row.names=FALSE, 		col.names=FALSE, sep=",")
  dev.off()
  cat(c("Writing:" , paste(c(str1,"/PloidyRFiles/InformativeLociGraph.pdf"), collapse = ""), "\n")) 
}