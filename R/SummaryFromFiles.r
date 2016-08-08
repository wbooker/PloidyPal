#' @export
SummaryFromFiles <- function(){
  BEG1 <- as.numeric(infoTable[1,2])
  END1 <- as.numeric(infoTable[2,2])
  str1 <- toString(infoTable[4,2])
  assign("str1", str1, globalenv())
  dir.create(paste(c(str1, "/PloidyRFiles"), collapse= ""), showWarnings = FALSE)
  str1.1 <- "/I"
  str2 <- "/I"
  str3 <- "_allelesFromPost_4_diffs.txt"
  infoFile <- toString(infoTable[6,2])
  info <- read.csv(infoFile)
  Notes <- table(info$PrePloidy)
  LocusRef <- scan(infoTable[11,2])
  testLoci <- sort(LocusRef)
  numLoci <- length(testLoci)
  refCheck <- as.numeric(infoTable[7,2])
  alleleMat <- matrix(ncol=8)
  if(refCheck == 1){
    allMat <- as.matrix(read.csv(infoTable[8,2], header = FALSE) )
    allDistancesMat <- as.matrix(read.csv(infoTable[9,2], header = FALSE))
    prevInfo <- read.csv(infoTable[10,2], header = TRUE)
    info <- rbind(info, prevInfo)
  } else{
    allMat <- matrix(ncol=5)
    allDistancesMat <- matrix(data=c("Locus",testLoci),nrow=numLoci+1,ncol=1)
  }
  assign("info",info,globalenv())
  ########## This part summarizes all allele data into allMat ################
  for (i in BEG1:END1){
    filePath <- paste(c(str1,str1.1,i,str2,i,str3), collapse = "")
    iNum <- paste(c("I",i), collapse = "")
    if (file.exists(filePath) == TRUE){
      dat <- read.table(filePath)	
      matLen <- length(dat[,7])/6
      valuesMat <- matrix(ncol=5)
      rowRef <- which(grepl(paste(c("^",iNum,"$"), collapse = ""), info$Individual))
      alleleCounter <- 1
      for (j in 1:matLen){
        if(dat[(alleleCounter),1] <= testLoci[length(testLoci)] & dat[(alleleCounter),1] %in% testLoci){
          alleles <- sort(c((dat[(alleleCounter),7]/dat[alleleCounter,8]), (dat[(alleleCounter+1),7]/dat[alleleCounter,8]), 					(dat[(alleleCounter+2),7]/dat[alleleCounter,8]), (dat[(alleleCounter+3),7]/dat[alleleCounter,8]), (dat[(alleleCounter+4),				7]/dat[alleleCounter,8]), (dat[(alleleCounter+5),7]/dat[alleleCounter,8])))
          allAlleleDists <- alleles		
          sameAlleles <- alleles[3:6]
          alleles <- alleles[1:2]
          
          valuesMat <- rbind(valuesMat, c(dat[alleleCounter,1],i,info$PrePloidy[rowRef], mean(alleles), sd(sameAlleles)))
          alleleMat <- rbind(alleleMat,c(i,dat[alleleCounter,1],allAlleleDists[1],allAlleleDists[2], 										allAlleleDists[3],allAlleleDists[4], allAlleleDists[5], allAlleleDists[6]))
        }
        else{}
        alleleCounter <- alleleCounter+6
        
      }
      valuesMat <- valuesMat[2:length(valuesMat[,1]),]
      allMat <- rbind(allMat, valuesMat)
      ## tempVal <- mean(filtMat)
      indDist <- matrix(nrow=numLoci+1, ncol = 1)
      indDist[1,] = iNum
      for (k in 1:numLoci){
        M <- valuesMat[valuesMat[,1] == testLoci[k], 4]
        indDist[k+1,] <-mean(M)
      }	
      allDistancesMat <- cbind(allDistancesMat,indDist)
      print(paste(c(iNum, " COMPLETED"), collapse = ""), quote = FALSE)
    }
    
    else {
      print(paste(c(iNum," FILE NOT FOUND"), collapse = ""), quote = FALSE)
    }
  }
  
  
  
  
  
  
  write.table(allDistancesMat, file = paste(c(str1,"/PloidyRFiles/allDistancesMat.csv"), collapse=""), row.names=FALSE, col.names=FALSE, sep=",")
  cat(c("Writing:" ,paste(c(str1,"/PloidyRFiles/allDistancesMat.csv"), collapse = "") ,"\n"))
  
  allMat <- allMat[2:length(allMat[,1]),]  ####Important!
  
  write.table(allMat, file = paste(c(str1,"/PloidyRFiles/allMat.csv"), collapse=""), row.names=FALSE, col.names=FALSE, sep=",")
  cat(c("Writing:" , paste(c(str1,"/PloidyRFiles/allMat.csv"), collapse = ""), "\n"))
  alleleMat <- alleleMat[2:length(alleleMat[,1]),] 
  write.table(alleleMat, file = paste(c(str1,"/PloidyRFiles/alleleMat.csv"), collapse=""), row.names=FALSE, col.names=FALSE, sep=",")
  cat(c("Writing:" , paste(c(str1,"/PloidyRFiles/alleleMat.csv"), collapse = ""), "\n"))
  assign("allMat",allMat,globalenv())
  assign("alleleMat", alleleMat, globalenv())
  assign("testLoci",testLoci,globalenv())
  assign("numLoci", numLoci, globalenv())
  avgDist <- matrix(nrow=numLoci, ncol = 1)
  devDist <- matrix(nrow=numLoci, ncol = 1)
  for (k in 1:numLoci){
    M <- allMat[allMat[,1] == testLoci[k], 3]
    avgDist[k,] <- mean(M)
    devDist[k,] <- sd(M)
  }
}