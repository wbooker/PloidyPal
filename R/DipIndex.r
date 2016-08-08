#' @export
DipIndex <- function(){
  preInums <- unique(alleleMat[,1])
  tetNums <- array()
  dipNums <- array()	
  for(i in 1:length(preInums)){
    rowRef <- which(grepl(paste(c("^I",preInums[i],"$"), collapse = ""), info$Individual))
    if (info$PrePloidy[rowRef]==4){
      tetNums <- c(tetNums, preInums[i])
    }
    else if (info$PrePloidy[rowRef]==2){
      dipNums <- c(dipNums, preInums[i])
    }
    else{}
  }
  tetNums <- tetNums[2:length(tetNums)]
  dipNums <- dipNums[2:length(dipNums)]
  alleleAgreeMat <- matrix(nrow=length(tetNums)+6, ncol=numLoci+1)
  alleleAgreeMat[1,] <- c("Individual",testLoci)
  alleleAgreeMat[,1] <- c("Individual", tetNums, "Individuals with 4 alleles", "Percent Individuals Retained", "Retention against reference", "Avg", "Y/N Non-Dip")
  nonDipLoci <- array()
  
  refAgreeMat <- matrix(nrow=length(dipNums)+4, ncol=numLoci+1)
  refAgreeMat[1,] <- c("Individual",testLoci)
  refAgreeMat[,1] <- c("Individual", dipNums, "Individuals with 4 alleles", "Percent Individuals Retained", "Avg")
  library(Matrix)
  for (i in 1:length(tetNums)){
    tempAllele <- alleleMat[alleleMat[,1]==tetNums[i],]
    for (j in 1:numLoci){
      allAgree <- 0
      if (testLoci[j] %in% tempAllele[,2] == TRUE){
        tempLoc <- matrix(ncol = length(tempAllele[1,]))
        tempLoc <- rbind(tempLoc, tempAllele[tempAllele[,2]==testLoci[j],])
        for(k in 1:(length(tempLoc[,1])-1)){
          #if(nnzero(tempLoc[k+1,3:8]) == 6){
          #	allAgree <- 1
          #}
          #else{}
          allAgree <- nnzero(tempLoc[k+1,3:8])
          if(allAgree == 0){
            
          }
          else{allAgree <- allAgree - 2}
        }
      }
      else{}
      alleleAgreeMat[i+1,j+1] <- allAgree
    }
    
    
  }
  
  
  for (i in 1:length(dipNums)){
    tempAllele <- alleleMat[alleleMat[,1]==dipNums[i],]
    for (j in 1:numLoci){
      allAgree <- 0
      if (testLoci[j] %in% tempAllele[,2] == TRUE){
        tempLoc <- matrix(ncol = length(tempAllele[1,]))
        tempLoc <- rbind(tempLoc, tempAllele[tempAllele[,2]==testLoci[j],])
        for(k in 1:(length(tempLoc[,1])-1)){
          if(nnzero(tempLoc[k+1,3:8]) >= 4){
            allAgree <- 2
          }
          else{allAgree <- nnzero(tempLoc[k+1,3:8])
          if(allAgree == 0){}
          else{allAgree <- allAgree - 2}
          }
        }
      }
      else{}
      refAgreeMat[i+1,j+1] <- allAgree
    }
    
    
  }
  
  for (i in 1:numLoci){
    
    tempCol <- refAgreeMat[2:(length(refAgreeMat[,1])-3), i+1] 
    refAgreeMat[length(refAgreeMat[,1])-2,i+1] <- sum(as.numeric(tempCol==2))
    refDiv <- sum(as.numeric(tempCol==2))/length(tempCol)
    refAgreeMat[length(refAgreeMat[,1])-1,i+1] <- refDiv
    refAgreeMat[length(refAgreeMat[,1]),i+1] <- mean(as.numeric(tempCol))
    projAlleleNum <- mean(as.numeric(tempCol))*2 ##### Really div by 2 and times 4 but you know
    
    tempCol <- alleleAgreeMat[2:(length(alleleAgreeMat[,1])-5), i+1] 
    alleleAgreeMat[length(alleleAgreeMat[,1])-4,i+1] <- sum(as.numeric(tempCol==4))
    tetDiv <- sum(as.numeric(tempCol==4))/length(tempCol)
    alleleAgreeMat[length(alleleAgreeMat[,1])-3,i+1] <- tetDiv
    retPercent <- if(tetDiv==0 & refDiv==0) {"Conserved Across"} else{tetDiv/refDiv}
    alleleAgreeMat[length(alleleAgreeMat[,1])-2,i+1] <- retPercent
    meanAlleleComps <- mean(as.numeric(tempCol))
    alleleAgreeMat[length(alleleAgreeMat[,1])-1,i+1] <- meanAlleleComps
    if(meanAlleleComps <= projAlleleNum+-0.1){
      alleleAgreeMat[length(alleleAgreeMat[,1]),i+1]	<- 0
    }
    else{
      alleleAgreeMat[length(alleleAgreeMat[,1]),i+1]	<- 1
      nonDipLoci <- c(nonDipLoci, alleleAgreeMat[1,i+1])
      
    }
  }
  nonDipLoci <- as.numeric(nonDipLoci[2:length(nonDipLoci)])
  assign("nonDipLoci", nonDipLoci, globalenv())
  write.table(nonDipLoci, file = paste(c(str1,"/PloidyRFiles/nonDipLoci.txt"), collapse = ""), row.names=FALSE, 		col.names=FALSE, sep=",")
  assign("alleleAgreeMat", alleleAgreeMat, globalenv())
  write.table(alleleAgreeMat, file = paste(c(str1,"/PloidyRFiles/DiploidizationIndex.csv"), collapse = ""), row.names=FALSE, 		col.names=FALSE, sep=",")	
  
  assign("refAgreeMat", refAgreeMat, globalenv())
  write.table(refAgreeMat, file = paste(c(str1,"/PloidyRFiles/RefDipIndex.csv"), collapse = ""), row.names=FALSE, 			   col.names=FALSE, sep=",")		
}