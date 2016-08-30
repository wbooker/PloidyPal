#' @export
CalcAlleleDiffs <- function(f){
  infoTable <- as.matrix(read.csv(f, header=TRUE))
  BEG1 <- as.numeric(infoTable[1,2])
  END1 <- as.numeric(infoTable[2,2])
  str1 <- toString(infoTable[4,2])
  
  for(j in BEG1:END1){
    if (file.exists(filePath) == TRUE){
    
    filePath <- paste(c(str1,"/I",j,"/I",j,"_allelesFromPost_4.txt"), collapse = "")	
    outFilePath <- paste(c(str1,"/I",j,"/I",j,"_allelesFromPost_4_diffs.txt"), collapse = "")	
    outTable <- matrix(nrow = 1, ncol = 8)	
    k <-scan(filePath, sep = ">", what = "complex")	
    for(i in seq(2,length(k),12)){
      textMat <- matrix(nrow = 6, ncol = 8)
      allele <- strsplit(k[i],".", fixed = TRUE)
      locus <- strsplit(allele[[1]][1],"L")[[1]][2]
      copy <- allele[[1]][2]
      alleleNumber <- allele[[1]][3]
      a1 <- k[i+1]
      a2 <- k[i+4]
      a3 <- k[i+7]
      a4 <- k[i+10]
      textMat[,c(1,4)] <- locus
      textMat[,c(2,5)] <- copy
      textMat[,3] <- c(1,1,1,2,2,3)
      textMat[,6] <- c(2,3,4,3,4,4)
      textMat[,7] <- c(adist(a1,a2),adist(a1,a3), adist(a1,a4), adist(a2,a3), adist(a2,a4), adist(a3,a4))
      textMat[,8] <- nchar(a1)
      outTable <- rbind(outTable,textMat)
    }
    outTable <- outTable[2:length(outTable[,1]),]
    write.table(outTable,file = outFilePath, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    print(j)
    
    }
    else {
     print(paste(c("I",j," FILE NOT FOUND"), collapse = ""), quote = FALSE)
    }
    
  }
  assign("infoTable",infoTable,globalenv())
}
