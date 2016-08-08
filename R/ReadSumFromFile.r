ReadSumFromFile <- function (z){
  ######### Sample Input:
  ######### ReadSumFromFile("/Volumes/P0027_Gerhardt/P0027")
  allMat <- as.matrix(read.csv(paste(c(z, "_allMat.csv"), collapse = ""), header = FALSE))
  info <- read.csv((paste(c(z, "_RefTable.csv"), collapse = "")))
  assign("allMat",allMat,globalenv())	
  assign("info",info,globalenv())	
}