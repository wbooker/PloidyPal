CutInfs <- function(len){
  cutMat <- rankMat[1:len,]
  assign("cutMat", cutMat, globalenv())
}