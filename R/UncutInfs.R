#' @export
UncutInfs <- function(){
  cutMat <- rankMat
  assign("cutMat", cutMat, globalenv())
}