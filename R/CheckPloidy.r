#' @export
CheckPloidy <- function(f){
  infoTable <- as.matrix(read.csv(f, header=TRUE))
  assign("infoTable",infoTable,globalenv())
  SummaryFromFiles()
}