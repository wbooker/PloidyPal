MakeGraphs <- function(){
  somePDFPath = paste (c(str1,"/PloidyRFiles/Allelegraphs.pdf"), collapse = "")
  pdf(file=somePDFPath)  
  for (i in 1:numLoci)   
  {   
    LocusPlot(testLoci[i])  
  } 
  dev.off()
  cat(c("Writing:" , paste(c(str1,"/PloidyRFiles/Allelegraphs.pdf"), collapse = ""), "\n")) 
}