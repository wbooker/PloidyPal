\name{SummaryFromFiles}
\title{
Get Summary of Allele Info from Files
}
\description{
Gathers and summarizes allelic distance info from all _allelesFromPost_4_diffs.txt files and summarizes them into a few tables. One is allDistancesMat, which is a locus by individual matrix with the mean of the two smallest allelic distances per individual per locus. \cr \cr

Second is allMat which is a line by line summary of the same information \cr \cr

Third is alleleMat, which is a line by line summary of the other 4 allelic distances. \cr \cr

---All of these files are stored globally for use as well as written to .csv files in a created 
/PloidyRFiles directory.
}

\note{
	SummaryFromFiles(file)
	The filepath is to a paramfile .csv find a template and example at https://github.com/wbooker/PloidyPal
}

\author{
William Booker, FSU
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PloidyPal}
