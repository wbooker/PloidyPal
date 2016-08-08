\name{FindInformativeLoci}
\alias{FindInformativeLoci}
\alias{CutInfs}
\alias{UncutInfs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find Informative Loci}
\description{
ses the info from CheckPloidy() to find which loci are informative based on some previously known individuals.}
\usage{
FindInformativeLoci(outlierThresh, numOutliers, informativeThresh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outlierThresh}{
Threshold for what allelic distance constitues an outlier
}
  \item{numOutliers}{
number of outliers allowed per locus}
  \item{informativeThresh}{
Distance between average tetraploid vs diploid that makes a loci count as informative}
}
\details{
CutInfs takes a number and shortens the informative loci to that number.
UncutInfs undoes that process
}
\value{
cutMat is a sorted matrix of informative loci and their distance differences.
}

\author{
William Booker}






\keyword{ PloidyPal }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
