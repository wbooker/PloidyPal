\name{DeterminePloidy}
\alias{DeterminePloidy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine Ploidy}
\description{
Graphs individuals for agreeing informative sites by ploidy 
}
\usage{
DeterminePloidy(meanThresh, sdThresh, inds)
}
\details{
	meanThresh is your plus minus threshold (between 0 and 1) multiplier to mean allelic distance for comparison of informative loci. Default is 0.75\cr
	sdThresh is your plus minus threshold (between 0 and 1) multiplier to mean tetraploid standard deviation for comparison of informative loci. Default is 0.8\cr
	inds is number of informative loci to be considered ranked in descending order. Default is all informative loci 
}
\value{
badDip is a matrix of previously assumed diploids that are potentially mislabeled and their number of agreeing sites

badTet is similar for tetraploids

unknowns is a similar matrix for unknown individuals

infMat is a matrix of all individuals with ploidy level and number of informative sites
}

\author{
William Booker, FSU}
\examples{
DeterminePloidy(meanThresh = 0.6, sdThresh = 0.7, inds = 12)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PloidyPal }
