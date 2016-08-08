\name{PloidyPal-package}
\alias{PloidyPal-package}
\alias{PloidyPal}
\docType{package}
\title{
\packageTitle{PloidyPal}
Ploidy Pal
}
\description{
A package for running analyses of polyploid data. 
\packageDescription{PloidyPal}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{PloidyPal}
\packageIndices{PloidyPal}
CMPloidyPal completes a number of functions that can assist in and polyploid genetic analysis. As of 08/05/2016, only diploid and tetraploid datasets are functional. \cr \cr



CalcAlleleDiffs() takes files with 4 alleles per locus and calculates the differences between alleles for all alleles at all loci for each individual. These data are then stored in text files that is used for many more analyses. \cr \cr

CheckPloidy() summarizes all allele distance info into global matrices and .csv files for further analysis. \cr \cr

FindInformativeLoci() uses the info from CheckPloidy() to find which loci are informative based on some previously known individuals. \cr \cr

DeterminePloidy() is the final step which graphs agreeing informative loci by ploidy for each individual and stores that information in matrices. \cr \cr

MakeGraphs() graphs allelic distances for the two least smallest distances for each individual locus. points on those lines are individuals, green individuals are unknown, blue are tetraploids, and red are diploids. \cr \cr

DipIndex() calculates the relative levels of diploidization at each locus for each individual and stores this information in a .csv file. This can be useful for downstream analyses.
 }
\author{
\packageAuthor{PloidyPal}

Maintainer: \packageMaintainer{PloidyPal}\cr
	William Booker \cr Email: wbooker14@gmail.com
}

\keyword{ package }
\seealso{
See github site: https://github.com/wbooker/PloidyPal for examples of parameter files and to download example data

~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}


