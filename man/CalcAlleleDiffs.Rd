\name{CalcAlleleDiffs}
\alias{CalcAlleleDiffs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate allele differences
}
\description{
This goes through individual text files for each specimen and calculates the distances between all alleles at a given locus or copy of that locus.
}
\usage{
CalcAlleleDiffs("/CalcParamFile.csv")
}
\details{
Parameter file should be a .csv with a header with four lines and two columns (see https://github.com/wbooker/PloidyPal for examples), names are in the left column and vlues in the right. rows are in order as follows: Begin number, end number, ploidy level, and root filepath. 
}

\note{
	Files should be in this format: \cr \cr
	---Individuals should be in seperate folders inside the root filepath with numbers followed by the letter I, as in 'I3756' \cr \cr
	---Files should follow a similar format and be titled as thus: "I3756_allelesFromPost_4.txt"\cr \cr
	---Each file should have alleles formatted thusly (see github for clear example): \cr \cr
	
>L1.1.1\cr
AATTCCGGAA\cr
>L1.1.2\cr
AATTCGCCAA\cr

Where the first number is the locus, second is the copy, and third is the allele number

}

\author{
William Booker
Florida State University Dept. of Biological science
}


