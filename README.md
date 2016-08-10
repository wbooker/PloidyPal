# PloidyPal
R Package and Files for polyploid analysis

Install in R by running install_github("wbooker/PloidyPal") with devtools package installed 

PloidyPal completes a number of functions that can assist in and polyploid genetic analysis. Mainly, it helps determine the ploidy level of individuals in a dataset using allele phase data. It can also give a rough estimate of the level of diploidization that has occured at each locus for each individual. The data required are individual .txt files of phased alleles for individuals--with 4 phased alleles for individuals regardless of ploidy level. This is cruicial to the functioning of the program.

All data should be stored in a root directory which is referred to throughout use of the program.

You can download the example folder in the github directory which should allow running of an example dataset instantly (once the paths are updated in the csv files).

Files should be stored in this heirarchy:


    >/rootpath
      >/rootpath/I001         ###############Folder for each individual, predicated by I and followed by a number
          >/rootpath/I001/I001_allelesFromPost_4.txt   ############## Allele phase data for four alleles for that individual

Inside the root path there should be one or two parameter files in .csv format (if necessary for different functions) as well as two other files. A .csv or .txt of the loci in the data set you want to use in analysis, and a .csv or .txt file of the individual numbers (without the I) and their ploidy level (2 for diploid, 4 for tetraploid, and 0 for unknown).

As of 08/05/2016, only diploid and tetraploid datasets are functional. 

Currently, the package requires some a priori knowledge of at least a few individuals in the data set. It is currently tailored for  Anchored Phylogenomics data simply because that is the format I had available. Tailoring to individual datasets should be relatively easy though, or, you can reformat to the examples shown.



CalcAlleleDiffs() takes files with 4 alleles per locus and calculates the differences between alleles for all alleles at all loci for each individual. These data are then stored in text files that is used for many more analyses. 

SummaryFromFiles() summarizes all allele distance info into global matrices and .csv files for further analysis. 

FindInformativeLoci() uses the info from CheckPloidy() to find which loci are informative based on some previously known individuals. 

DeterminePloidy() is the final step which graphs agreeing informative loci by ploidy for each individual and stores that information in matrices. 

MakeGraphs() graphs allelic distances for the two least smallest distances for each individual locus. points on those lines are individuals, green individuals are unknown, blue are tetraploids, and red are diploids. 

DipIndex() calculates the relative levels of diploidization at each locus for each individual and stores this information in a .csv file. This can be useful for downstream analyses.


If you would like to use this package on your data, you can either get everything into the correct compatable filetypes and formats which should allow you to run it instantly, or, you can email wbooker14@gmail.com to request tailoring the package to work with your files. Availability for this is unknown and dependent on how much your data deviate from the example, but I will try to work on these when I can. 

-William Booker
