CheckPloidy <- function(f){
infoTable <- as.matrix(read.csv(f, header=TRUE))
assign("infoTable",infoTable,globalenv())
SummaryFromFiles()
}


CalcAlleleDiffs <- function(f){
	infoTable <- as.matrix(read.csv(f, header=TRUE))
	BEG1 <- as.numeric(infoTable[1,2])
	END1 <- as.numeric(infoTable[2,2])
	str1 <- toString(infoTable[4,2])
	
	for(j in BEG1:END1){
	filePath <- paste(c(str1,"/I",j,"/I",j,"_allelesFromPost_4.txt"), collapse = "")	
	outFilePath <- paste(c(str1,"/I",j,"/I",j,"_allelesFromPost_4_diffs.txt"), collapse = "")	
	outTable <- matrix(nrow = 1, ncol = 8)	
	k <-scan(filePath, sep = ">", what = "complex")	
	for(i in seq(2,length(k),12)){
		textMat <- matrix(nrow = 6, ncol = 8)
		allele <- strsplit(k[i],".", fixed = TRUE)
		locus <- strsplit(allele[[1]][1],"L")[[1]][2]
		copy <- allele[[1]][2]
		alleleNumber <- allele[[1]][3]
		a1 <- k[i+1]
		a2 <- k[i+4]
		a3 <- k[i+7]
		a4 <- k[i+10]
		textMat[,c(1,4)] <- locus
		textMat[,c(2,5)] <- copy
		textMat[,3] <- c(1,1,1,2,2,3)
		textMat[,6] <- c(2,3,4,3,4,4)
		textMat[,7] <- c(adist(a1,a2),adist(a1,a3), adist(a1,a4), adist(a2,a3), adist(a2,a4), adist(a3,a4))
		textMat[,8] <- nchar(a1)
		outTable <- rbind(outTable,textMat)
	}
	outTable <- outTable[2:length(outTable[,1]),]
	write.table(outTable,file = outFilePath, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
	print(j)
	
	}
	assign("infoTable",infoTable,globalenv())
}


SummaryFromFiles <- function(){
BEG1 <- as.numeric(infoTable[1,2])
END1 <- as.numeric(infoTable[2,2])
str1 <- toString(infoTable[4,2])
assign("str1", str1, globalenv())
dir.create(paste(c(str1, "/PloidyRFiles"), collapse= ""), showWarnings = FALSE)
str1.1 <- "/I"
str2 <- "/I"
str3 <- "_allelesFromPost_4_diffs.txt"
infoFile <- toString(infoTable[6,2])
info <- read.csv(infoFile)
Notes <- table(info$PrePloidy)
LocusRef <- scan(infoTable[11,2])
testLoci <- sort(LocusRef)
numLoci <- length(testLoci)
refCheck <- as.numeric(infoTable[7,2])
alleleMat <- matrix(ncol=8)
if(refCheck == 1){
allMat <- as.matrix(read.csv(infoTable[8,2], header = FALSE) )
allDistancesMat <- as.matrix(read.csv(infoTable[9,2], header = FALSE))
prevInfo <- read.csv(infoTable[10,2], header = TRUE)
info <- rbind(info, prevInfo)
} else{
allMat <- matrix(ncol=5)
allDistancesMat <- matrix(data=c("Locus",testLoci),nrow=numLoci+1,ncol=1)
}
assign("info",info,globalenv())
########## This part summarizes all allele data into allMat ################
for (i in BEG1:END1){
	filePath <- paste(c(str1,str1.1,i,str2,i,str3), collapse = "")
	iNum <- paste(c("I",i), collapse = "")
	if (file.exists(filePath) == TRUE){
		dat <- read.table(filePath)	
		matLen <- length(dat[,7])/6
		valuesMat <- matrix(ncol=5)
		rowRef <- which(grepl(paste(c("^",iNum,"$"), collapse = ""), info$Individual))
		alleleCounter <- 1
		for (j in 1:matLen){
			if(dat[(alleleCounter),1] <= testLoci[length(testLoci)] & dat[(alleleCounter),1] %in% testLoci){
			alleles <- sort(c((dat[(alleleCounter),7]/dat[alleleCounter,8]), (dat[(alleleCounter+1),7]/dat[alleleCounter,8]), 					(dat[(alleleCounter+2),7]/dat[alleleCounter,8]), (dat[(alleleCounter+3),7]/dat[alleleCounter,8]), (dat[(alleleCounter+4),				7]/dat[alleleCounter,8]), (dat[(alleleCounter+5),7]/dat[alleleCounter,8])))
			allAlleleDists <- alleles		
			sameAlleles <- alleles[3:6]
			alleles <- alleles[1:2]
			
			valuesMat <- rbind(valuesMat, c(dat[alleleCounter,1],i,info$PrePloidy[rowRef], mean(alleles), sd(sameAlleles)))
			alleleMat <- rbind(alleleMat,c(i,dat[alleleCounter,1],allAlleleDists[1],allAlleleDists[2], 										allAlleleDists[3],allAlleleDists[4], allAlleleDists[5], allAlleleDists[6]))
			}
			else{}
			alleleCounter <- alleleCounter+6
			
		}
		valuesMat <- valuesMat[2:length(valuesMat[,1]),]
		allMat <- rbind(allMat, valuesMat)
		## tempVal <- mean(filtMat)
		indDist <- matrix(nrow=numLoci+1, ncol = 1)
		indDist[1,] = iNum
		for (k in 1:numLoci){
			M <- valuesMat[valuesMat[,1] == testLoci[k], 4]
				indDist[k+1,] <-mean(M)
			}	
		allDistancesMat <- cbind(allDistancesMat,indDist)
		print(paste(c(iNum, " COMPLETED"), collapse = ""), quote = FALSE)
	}
	
	else {
		print(paste(c(iNum," FILE NOT FOUND"), collapse = ""), quote = FALSE)
		}
		}
	
  




	write.table(allDistancesMat, file = paste(c(str1,"/PloidyRFiles/allDistancesMat.csv"), collapse=""), row.names=FALSE, col.names=FALSE, sep=",")
	cat(c("Writing:" ,paste(c(str1,"/PloidyRFiles/allDistancesMat.csv"), collapse = "") ,"\n"))

	allMat <- allMat[2:length(allMat[,1]),]  ####Important!

	write.table(allMat, file = paste(c(str1,"/PloidyRFiles/allMat.csv"), collapse=""), row.names=FALSE, col.names=FALSE, sep=",")
     cat(c("Writing:" , paste(c(str1,"/PloidyRFiles/allMat.csv"), collapse = ""), "\n"))
	alleleMat <- alleleMat[2:length(alleleMat[,1]),] 
	write.table(alleleMat, file = paste(c(str1,"/PloidyRFiles/alleleMat.csv"), collapse=""), row.names=FALSE, col.names=FALSE, sep=",")
	cat(c("Writing:" , paste(c(str1,"/PloidyRFiles/alleleMat.csv"), collapse = ""), "\n"))
	assign("allMat",allMat,globalenv())
	assign("alleleMat", alleleMat, globalenv())
	assign("testLoci",testLoci,globalenv())
	assign("numLoci", numLoci, globalenv())
	avgDist <- matrix(nrow=numLoci, ncol = 1)
	devDist <- matrix(nrow=numLoci, ncol = 1)
	  for (k in 1:numLoci){
		M <- allMat[allMat[,1] == testLoci[k], 3]
		avgDist[k,] <- mean(M)
		devDist[k,] <- sd(M)
	   }
}

ReadSumFromFile <- function (z){
	######### Sample Input:
	######### ReadSumFromFile("/Volumes/P0027_Gerhardt/P0027")
allMat <- as.matrix(read.csv(paste(c(z, "_allMat.csv"), collapse = ""), header = FALSE))
info <- read.csv((paste(c(z, "_RefTable.csv"), collapse = "")))
assign("allMat",allMat,globalenv())	
assign("info",info,globalenv())	
}

DipIndex <- function(){
	preInums <- unique(alleleMat[,1])
	tetNums <- array()
	dipNums <- array()	
	for(i in 1:length(preInums)){
		rowRef <- which(grepl(paste(c("^I",preInums[i],"$"), collapse = ""), info$Individual))
		if (info$PrePloidy[rowRef]==4){
		tetNums <- c(tetNums, preInums[i])
		}
		else if (info$PrePloidy[rowRef]==2){
		dipNums <- c(dipNums, preInums[i])
		}
		else{}
	}
	tetNums <- tetNums[2:length(tetNums)]
	dipNums <- dipNums[2:length(dipNums)]
	alleleAgreeMat <- matrix(nrow=length(tetNums)+6, ncol=numLoci+1)
	alleleAgreeMat[1,] <- c("Individual",testLoci)
	alleleAgreeMat[,1] <- c("Individual", tetNums, "Individuals with 4 alleles", "Percent Individuals Retained", "Retention against reference", "Avg", "Y/N Non-Dip")
	nonDipLoci <- array()
	
	refAgreeMat <- matrix(nrow=length(dipNums)+4, ncol=numLoci+1)
	refAgreeMat[1,] <- c("Individual",testLoci)
	refAgreeMat[,1] <- c("Individual", dipNums, "Individuals with 4 alleles", "Percent Individuals Retained", "Avg")
	library(Matrix)
	for (i in 1:length(tetNums)){
		tempAllele <- alleleMat[alleleMat[,1]==tetNums[i],]
		for (j in 1:numLoci){
			allAgree <- 0
 			if (testLoci[j] %in% tempAllele[,2] == TRUE){
			tempLoc <- matrix(ncol = length(tempAllele[1,]))
			tempLoc <- rbind(tempLoc, tempAllele[tempAllele[,2]==testLoci[j],])
			for(k in 1:(length(tempLoc[,1])-1)){
				#if(nnzero(tempLoc[k+1,3:8]) == 6){
				#	allAgree <- 1
				#}
				#else{}
				allAgree <- nnzero(tempLoc[k+1,3:8])
				if(allAgree == 0){
					
				}
				else{allAgree <- allAgree - 2}
			}
			}
			else{}
			alleleAgreeMat[i+1,j+1] <- allAgree
		}
		
		
	}
	
	
		for (i in 1:length(dipNums)){
		tempAllele <- alleleMat[alleleMat[,1]==dipNums[i],]
		for (j in 1:numLoci){
			allAgree <- 0
 			if (testLoci[j] %in% tempAllele[,2] == TRUE){
			tempLoc <- matrix(ncol = length(tempAllele[1,]))
			tempLoc <- rbind(tempLoc, tempAllele[tempAllele[,2]==testLoci[j],])
			for(k in 1:(length(tempLoc[,1])-1)){
				if(nnzero(tempLoc[k+1,3:8]) >= 4){
					allAgree <- 2
				}
				else{allAgree <- nnzero(tempLoc[k+1,3:8])
					if(allAgree == 0){}
			        else{allAgree <- allAgree - 2}
				}
			}
			}
			else{}
			refAgreeMat[i+1,j+1] <- allAgree
		}
		
		
	}
	
	for (i in 1:numLoci){
	
		tempCol <- refAgreeMat[2:(length(refAgreeMat[,1])-3), i+1] 
		refAgreeMat[length(refAgreeMat[,1])-2,i+1] <- sum(as.numeric(tempCol==2))
		refDiv <- sum(as.numeric(tempCol==2))/length(tempCol)
		refAgreeMat[length(refAgreeMat[,1])-1,i+1] <- refDiv
		refAgreeMat[length(refAgreeMat[,1]),i+1] <- mean(as.numeric(tempCol))
		projAlleleNum <- mean(as.numeric(tempCol))*2 ##### Really div by 2 and times 4 but you know
		
		tempCol <- alleleAgreeMat[2:(length(alleleAgreeMat[,1])-5), i+1] 
		alleleAgreeMat[length(alleleAgreeMat[,1])-4,i+1] <- sum(as.numeric(tempCol==4))
		tetDiv <- sum(as.numeric(tempCol==4))/length(tempCol)
		alleleAgreeMat[length(alleleAgreeMat[,1])-3,i+1] <- tetDiv
		retPercent <- if(tetDiv==0 & refDiv==0) {"Conserved Across"} else{tetDiv/refDiv}
		alleleAgreeMat[length(alleleAgreeMat[,1])-2,i+1] <- retPercent
		meanAlleleComps <- mean(as.numeric(tempCol))
		alleleAgreeMat[length(alleleAgreeMat[,1])-1,i+1] <- meanAlleleComps
		if(meanAlleleComps <= projAlleleNum+-0.1){
			alleleAgreeMat[length(alleleAgreeMat[,1]),i+1]	<- 0
		}
		else{
			alleleAgreeMat[length(alleleAgreeMat[,1]),i+1]	<- 1
			nonDipLoci <- c(nonDipLoci, alleleAgreeMat[1,i+1])
			
		}
	}
	nonDipLoci <- as.numeric(nonDipLoci[2:length(nonDipLoci)])
	assign("nonDipLoci", nonDipLoci, globalenv())
	write.table(nonDipLoci, file = paste(c(str1,"/PloidyRFiles/nonDipLoci.txt"), collapse = ""), row.names=FALSE, 		col.names=FALSE, sep=",")
	assign("alleleAgreeMat", alleleAgreeMat, globalenv())
	write.table(alleleAgreeMat, file = paste(c(str1,"/PloidyRFiles/DiploidizationIndex.csv"), collapse = ""), row.names=FALSE, 		col.names=FALSE, sep=",")	

	assign("refAgreeMat", refAgreeMat, globalenv())
	write.table(refAgreeMat, file = paste(c(str1,"/PloidyRFiles/RefDipIndex.csv"), collapse = ""), row.names=FALSE, 			   col.names=FALSE, sep=",")		
}

########### function for plotting individual allele distances, colored by ploidy ########
LocusPlot <- function(x){
	par(mfrow=c(1,1))
	locus <- x
	if (locus %in% allMat[,1]==TRUE){
		locusMat <- allMat[allMat[,1]==locus,]
		plot(1:length(locusMat[,1]),rep(0,length(locusMat[,1])), ylim=c(-.5,.5),cex=.01, main=toString(locus), ylab="Allele Distance", xlab = "Individual")
		for (i in 1:length(locusMat[,1])){
			points(i,locusMat[i,4], col = if(locusMat[i,3]==0){"green"} else if(locusMat[i,3]==2){"blue"} 						else{"red"},cex=1)
		}
		plot(1:length(locusMat[,1]),rep(0,length(locusMat[,1])), ylim=c(-.1,.1),cex=.01,main=toString(locus), ylab="Allele Distance", xlab = "Individual")
		for (i in 1:length(locusMat[,1])){
			points(i,locusMat[i,4], col = if(locusMat[i,3]==0){"green"} else if(locusMat[i,3]==2){"blue"} 						else{"red"},cex=1)
		}
	}

}

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

FindInformativeLoci <- function(outlierThresh, numOutliers, informativeThresh) {
######## Function to plot and calculate informative loci based on outlier threshold and difference threshold #####

######## Defaults ##########
	#  outlierThresh <- 0.1  ### level of divergence determining that individual as outlier
	#  numOutliers <- 4      ### number of outliers allowed per locus
	#  informativeThresh <- 0.0055  

	if(missing(outlierThresh)){outThresh <- 0.1} else {outThresh <- outlierThresh}
	if(missing(numOutliers)){numOutThresh <- 4} else {numOutThresh <- numOutliers}
	if(missing(informativeThresh)){infoThresh <- 0.0055} else {infoThresh <- informativeThresh}
	alleleCheckMat <- matrix(nrow=numLoci,ncol=2)
	sdOutsCheckMat <- matrix(nrow=numLoci, ncol=2)
	sdDifs <- matrix(nrow=numLoci,ncol=1)  
	outMat <- matrix(nrow=numLoci,ncol=2)
	dipTemp <- c()
	tetTemp <- c()
	for (i in 1:numLoci){
		locusMat <- allMat[allMat[,1]==testLoci[i],]
		sdDifs[i] <- sd(locusMat[,4])
		dipTemp <- locusMat[locusMat[,3]==2,4]
		dipSdTemp <- locusMat[locusMat[,3]==2,5]
		outMat[i,1] <- sum(dipTemp > outThresh)
		dipTemp <- dipTemp[dipTemp <= outThresh]	 
		tetTemp <- locusMat[locusMat[,3]==4,4]
		tetSdTemp <- locusMat[locusMat[,3]==4,5]
		outMat[i,2] <- sum(tetTemp > outThresh)
		tetTemp <- tetTemp[tetTemp <= outThresh]
		alleleCheckMat[i,1] <- mean(dipTemp)
		alleleCheckMat[i,2] <- mean(tetTemp)
		sdOutsCheckMat[i,1] <- mean(dipSdTemp)
		sdOutsCheckMat[i,2] <- mean(tetSdTemp)
	}
	alleleDifMat <- alleleCheckMat[,2] - alleleCheckMat[,1]
	iLoci <- c()	
	InfsPDFPath = paste (c(str1,"/PloidyRFiles/InformativeLocigraph.pdf"), collapse = "")
     pdf(file = InfsPDFPath )
	plot(1:testLoci[length(testLoci)],rep(0,testLoci[length(testLoci)]), ylim=c(-.025,.03),cex=.01, xlab = "Locus", ylab = "Level of Divergence")
	for (i in 1:numLoci){
		if(sum(outMat[i,]) <= numOutThresh){
			if(alleleDifMat[i] != "NaN"){
				if(alleleDifMat[i] >= infoThresh){
					color <- "blue"
					iLoci<- rbind(iLoci, as.numeric(testLoci[i]))
					text(testLoci[i]+9,alleleDifMat[i],testLoci[i])
				}
				else{
					color <- "black"
				}
				}
			else{}
		
		} 
		else{
			color = "red"
		}
	points(testLoci[i], alleleDifMat[i], col = color)
	}

##################################   RANKING   ##############################################################		
rankMat	<- matrix(nrow=length(iLoci),ncol=2)
for (i in 1:length(iLoci)){
	rankMat[i,1] <- iLoci[i,]
	rankMat[i,2] <- alleleDifMat[which(testLoci==iLoci[i,])]	
}
rankMat <- rankMat[sort.list(rankMat[,2],decreasing=TRUE),]
cutMat <- rankMat	
assign("iLoci", iLoci, globalenv())
assign("alleleCheckMat",alleleCheckMat,globalenv())	
assign("alleleDifMat",alleleDifMat,globalenv())	
assign("cutMat", cutMat, globalenv())
assign("rankMat", cutMat, globalenv())
assign("sdOutsCheckMat", sdOutsCheckMat, globalenv())
write.table(sdOutsCheckMat, file = paste(c(str1,"/PloidyRFiles/sdOutsCheckMat.csv"), collapse = ""), row.names=FALSE, 		col.names=FALSE, sep=",")
dev.off()
cat(c("Writing:" , paste(c(str1,"/PloidyRFiles/InformativeLociGraph.pdf"), collapse = ""), "\n")) 
}



CutInfs <- function(len){
cutMat <- rankMat[1:len,]
assign("cutMat", cutMat, globalenv())
}

UncutInfs <- function(){
cutMat <- rankMat
assign("cutMat", cutMat, globalenv())
}



DeterminePloidy <- function(meanThresh, sdThresh, inds){
########### Ind x Locus agreement with informative ######################################
if(missing(inds)){} else {cutMat <- rankMat[1:inds,]}
if(missing(meanThresh)){indThresh <- .75} else {indThresh <- meanThresh}
if(missing(sdThresh)){sThresh <- .8} else {sThresh <- sdThresh}

meanLocusDif <- (alleleCheckMat[,2]+alleleCheckMat[,1])/2
individuals <- unique(allMat[,2])
infTestMat <- matrix(nrow=length(individuals)+1, ncol=length(cutMat[,1])+1)
infTestMat[,1] <- c(0,individuals)
infTestMat[1,] <- c(0,cutMat[,1])
ratioMat <- matrix (nrow=length(individuals), ncol = 3)
totalsMat <- matrix(nrow=length(individuals),ncol=3)
for (i in 1:length(individuals)){
	tempMeanSum <- 0

	rowRef <- which(grepl(paste(c("^I",individuals[i],"$"), collapse = ""), info$Individual))

	tSum <- 0
	fSum <- 0
	iCount = i+1
	indMat <- allMat[allMat[,2]==individuals[i],]
	for (j in 1:length(cutMat[,1])){
		jCount = j+1
         locusMean <- meanLocusDif[which(testLoci==cutMat[j,1])]
          tetSd <- sdOutsCheckMat[which(testLoci==cutMat[j,1]),2]
		indLocMean <- mean(indMat[indMat[,1]==cutMat[j,1],4])
		indSdofOuts <- mean(indMat[indMat[,1]==cutMat[j,1],5])
		tempMeanSum <- tempMeanSum + indLocMean
		if (indLocMean != "NaN"){
			if (indLocMean >= locusMean+-(locusMean*indThresh) & indSdofOuts >=tetSd+-(tetSd*sThresh)){
				infTestMat[iCount,jCount] <- 1
				tSum <- tSum +1
			}
			else{
				infTestMat[iCount,jCount] <- 0
				fSum <- fSum +1
			}
		}	
		
		else{}
	}
	
			ratioMat[i,] <- c(individuals[i],tSum,info$PrePloidy[rowRef])
		
	
}

plot(ratioMat[,3]+runif(length(ratioMat[,3]),-0.1,0.1),ratioMat[,2], xlab="Ploidy", ylab="Number of Agreeing Informative Loci")

tetMin <- 1
dipMax <-1
tetMean <-1
dipMean <-1

tets <- ratioMat[ratioMat[,3]==4,]
avgAgree <- mean(tets[,2])
sdAgree <- sd(tets[,2])
badTet <- tets[tets[,2]<(avgAgree-(1.96*sdAgree)),]

dips <- ratioMat[ratioMat[,3]==2,]
badDip <- dips[dips[,2]>(avgAgree-(1.96*sdAgree)),]

unknowns <- ratioMat[ratioMat[,3]==0,]

assign("badTet", badTet, globalenv())
assign("badDip", badDip, globalenv())
assign("unknowns", unknowns, globalenv())
assign("infMat", ratioMat, globalenv())
}





