# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#--------------------------- rateGenes -----------------------------------------
# Methods and functions useful to generate the microarray-like image where
# a dot represent a gene and the brighter is the dot the more frequent it is.
#
# Author: Camille Maumet
# Creation: Jun. 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

# Create a palette with colors from bright yellow to black
myPalette <- rgb(red=1,green=1, blue=0.76)
for (i in 13:50){
  myPalette <- c(myPalette, rgb(red=1, blue=(50-i)/50, green=1))
}
for (i in 0:50){
  myPalette <- c(myPalette, rgb(red=(50-i)/50, blue=0, green=(50-i)/50))
}
myPalette <- rev(myPalette)

#-------------------------------------------------------------------------------
# Generate a matrix corresponding to the dot that will be displayed for gene
# with a mark <geneRate> on the microarray-like image
#
# @param    geneRate (numeric) mark obtained by the gene
#
# @return   matImage (matrix) to be plotted to display the gene
#-------------------------------------------------------------------------------
imageGene <- function(geneRate){
  # 81 pixels per imageGene
  nbPixels <- 9
  matImage <- matrix(c( rep(0, nbPixels),
                        rep(0, 3), rep(geneRate,3), rep(0, 3),
                        rep(0, 2), rep(geneRate,5), rep(0, 2),
                        rep(c(0, rep(geneRate, 7), 0),3),
                        rep(0, 2), rep(geneRate,5), rep(0, 2),
                        rep(0, 3), rep(geneRate,3), rep(0, 3),
                        rep(0, nbPixels)
                        ), dim<-c(nbPixels,nbPixels))
  return(matImage)
}

#-------------------------------------------------------------------------------
# Determine the frequency of each gene and generate a microarray-like image if
# <disp> is TRUE, in this case the image is stored at <storagePath>.
#
# @param    object (assessment) assessment of interest
#           storagePath (character) URL where the image must be stored, ignored
#               if disp is FALSE
#           disp (logical) TRUE if the microarray-image must be generated
#
# @return   rankedGenes (list) For each model, frequency of the genes ordered as
#           in geneNames (vector(character)) Names of the genes
#-------------------------------------------------------------------------------
setMethod(rankedGenesImg, signature="assessment",
    function(object, storagePath, disp=FALSE){
        # Names of the genes
        geneNames <- featureNames(getDataset(object))
        # Number of genes
        noGenes <- length(geneNames)
        # Values of the options
        optionValues <-  getFeatureSelectionOptions(object, topic='optionValues')
        rankedGenes <- alist()
        # Number of options
        noModels <- length(optionValues)

        # Generate an image for each option
        for (i in 1:noModels) {
            rankedGenes[[i]] <- vector("numeric", noGenes)
            names(rankedGenes[[i]]) <- geneNames
            frequencyTopGenes <- getResults(object, 1, 'genesSelected')[[i]]
            nbLists <- length(frequencyTopGenes)
            if (nbLists > 0){
                for (p in 1:nbLists){
                    genesList <- frequencyTopGenes[[p]]@genesList
                    noGenesInList <- length(genesList)
                    for (k in 1:noGenesInList){
                        currGene <- genesList[k]
                        rankedGenes[[i]][currGene] <- frequencyTopGenes[[p]]@frequ*100
                    }
                }

                if (disp){
                    printImageGenes(rankedGenes[[i]], storagePath, i)
                }
                rankedGenes[[i]] <- as.numeric(rankedGenes[[i]])
            }
        }
      return(list(ranks=rankedGenes, genesNames=geneNames))
})

#-------------------------------------------------------------------------------
# Generate the microarray-like image with genes of marks <rankedGenes> and save
# it at the path <storagePath$path> with the name:
# expe<storagePath$expeId>_chipImage<i>.png
#
# @param    rankedGenes (list) List of marks (frequency) of the genes
#           storagePath (list)
#               storagePath$path: URL where the image must be stored
#               storagePath$expeId: identifier of the current assessment
#           i (numeric) identifier of the current model
#-------------------------------------------------------------------------------
printImageGenes <- function(rankedGenes, storagePath, i){
	# By default display 30 genes per row
	defaultGenesPerRow <- 30
	geneImageSize <- 9

	sqrtNoGenes <- ceiling(sqrt(length(rankedGenes)))
	if ( sqrtNoGenes <= defaultGenesPerRow){
		nbGenePerRow <- defaultGenesPerRow
	} else {
		nbGenePerRow <- sqrtNoGenes
	}
	# Matrix must be squared => no of rows = no of columns
	nbRows <- nbGenePerRow
	matImage <- matrix(nrow=0, ncol=nbGenePerRow*geneImageSize)
	nbGenes <- length(rankedGenes)

	for (j in 1:nbRows) {
		# Matrix corresponding to the current row
		matImageRows <- matrix(nrow=geneImageSize, ncol=0)
		for (g in 1:nbGenePerRow){
		  geneNo <- g + ( (j-1)*nbGenePerRow)
		  if (geneNo <= nbGenes) {
		    oneGene <- imageGene(rankedGenes[geneNo])
		    matImageRows<- cbind(matImageRows, oneGene)
		  } else {
		    matImageRows<- cbind(matImageRows, matrix(rep(0, geneImageSize^2), dim<-c(geneImageSize,geneImageSize)))
		  }
		}
		matImage <- rbind(matImage, matImageRows)
	}
	matImage <- (t(matImage))[,dim(matImage)[2]:1]
	if (missing(storagePath)){
		par(oma=c(0,0,0,0))
		par(mar=c(0,0,0,0))
		image(matImage, col=myPalette, xaxt='n', yaxt='n', oma=c(0,0,0,0), mar=c(0,0,0,0))
	} else {
        png(paste(storagePath$path, "expe", storagePath$expeId, "_chipImage", i,".png",sep=""), width = 500, height = 500)
      		par(oma=c(0,0,0,0))
			par(mar=c(0,0,0,0))
			image(matImage, col=myPalette, xaxt='n', yaxt='n', oma=c(0,0,0,0), mar=c(0,0,0,0))
	    dev.off()
	}

}