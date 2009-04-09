# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#-------------------------- runTwoLayerExtCV-methods -------------------
# Run a two-layer external CV on a given assessment
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Run a two-layer external CV on the assessment <object>
#
# @param    object (assessment) assessment of interest
# @return   object (assessment) assessment of interest containing the results of
#               two-layer CV
#-------------------------------------------------------------------------------
setMethod("runTwoLayerExtCV", "assessment", function(object) {
    getResult2LayerCV(object) <- twoLayerExtCV(eset = getDataset(object),
                  noFoldsExtLayer = getNoFolds2ndLayer(object),
                  noFoldsIntLayer = getNoFolds1stLayer(object),
                  foldCreation = getTypeFoldCreation(object),
                  classifierName = getClassifierName(object),
                  featureSelectionMethod = getFeatureSelectionMethod(object),
                  #noTopGenes = getNoTopGene(object),
                  kernel = getSvmKernel(object),
                  noOfRepeats = getNoOfRepeats(object),
                  optionValues = getFeatureSelectionOptions(object, 'optionValues'),
                  verbose=getOption('verbose'))
    return(object)
    })