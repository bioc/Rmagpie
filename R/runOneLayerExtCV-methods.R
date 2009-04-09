# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#-------------------------- runOneLayerExtCV-methods -------------------
# Run a one-layer external CV on a given assessment
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Run a one-layer external CV on the assessment <object>
#
# @param    object (assessment) assessment of interest
# @return   object (assessment) assessment of interest containing the results of
#               one-layer CV
#-------------------------------------------------------------------------------
setMethod("runOneLayerExtCV", "assessment", function(object) {
    generalRes <-   oneLayerExtCV(eset = getDataset(object),
                    noFolds = getNoFolds2ndLayer(object),
                    foldCreation = getTypeFoldCreation(object),
                    classifierName = getClassifierName(object),
                    featureSelectionMethod = getFeatureSelectionMethod(object),
                    #noTopGenes = getnoTopGene(object),
                    kernel = getSvmKernel(object),
                    noOfRepeats = getNoOfRepeats(object),
                    optionValues = getFeatureSelectionOptions(object, 'optionValues'),
                    verbose=getOption('verbose'))
    getResult1LayerCV(object) <- new("resultRepeated1LayerCV", generalRes)
    return(object)
    })