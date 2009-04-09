# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#-------------------------- AllGenerics ----------------------------------------
# This file contains all the generic functions declaration definitions
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

############################ Private (not in namespace) ########################
# Display of an object (with print or show)
setGeneric("showWithPrefix", function(object, prefix="", short=FALSE) standardGeneric("showWithPrefix"))

# load data in a dataset
#setGeneric("loadData", function(object, pathGeneExpr, pathClasses, miame, annotation) standardGeneric("loadData"))

# Type checking
setGeneric("test.rateValue", function(attValue, attName, greaterThan0) standardGeneric("test.rateValue"))
setGeneric("test.integer", function(attValue, attName) standardGeneric("test.integer"))
setGeneric("test.positiveInteger", function(attValue, attName, greaterThan0) standardGeneric("test.positiveInteger"))
setGeneric("test.positiveFloat", function(attValue, attName, greaterThan0) standardGeneric("test.positiveFloat"))
setGeneric("test.orderedPositiveInteger", function(attValue, attName, greaterThan0) standardGeneric("test.orderedPositiveInteger"))
setGeneric("test.orderedPositiveFloat", function(attValue, attName, greaterThan0) standardGeneric("test.orderedPositiveFloat"))

# ---------------------- PRIVATE Accessors ------------------
# resultRepeated1LayerCV
setGeneric("getExecutionTime", function(object) standardGeneric("getExecutionTime"))
setGeneric("getSummaryErrorRate", function(object) standardGeneric("getSummaryErrorRate"))
setGeneric("getOriginal1LayerCV", function(object) standardGeneric("getOriginal1LayerCV"))
setGeneric("getBestOptionValue", function(object) standardGeneric("getBestOptionValue"))
setGeneric("getFrequencyTopGenes", function(object) standardGeneric("getFrequencyTopGenes"))

# resultSingle1LayerCV
setGeneric("getErrorRates", function(object) standardGeneric("getErrorRates"))
setGeneric("getExecutionTime", function(object) standardGeneric("getExecutionTime"))
setGeneric("getBestOptionValue", function(object) standardGeneric("getBestOptionValue"))

# selectedGenes1PerOneOption
setGeneric("getOptionValue", function(object) standardGeneric("getOptionValue"))
setGeneric("getGenesList", function(object) standardGeneric("getGenesList"))

# selectedGenes1stLayerCV
setGeneric("getSelectedGenesPerFold", function(object) standardGeneric("getSelectedGenesPerFold"))
setGeneric("getFrequencyTopGene", function(object) standardGeneric("getFrequencyTopGene"))

# errorRate1stLayerCV
setGeneric("getErrorRatePerFold", function(object) standardGeneric("getErrorRatePerFold"))
setGeneric("getCvErrorRate", function(object) standardGeneric("getCvErrorRate"))
setGeneric("getSeErrorRate", function(object) standardGeneric("getSeErrorRate"))
setGeneric("getClassErrorRates", function(object) standardGeneric("getClassErrorRates"))
setGeneric("getNoSamplesPerFold", function(object) standardGeneric("getNoSamplesPerFold"))


# frequencyGenes
setGeneric("getFrequency", function(object) standardGeneric("getFrequency"))
setGeneric("getGenesList", function(object) standardGeneric("getGenesList"))

# resultRepeated2LayerCV
setGeneric("getAvgBestOptionValue", function(object) standardGeneric("getAvgBestOptionValue"))
#setGeneric("getFinalErrorRate", function(object) standardGeneric("getFinalErrorRate"))
setGeneric("getOriginal2LayerCV", function(object) standardGeneric("getOriginal2LayerCV"))

# result2LayerCV
setGeneric("getResults1stLayer", function(object) standardGeneric("getResults1stLayer"))
setGeneric("getErrorRates", function(object) standardGeneric("getErrorRates"))
setGeneric("getSelectedGenes", function(object) standardGeneric("getSelectedGenes"))
setGeneric("getExecutionTime", function(object) standardGeneric("getExecutionTime"))
setGeneric("getAvgBestOptionValue", function(object) standardGeneric("getAvgBestOptionValue"))

# errorRate2ndLayerCV
setGeneric("getErrorRatePerFold", function(object) standardGeneric("getErrorRatePerFold"))
setGeneric("getFinalErrorRate", function(object) standardGeneric("getFinalErrorRate"))

# selectedGenes
setGeneric("getNoOfGenes", function(object) standardGeneric("getNoOfGenes"))
setGeneric("getGenesList", function(object) standardGeneric("getGenesList"))

################################ Public ########################################
# Run a one-layer CV
setGeneric("runOneLayerExtCV", function(object) standardGeneric("runOneLayerExtCV"))
# Run a two-layer CV
setGeneric("runTwoLayerExtCV", function(object) standardGeneric("runTwoLayerExtCV"))
# Plot results of repeated one-layer CV with all the repeats
setGeneric("plotErrorsRepeatedOneLayerCV", function(object) standardGeneric("plotErrorsRepeatedOneLayerCV"))
# Plot results of repeated one-layer CV (summary only)
setGeneric("plotErrorsSummaryOneLayerCV", function(object) standardGeneric("plotErrorsSummaryOneLayerCV"))
# Plot results of repeated two-layer CV (summary only)
setGeneric("plotErrorsFoldTwoLayerCV", function(object) standardGeneric("plotErrorsFoldTwoLayerCV"))
# Plot the microarray-like image
setGeneric("rankedGenesImg", function(object, storagePath, disp=FALSE) standardGeneric("rankedGenesImg"))
# Classify new samples
setGeneric("classifyNewSamples", function(object, newSamplesFile, optionValue=0) standardGeneric("classifyNewSamples"))
# Train the final classifier
setGeneric("findFinalClassifier", function(object) standardGeneric("findFinalClassifier"))

# ---------------------- Accessors ---------------------------------------------
# ** geneSubsets accessors
setGeneric("getMaxSubsetSize", function(object) standardGeneric("getMaxSubsetSize"))
setGeneric("getOptionValues", function(object) standardGeneric("getOptionValues"))
setGeneric("getSubsetsSizes", function(object) standardGeneric("getSubsetsSizes"))
setGeneric("getNoThresholds", function(object) standardGeneric("getNoThresholds"))
setGeneric("getNoOfOptions", function(object) standardGeneric("getNoOfOptions"))
setGeneric("getOptionValues", function(object) standardGeneric("getOptionValues"))
setGeneric("getNoModels", function(object) standardGeneric("getNoModels"))
setGeneric("getSpeed", function(object) standardGeneric("getSpeed"))
setGeneric("getMaxSubsetSize<-", function(object, value) standardGeneric("getMaxSubsetSize<-"))
setGeneric("getSubsetsSizes<-", function(object, value) standardGeneric("getSubsetsSizes<-"))
setGeneric("getNoThresholds<-", function(object, value) standardGeneric("getNoThresholds<-"))
setGeneric("getOptionValues<-", function(object, value) standardGeneric("getOptionValues<-"))
setGeneric("getSpeed<-", function(object, value) standardGeneric("getSpeed<-"))

# ** finalClassifier accessors
setGeneric("getGenesFromBestToWorst", function(object) standardGeneric("getGenesFromBestToWorst"))
setGeneric("getModels", function(object) standardGeneric("getModels"))
setGeneric("getGenesFromBestToWorst<-", function(object, value) standardGeneric("getGenesFromBestToWorst<-"))
setGeneric("getModels<-", function(object, value) standardGeneric("getModels<-"))

# ** assessment accessors
setGeneric("getDataset", function(object) standardGeneric("getDataset"))
setGeneric("getFeatureSelectionOptions", function(object, topic) standardGeneric("getFeatureSelectionOptions"))
setGeneric("getNoFolds1stLayer", function(object) standardGeneric("getNoFolds1stLayer"))
setGeneric("getNoFolds2ndLayer", function(object) standardGeneric("getNoFolds2ndLayer"))
setGeneric("getClassifierName", function(object) standardGeneric("getClassifierName"))
setGeneric("getFeatureSelectionMethod", function(object) standardGeneric("getFeatureSelectionMethod"))
setGeneric("getSvmKernel", function(object) standardGeneric("getSvmKernel"))
setGeneric("getTypeFoldCreation", function(object) standardGeneric("getTypeFoldCreation"))
setGeneric("getNoOfRepeats", function(object) standardGeneric("getNoOfRepeats"))
setGeneric("getFinalClassifier", function(object, topic) standardGeneric("getFinalClassifier"))
setGeneric("getResult1LayerCV", function(object) standardGeneric("getResult1LayerCV"))
setGeneric("getResult2LayerCV", function(object) standardGeneric("getResult2LayerCV"))
setGeneric("getResults", function(object, layer, ...) standardGeneric("getResults"))
setGeneric("getResultsLayer", function(object, topic, genesType, errorType) standardGeneric("getResultsLayer"))

setGeneric("getDataset<-", function(object, value) standardGeneric("getDataset<-"))
setGeneric("getFeatureSelectionOptions<-", function(object, topic, value) standardGeneric("getFeatureSelectionOptions<-"))
setGeneric("getNoFolds1stLayer<-", function(object, value) standardGeneric("getNoFolds1stLayer<-"))
setGeneric("getNoFolds2ndLayer<-", function(object, value) standardGeneric("getNoFolds2ndLayer<-"))
setGeneric("getClassifierName<-", function(object, value) standardGeneric("getClassifierName<-"))
setGeneric("getSvmKernel<-", function(object, value) standardGeneric("getSvmKernel<-"))
setGeneric("getTypeFoldCreation<-", function(object, value) standardGeneric("getTypeFoldCreation<-"))
setGeneric("getNoOfRepeats<-", function(object, value) standardGeneric("getNoOfRepeats<-"))
setGeneric("getFinalClassifier<-", function(object, value) standardGeneric("getFinalClassifier<-"))
setGeneric("getResult1LayerCV<-", function(object, value) standardGeneric("getResult1LayerCV<-"))
setGeneric("getResult2LayerCV<-", function(object, value) standardGeneric("getResult2LayerCV<-"))

