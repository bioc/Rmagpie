# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#-------------------------- AllClasses -----------------------------------------
# This file contains all the class definitions
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

############################ Private (not in namespace) ########################
# Store the cross-validated error rate of a one-layer CV
setClass("cvErrorRate", representation(cvErrorRate = "numeric",
                                       seErrorRate = "numeric",
                                       classErrorRates = "matrix"))
# Store the error rates of a one-layer CV
setClass("errorRate1stLayerCV",
              representation( errorRatePerFold = "matrix",
                              noSamplesPerFold = "numeric",
                              cvErrorRate = "cvErrorRate"))

# Store the genes selected in a one-layer CV for each option (threshold or size of subset)
setClass("selectedGenesPerOneOption", representation("list",
                                                    optionValue="numeric"))

# Store a frequency and the corresponding list of genes
setClass("frequencyGenes",
          representation( frequ="numeric",
                          genesList="character"))

# Store a list frequencies with their corresponding list of genes
setClass("frequencyTopGenePerOneModel", "list")

# Store the genes selected in a one-layer CV (per model and fold and per frequency)
setClass("selectedGenes1stLayerCV",
              representation( selectedGenesPerFold = "list",
                              frequencyTopGene = "list"))

# Store the results of a single one-layer CV
setClass("resultSingle1LayerCV",
              representation( errorRates = "errorRate1stLayerCV",
                              selectedGenes = "selectedGenes1stLayerCV",
                              executionTime = "numeric",
                              bestOptionValue = "numeric"))

# Store the results of a repeated one-layer CV
setClass("resultRepeated1LayerCV",
              representation(   original1LayerCV = "list",
                                summaryErrorRate = "cvErrorRate",
                                summaryFrequencyTopGenes = "list",
                                bestOptionValue = "numeric",
                                executionTime = "numeric"))

# Store the cross-validated error rate of a two-layer CV
setClass("cvErrorRate2ndLayer",
            representation( seFinalErrorRate = "numeric",
                            finalErrorRate = "numeric",
                            classErrorRates = "numeric"))

# Store the error rates of a two-layer CV
setClass("errorRate2ndLayerCV",
            representation( errorRatePerFold = "numeric",
                            noSamplesPerFold = "numeric",
                            cvErrorRate = "cvErrorRate2ndLayer"))

# Store the genes selected in a fold of a two-layer CV and the corresponding option and number of genes
setClass("selectedGenes",
          representation( optionValue="numeric",
                          noOfGenes="numeric",
                          genesList="matrix"))

# Store the genes selected in a single two-layer CV
setClass("selectedGenes2ndLayerCV", "list")

# Store the genes selected in a single two-layer CV
setClass("result2LayerCV",
          representation(   results1stLayer = "list",
                            errorRates = "errorRate2ndLayerCV",
                            selectedGenes = "selectedGenes2ndLayerCV",
                            avgBestOptionValue = "numeric",
                            executionTime ="numeric"))

# Store the genes selected in a repeated two-layer CV
setClass("resultRepeated2LayerCV",
          representation(   original2LayerCV = "list",
                            summaryErrorRate = "cvErrorRate2ndLayer",
                            avgBestOptionValue = "numeric",
                            executionTime ="numeric"))


## Store the microarray data or NULL if it has not been loaded yet
#setClassUnion("ExpressionSetOrNull", c("ExpressionSet", "NULL"))
#
################################# Public ######################################
## Store the microarray data and its related files
#setClass("dataset", representation( dataId = "character",
#                                    dataPath = "character",
#                                    geneExprFile = "character",
#                                    classesFile = "character",
#                                    eset = "ExpressionSetOrNull"
#                                    ),
#                    prototype( dataPath="." ) )

# Store the results of a repeated one-layer CV or NULL if it has not been computed yet
setClassUnion("resultRepeated1LayerCVOrNULL",c("resultRepeated1LayerCV", "NULL"))

# Store the results of a repeated two-layer CV or NULL if it has not been computed yet
setClassUnion("resultRepeated2LayerCVorNULL",c("resultRepeated2LayerCV", "NULL"))

# Store the options related to an assessment (sizes of subset of thresholds)
setClass("featureSelectionOptions", representation( optionValues = "numeric",
                                        noOfOptions = "numeric", "VIRTUAL"))

# Store the size of subsets for RFE
setClass("geneSubsets", representation( "featureSelectionOptions",
                                        maxSubsetSize = "numeric",
                                        speed = "character"),
                       prototype( speed = "high"))

# Store the thresholds for NSC
setClass("thresholds", representation( "featureSelectionOptions"))

# Store the final classifier corresponding to an assessment (can then be used to
# classify new samples)
setClass("finalClassifier", representation( genesFromBestToWorst = "character",
                                            models = "list" ))

setIs("geneSubsets", "featureSelectionOptions")
setIs("thresholds", "featureSelectionOptions")

# Final classifier or NULL if it has not been computed yet
setClassUnion("finalClassifierOrNULL",c("finalClassifier", "NULL"))

# Store the options and results of an assessment
setClass("assessment", representation(  dataset = "ExpressionSet",
                                        noFolds1stLayer = "numeric",
                                        noFolds2ndLayer = "numeric",
                                        #noTopGene = "numeric",
                                        classifierName = "character",
                                        featureSelectionMethod = "character",
                                        typeFoldCreation = "character",
                                        svmKernel = "character",
                                        noOfRepeats = "numeric",
                                        featureSelectionOptions = "featureSelectionOptions",
                                        resultRepeated1LayerCV = "resultRepeated1LayerCVOrNULL",
                                        resultRepeated2LayerCV = "resultRepeated2LayerCVorNULL",
                                        finalClassifier = "finalClassifierOrNULL"
                                       ),
                       prototype(   noFolds1stLayer=10,
                                    noFolds2ndLayer=9,
                                    classifierName="svm",
                                    featureSelectionMethod="rfe",
                                    typeFoldCreation="original",
                                    svmKernel="linear",
                                    noOfRepeats=2,
                                    resultRepeated1LayerCV=NULL,
                                    resultRepeated2LayerCV=NULL,
                                    finalClassifier=NULL ))