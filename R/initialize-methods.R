# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#--------------------------- initialize-methods --------------------------------
# Initialization methods of the objects
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

# geneSubsets
setMethod("initialize", "geneSubsets",
    function(.Object, optionValues, maxSubsetSize, speed=.Object@speed) {
        if (missing(optionValues)){
            if (missing("maxSubsetSize")){
                stop("In 'intialize', With 'rfe', 'optionValues' or 'maxSubsetSize' must be given")
            } else {
                    if (length(maxSubsetSize) != 1 || maxSubsetSize <= 1){
                        stop("'maxSubsetSize' must be > 1")
                    }
                    .Object@maxSubsetSize <- maxSubsetSize
                    .Object@speed <- speed
                    if (speed == "high"){
                        noOfOptions <- ceiling(log2(maxSubsetSize)+1)
                        # Number of model (subset of genes) to be tested
                        .Object@noOfOptions <- noOfOptions
                        # noSelectedFeature contains the different number of genes that we will try
                        .Object@optionValues <- 2^(0:(noOfOptions-1))
                        .Object@optionValues[noOfOptions] <- maxSubsetSize

                    } else{
                        .Object@noOfOptions <- maxSubsetSize
                        .Object@optionValues <- seq(from = 1, to = maxSubsetSize)
                    }
            }
        } else {
            if ( ! missing("maxSubsetSize")){
                stop("'optionValues' and 'maxSubsetSize' can't be both defined")
            }
            .Object@optionValues <- sort(optionValues)
            .Object@noOfOptions <- length(optionValues)
            if (length(optionValues) != 0){
                .Object@maxSubsetSize <- max(optionValues)
            }
            else {
                # Wrong value since optionValues is invalid
                .Object@maxSubsetSize <- 0
            }
            if (.Object@noOfOptions < .Object@maxSubsetSize ){
                .Object@speed <- "high"
            } else {
                .Object@speed <- "slow"
            }
        }
        validObject(.Object)
        .Object
    }
)

# thresholds
setMethod("initialize", "thresholds",
    function(.Object, optionValues) {
        if (missing(optionValues)){
            .Object@optionValues <- numeric(0)
        } else {
            if (length(optionValues) == 0){
                stop("thresholds: 'optionValue' must be a vector of length >= 1")
            }
            if (any(optionValues< 0)){
                stop("thresholds: 'optionValue' must be a vector of positive values")
            }
            compaValue <- 0
            for (i in 1:length(optionValues)){
                if (! optionValues[i] >= compaValue){
                    stop("thresholds: 'optionValue' must be in ascending orde")
                }
                compaValue <- optionValues[i]
            }
            .Object@optionValues <- optionValues
            .Object@noOfOptions <- length(optionValues)
        }
        .Object@noOfOptions <- length(.Object@optionValues)
        # Do not valid the object here, wait its integration in assessment
        # where the values of the thresholds will be determined
        .Object
    }
)

# assessment
setMethod("initialize", "assessment", 
    function(.Object, dataset, 
              featureSelectionOptions=.Object@featureSelectionOptions,
              noFolds1stLayer=.Object@noFolds1stLayer,
              noFolds2ndLayer=.Object@noFolds2ndLayer,
              classifierName=.Object@classifierName,
              featureSelectionMethod=.Object@featureSelectionMethod,
              typeFoldCreation=.Object@typeFoldCreation,
              svmKernel=.Object@svmKernel, 
              noOfRepeats=.Object@noOfRepeats,
              resultRepeated1LayerCV=.Object@resultRepeated1LayerCV,
              resultRepeated2LayerCV=.Object@resultRepeated2LayerCV) {
          .Object@dataset <- dataset
          .Object@noFolds1stLayer <- noFolds1stLayer
          .Object@noFolds2ndLayer <- noFolds2ndLayer
          .Object@classifierName <- classifierName
          .Object@featureSelectionMethod <- featureSelectionMethod
          .Object@typeFoldCreation <- typeFoldCreation
          if ( svmKernel=='' && classifierName!='svm'){
            # If classifierNames is not SVM the RFE is done with a linear SVM
            .Object@svmKernel <- 'linear'
          } else {
            .Object@svmKernel <- svmKernel
          }
          .Object@noOfRepeats <- noOfRepeats
          .Object@resultRepeated1LayerCV <- resultRepeated1LayerCV
          .Object@resultRepeated2LayerCV <- resultRepeated2LayerCV
          if (! missing("featureSelectionOptions")){
            .Object@featureSelectionOptions <- featureSelectionOptions
            # Set the default thresholds
            if (class(.Object@featureSelectionOptions) == 'thresholds' && length(.Object@featureSelectionOptions@optionValues) == 0){
                trained <-  pamr.train(data=list(  x=exprs(dataset),
                                            y=pData(dataset)[[1]]))
                .Object@featureSelectionOptions@optionValues <- trained$threshold
                .Object@featureSelectionOptions@noOfOptions <- length(trained$threshold)
                # Valid the 'thresholds' object
                validObject(.Object@featureSelectionOptions)
            }
          } else {
            if (length(featureSelectionMethod) > 0){
                if (.Object@featureSelectionMethod == 'rfe'){
                    if (! is.null(.Object@dataset) ){
                       .Object@featureSelectionOptions <- new("geneSubsets", maxSubsetSize=length(featureNames(.Object@dataset)))
                    } else { # Wrong dataset, just put random data in genSubsets
                      .Object@featureSelectionOptions <- new("geneSubsets", maxSubsetSize=10)
                    }
                } else {
                    # Set the default thresholds
                    trained <-  pamr.train(data=list(   x=exprs(dataset),
                                                        y=pData(dataset)[[1]]))
                    .Object@featureSelectionOptions <- new("thresholds", optionValues=trained$threshold)
                    validObject(.Object@featureSelectionOptions)
                }
            }
        }
          validObject(.Object)
          .Object })