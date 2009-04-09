# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#--------------------------- validation-methods --------------------------------
# Validation methods of the objects
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

# cvErrorRate
setValidity("cvErrorRate", function(object){
    test <- test.nonLength0(object@cvErrorRate, "cvErrorRate")
    if (test != TRUE){
      return(test)
    }
    if (any(is.na(object@cvErrorRate))){
      return(paste("'cvErrorRate' must not contain missing values (NA)"))
    }
    test <- test.rateValue(object@cvErrorRate, "cvErrorRate", greaterThan0=FALSE)
    if (test != TRUE){
      return(test)
    }
    test <- test.nonLength0Vector(object@seErrorRate, "seErrorRate")
    if (test != TRUE){
      return(test)
    }
    if (any(is.na(object@seErrorRate))){
      return(paste("'seErrorRate' must not contain missing values (NA)"))
    }
    if (length(object@cvErrorRate) != length(object@seErrorRate)){
      return(paste("'cvErrorRate' and 'seErrorRate' must have the same length ",
                   "(", length(object@cvErrorRate), "!=", length(object@seErrorRate), ")",
                   sep=""))
    }
    test <- test.nonDim0NumericMatrix(object@classErrorRates, "classErrorRates")
    if (test != TRUE){
      return(test)
    }
    if (any(is.na(object@classErrorRates))){
      return(paste("'classErrorRates' must not contain missing values (NA)"))
    }
    test <- test.rateValue(object@classErrorRates, "classErrorRates", greaterThan0=FALSE)
    if (test != TRUE){
      return(test)
    }
    if (length(object@cvErrorRate) != dim(object@classErrorRates)[2]){
      return(paste("'classErrorRates' must have a number of columns equal to the length of 'cvErrorRate'",
                   "(", length(object@cvErrorRate), "!=", dim(object@classErrorRates)[2], ")",
                   sep=""))
    }
    return(TRUE)
  })

# errorRate1stLayerCV
setValidity("errorRate1stLayerCV", function(object) {
  test <- test.nonDim0NumericMatrix(object@errorRatePerFold, "errorRatePerFold")
  if (test != TRUE){
    return(test)
  }
  if (any(is.na(object@errorRatePerFold))){
    return(paste("'errorRatePerFold' must not contain missing values (NA)"))
  }
  test <- test.rateValue(object@errorRatePerFold, "errorRatePerFold", greaterThan0=FALSE)
  if (test != TRUE){
    return(test)
  }
  if (length(object@cvErrorRate@cvErrorRate) != dim(object@errorRatePerFold)[2]){
    return(paste("'errorRatePerFold' must provide a value per fold and per model (", length(object@cvErrorRate@cvErrorRate), "models )"))
  }
  return(TRUE)
  })
    
# ************ frequencyTopGenePerOneModel *************
setValidity("frequencyTopGenePerOneModel", function(object) {
    # Can be 0 with NSC
    #resTest <- test.nonLength0(object, "frequencyTopGenePerOneModel")
#    if ( resTest != TRUE ){
#      return(resTest)
#    }
    if (length(object) > 0){
        for (i in 1:length(object)){
          if ( class(object[[i]]) != "frequencyGenes" ){
            return(paste( "'frequencyTopGenePerOneModel' must be a list of element of class 'frequencyGenes'",
                              "(element", i, "is not)"))
          }
        }
    }
    return(TRUE)
  })
 
# *********** frequencyGenes ********************
setValidity("frequencyGenes", function(object) {  
    resTest <- test.rateValue(object@frequ, "frequ", greaterThan0=TRUE)
    if ( resTest != TRUE ){
      return(resTest)  
    }
    if (object@frequ == 0){
      return(paste("'frequ', must be >0 (instead of ",
            object@frequ,
            ")",
            sep=""))
    }
    if ( length(object@genesList) < 1){
      return("'genesList', must contains at least 1 gene (instead of 0)")
    }
    return(TRUE)
})
  
# ********** selectedGenes1stLayerCV ***********
setValidity("selectedGenes1stLayerCV", function(object) {
        for (i in 1:length(object@selectedGenesPerFold)){
          if (! class(object@selectedGenesPerFold[[i]]) == "selectedGenesPerOneOption"){
            return(paste( "selectedGenesPerFold must be a list of element of class 'selectedGenesPerOneOption'",
                          "(element", i, "is not)"))
          }
        }
        for (i in 1:length(object@frequencyTopGene)){
          if (! class(object@frequencyTopGene[[i]]) == "frequencyTopGenePerOneModel"){
            return( paste( "frequencyTopGene must be a list of elements of class 'frequencyTopGenePerOneModel'",
                           "(element", i, "is not)") )
          }
        }
        if (length(object@selectedGenesPerFold) != length(object@frequencyTopGene)){
          return(paste("selectedGenesPerFold and frequencyTopGene must have the same no of models (",
                        length(object@selectedGenesPerFold),
                        "!=",
                        length(object@frequencyTopGene),
                        ")",
                        sep="" ))
        }

      return(TRUE)
    })

# ********** resultSingle1LayerCV **************
setValidity("resultSingle1LayerCV", function(object) {
    # Can be = 0 if threshold in nsc  (and a float with nsc)
  test <- test.positiveFloat(object@bestOptionValue, "bestOptionValue", FALSE)
  if (test != TRUE){
    return(test)
  }
  return(TRUE)
})


# ********* resultRepeated1LayerCV *********
setValidity("resultRepeated1LayerCV", function(object) {
    for (i in 1:length(object@original1LayerCV)){
      if (! class(object@original1LayerCV[[i]]) == "resultSingle1LayerCV" ){  
        return( paste( "'original1LayerCV' must be a list of elements of class 'resultSingle1LayerCV'",
                           "(element", i, "is not)") )  
      }
    }
    # Can be a float Can be = 0 for nsc threshold
    test <- test.positiveFloat(object@bestOptionValue, "bestOptionValue", FALSE)
    if (test != TRUE){
      return(test)
    }
    return(TRUE) 
  })

# ********** errorRate2ndLayerCV *********
setValidity("errorRate2ndLayerCV", function(object) {
    test <- test.nonLength0(object@errorRatePerFold, "errorRatePerFold")
    if (test != TRUE){
        return(test)
    }
    if (any(is.na(object@errorRatePerFold))){
        return("'errorRatePerFold' must not contains missing values (NA)")
    }
    test <- test.rateValue(object@errorRatePerFold, "errorRatePerFold", greaterThan0=FALSE)
    if (test!=TRUE){
        return(test)
    }
    return(TRUE)
})

setValidity("cvErrorRate2ndLayer", function(object) {
  test <- test.oneElement(object@finalErrorRate, "finalErrorRate")
  if (test != TRUE){
    return(test)
  }
  if (is.na(object@finalErrorRate)){
    return("'finalErrorRate' must not be a missing value (NA)")
  }
  resTestRate <- test.rateValue(object@finalErrorRate, "finalErrorRate", greaterThan0=FALSE)
  if (resTestRate!= TRUE){
    return(resTestRate)
  }
  test <- test.oneElement(object@seFinalErrorRate, "seFinalErrorRate")
  if (test != TRUE){
    return(test)
  }
  if (is.na(object@seFinalErrorRate)){
    return("'seFinalErrorRate' must not be a missing value (NA)")
  }
  #if (length(object@errorRatePerFold) < 1){
#    return(paste("'errorRatePerFold', must be of length >=1 (instead of ",
#            length(object@errorRatePerFold),
#            ")",
#            sep=""))
#  }
  test <- test.nonLength0(object@classErrorRates, "classErrorRates")
  if (test != TRUE){
    return(test)
  }
  if (any(is.na(object@classErrorRates))){
    return(paste("'classErrorRates' must not contain missing values (NA)"))
  }
  test <- test.rateValue(object@classErrorRates, "classErrorRates", greaterThan0=FALSE)
  if (test!=TRUE){
    return(test)
  }
  
  return(TRUE)
})

# ********* selectedGenes ********
setValidity("selectedGenes", function(object) {
  resPosInt <- test.positiveInteger(object@noOfGenes, "noOfGenes", TRUE)
  if (resPosInt != TRUE){
    return(resPosInt)  
  }
  if ( length(object@genesList) < 1){
      return("'genesList', must contains at least 1 gene (instead of 0)")
  }
  if ( length(object@genesList)!=object@noOfGenes ) {
    return(paste( "'genesList', must contains 'noOfGenes' genes ",
                  "(", object@noOfGenes, " instead of ",
                  length(object@genesList), ")",
                  sep=""))
  }
  return(TRUE)
})


# ********* selectedGenes2ndLayerCV *****************
setValidity("selectedGenes2ndLayerCV", function(object) {
    if (length(object)==0){
      return(paste("'selectedGenes2ndLayerCV', must be of length >=1 (instead of ",
            length(object),
            ")",
            sep=""))  
    }
    for (i in 1:length(object)){
      if (class(object[[i]]) != "selectedGenes"){
        return( paste( "'selectedGenes2ndLayerCV' must be a list of elements of class 'selectedGenes'",
                             "(element", i, "is not)") )
      }
    }
  })

## ********** result2LayerCV *************
setValidity("result2LayerCV", function(object) {
  if (! object@avgBestOptionValue >= 0){
    return(paste("'avgBestOptionValue', must be >=0 (instead of ",
            object@avgBestOptionValue,
            ")",
            sep=""))  
  }
  if (length(object@results1stLayer)==0){
      return(paste("'selectedGenes2ndLayerCV', must be of length >=1 (instead of ",
            length(object@results1stLayer),
            ")",
            sep=""))  
    }
    for (i in 1:length(object@results1stLayer)){
      if (class(object@results1stLayer[[i]]) != "resultRepeated1LayerCV"){
        return( paste( "'results1stLayer' must be a list of elements of class 'resultRepeated1LayerCV'",
                             "(element", i, "is not)") )
      }
    }
    return(TRUE)
})

## ********** dataset *************
#setValidity("dataset", function(object) {
#    test <- test.oneElement(object@dataId, "dataId")
#    if (test != TRUE){
#    return(test)
#    }
#    test <- test.oneElement(object@dataPath, "dataPath")
#    if (test != TRUE){
#    return(test)
#    }
#    test <- test.oneElement(object@geneExprFile, "geneExprFile")
#    if (test != TRUE){
#    return(test)
#    }
#    test <- test.oneElement(object@classesFile, "classesFile")
#    if (test != TRUE){
#    return(test)
#    }
#    # Path to the gene expression file must be accessible
#    pathGeneExpr <- file.path(object@dataPath, object@geneExprFile)
#    if (file.access(pathGeneExpr, 4) == -1){
#    return(paste("'geneExprFile' inexistant or inaccessible (",
#                  pathGeneExpr,
#                  ")",
#                  sep="" ))
#    }
#    # Path to the classes file must be accessible
#    pathClasses <- file.path(object@dataPath, object@classesFile)
#    if (file.access(pathClasses, 4) == -1){
#    return(paste("'classesFile' inexistant or inaccessible (",
#                  pathClasses,
#                  ")",
#                  sep="" ))
#    }
#})

# ********** geneSubsets *********
setValidity("geneSubsets", function(object) {
    possibleSpeed <- c("high", "slow")
    position <- pmatch(object@speed, possibleSpeed, nomatch=0)
    if (position == 0 || possibleSpeed[position] != object@speed){
        return(paste("'speed' must be 'high' or 'slow' (instead of ", object@speed, ")", sep=""))
    }
    test <- test.nonLength0Vector(object@optionValues, "optionValues")
    if (test != TRUE){
        return(test)
    }
    test <- test.oneElement(object@maxSubsetSize, "maxSubsetSize")
    if (test != TRUE){
        return(test)
    }
    test <- test.oneElement(object@noOfOptions, "noOfOptions")
    if (test != TRUE){
        return(test)
    }
    test <- test.orderedPositiveInteger(object@optionValues, "optionValues", TRUE)
    if (test != TRUE){
        return(test)
    }
    test <- test.positiveInteger(object@maxSubsetSize, "maxSubsetSize", TRUE)
    if (test != TRUE){
        return(test)
    }
    test <- test.positiveInteger(object@noOfOptions, "noOfOptions", TRUE)
    if (test != TRUE){
        return(test)
    }
    return(TRUE)
})

# ********** thresholds *********
setValidity("thresholds", function(object) {
    test <- test.nonLength0Vector(object@optionValues, "optionValues")
    if (test != TRUE){
        return(test)
    }
    test <- test.oneElement(object@noOfOptions, "noOfOptions")
    if (test != TRUE){
        return(test)
    }
    test <- test.orderedPositiveFloat(object@optionValues, "optionValues", FALSE)
    if (test != TRUE){
        return(test)
    }
    test <- test.positiveInteger(object@noOfOptions, "noOfOptions", TRUE)
    if (test != TRUE){
        return(test)
    }
    return(TRUE)
})

# ********** assessment **********
setValidity("assessment", function(object) {
#    if ( is.null(object@dataset@eset) ){
#        return("dataset must contain a non-NULL eset attributes, please use loadData on your dataset before assigning it to the assessment")
#    }
    test <- test.oneElement(object@noFolds1stLayer, "noFolds1stLayer")
    if (test != TRUE){
        return(test)
    }
    test <- test.positiveInteger(object@noFolds1stLayer, "noFolds1stLayer", TRUE)
    if (test != TRUE){
        return(test)
    }
    if (object@noFolds1stLayer < 1){
        return(paste("'noFolds1stLayer' must be >=1 (instead of ", object@noFolds1stLayer,")", sep=""))
    }
    test <- test.oneElement(object@noFolds2ndLayer, "noFolds2ndLayer")
        if (test != TRUE){
        return(test)
    }
    test <- test.positiveInteger(object@noFolds2ndLayer, "noFolds2ndLayer", TRUE)
    if (test != TRUE){
        return(test)
    }
    if (object@noFolds2ndLayer < 1){
        return(paste("'noFolds2ndLayer' must be >=1 (instead of ", object@noFolds2ndLayer,")", sep=""))
    }
    possibleFoldCreation <- c("original", "balanced", "naive")
    position <- pmatch(object@typeFoldCreation, possibleFoldCreation, nomatch=0)
    if (position == 0 || possibleFoldCreation[position] != object@typeFoldCreation){
        return(paste("typeFoldCreation must be 'original', 'balanced' or 'naive' (instead of ", object@typeFoldCreation, ")", sep=""))
    }
#  if (object@classifierName == "svm"){
#   Always since the RFE is based on SVM
    possibleKernel <- c("linear", "polynomial", "radial")
    position <- pmatch(object@svmKernel, possibleKernel, nomatch=0)
    if (position == 0 || possibleKernel[position] != object@svmKernel){
        return(paste("svmKernel must be 'linear', 'polynomial' or 'radial' (instead of ", object@svmKernel, ")", sep=""))
    }
  #}
    # Code useful for lda or naiveBayes, unused in current version ################
    ###############################################################################
    #possibleClassifier <- c("svm", "lda", "naiveBayes", "nsc")
    # Code useful for lda or naiveBayes, unused in current version ################
    ###############################################################################
    possibleClassifier <- c("svm", "nsc")
    position <- pmatch(object@classifierName, possibleClassifier, nomatch=0)
    if (position == 0 || possibleClassifier[position] != object@classifierName){
        return(paste("classifierName must be ", paste(possibleClassifier, collapse=", ")," (instead of ", object@classifierName, ")", sep=""))
    }
    possibleFeatSelection <- c("rfe", "nsc")
    position <- pmatch(object@featureSelectionMethod, possibleFeatSelection, nomatch=0)
        if (position == 0 || possibleFeatSelection[position] != object@featureSelectionMethod){
            return(paste("featureSelectionMethod must be ", paste(possibleFeatSelection, collapse=", ")," (instead of ", object@featureSelectionMethod, ")", sep=""))
    }
    if (object@classifierName == 'svm' && object@featureSelectionMethod != 'rfe'){
        return("If 'classifierName' is 'svm' then 'featuresSelectionMethod' must be 'rfe'.")
    }
    if (object@classifierName == 'nsc' && object@featureSelectionMethod != 'nsc'){
        return("If 'classifierName' is 'nsc' then 'featuresSelectionMethod' must be 'nsc' too.")
    }
    if (object@classifierName != 'nsc' && object@featureSelectionMethod == 'nsc'){
        return("In the current implementation, nsc feature selection method ('featuresSelectionMethod') is restricted to the use of nsc as a classifier ('classifierName').")
    }
    test <- test.positiveInteger(object@noOfRepeats, "noOfRepeats", TRUE)
    if (test != TRUE){
        return(test)
    }
    if (class(object@featureSelectionOptions) == 'geneSubsets') {
        if (object@featureSelectionOptions@maxSubsetSize > length(featureNames(object@dataset)) ){
            return(paste("The maximum of genes in 'geneSubsets'(",
                      object@featureSelectionOptions@maxSubsetSize,
                      ") must not be greater than the number of features in 'dataset'(",
                      length(featureNames(object@dataset)),
                      ")",
                      sep=""))
        }
    }
    if (object@featureSelectionMethod == 'nsc' && class(object@featureSelectionOptions) != 'thresholds'){
        return("If 'featureSelectionMethod' is 'nsc' then 'subsetsType' in 'featureSelectionOptions' must be of class 'thresholds'.")
    }
    if (object@featureSelectionMethod == 'rfe' && class(object@featureSelectionOptions) != 'geneSubsets'){
        return("If 'featureSelectionMethod' is 'rfe' then 'subsetsType' in 'featureSelectionOptions' must be of class 'geneSubsets'.")
    }
  return(TRUE)
})