# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#--------------------------- geneSubsets-accessors -------------------------
# Class geneSubsets Accessors and replacement methods
#
# Object storing the sizes of subsets to be tried during the RFE.
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#---------------------------- Getters ------------------------------------------
setMethod("getMaxSubsetSize", "geneSubsets", function(object) object@maxSubsetSize)
setMethod("getNoModels", "geneSubsets", function(object) object@noOfOptions)
setMethod("getSubsetsSizes", "geneSubsets", function(object) object@optionValues)
setMethod("getSpeed", "geneSubsets", function(object) object@speed)

#---------------------------- Setters ------------------------------------------
setReplaceMethod("getMaxSubsetSize", "geneSubsets",
    function(object, value) {
        object@maxSubsetSize <- value
        if (object@speed == "high"){
          noOfOptions <- ceiling(log2(object@maxSubsetSize)+1)
          # Number of model (subset of genes) to be tested
          object@noOfOptions <- noOfOptions
          # noSelectedFeature contains the different number of genes that we will try
          object@optionValues <- 2^(0:(noOfOptions-1))
          object@optionValues[noOfOptions] <- object@maxSubsetSize
        } else{
          object@noOfOptions <- object@maxSubsetSize
          object@optionValues <- seq(from = 1, to = object@maxSubsetSize)
        }
        object
    }
)

setReplaceMethod("getSubsetsSizes", "geneSubsets",
    function(object, value) {
        object@optionValues <- value
        if (length(object@optionValues) < max(object@optionValues)){
            object@speed <- "high"
        } else {
           object@speed <- "slow"
        }
        object@noOfOptions <- length(object@optionValues)
        object
        }
)

setReplaceMethod("getSpeed", "geneSubsets",
    function(object, value) {
        object@speed <- value
        if (object@speed == "high"){
          noOfOptions <- ceiling(log2(object@maxSubsetSize)+1)
          # Number of model (subset of genes) to be tested
          object@noOfOptions <- noOfOptions
          # noSelectedFeature contains the different number of genes that we will try
          object@optionValues <- 2^(0:(noOfOptions-1))
          object@optionValues[noOfOptions] <- object@maxSubsetSize
        } else{
          object@noOfOptions <- object@maxSubsetSize
          object@optionValues <- seq(from = 1, to = object@maxSubsetSize)
        }
        object
        }
)
