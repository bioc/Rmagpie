# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#-------------------------- resultRepeated1LayerCV-accessors -------------------
# Class resultRepeated1LayerCV Accessors and replacement methods
#
# Object storing the result of a repeated one layer external CV
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#---------------------------- Getters ------------------------------------------
# frequencyGenes
setMethod("getFrequency", "frequencyGenes", function(object) object@frequ)
setMethod("getGenesList", "frequencyGenes", function(object) object@genesList)

# errorRate1stLayerCV
setMethod("getErrorRatePerFold", "errorRate1stLayerCV", function(object) object@errorRatePerFold)
setMethod("getNoSamplesPerFold", "errorRate1stLayerCV", function(object) object@noSamplesPerFold)
setMethod("getCvErrorRate", "errorRate1stLayerCV", function(object) object@cvErrorRate)

# cvErrorRate
setMethod("getCvErrorRate", "cvErrorRate", function(object) object@cvErrorRate)
setMethod("getSeErrorRate", "cvErrorRate", function(object) object@seErrorRate)
setMethod("getClassErrorRates", "cvErrorRate", function(object) object@classErrorRates)

# selectedGenesPerOneOption
setMethod("getOptionValue", "selectedGenesPerOneOption", function(object) object@optionValue)
setMethod("getGenesList", "selectedGenesPerOneOption", function(object) object)

# selectedGenes1stLayerCV
setMethod("getSelectedGenesPerFold", "selectedGenes1stLayerCV", function(object) object@selectedGenesPerFold)
setMethod("getFrequencyTopGene", "selectedGenes1stLayerCV", function(object) object@frequencyTopGene)

# resultSingle1LayerCV
setMethod("getErrorRates", "resultSingle1LayerCV", function(object) object@errorRates)
setMethod("getSelectedGenes", "resultSingle1LayerCV", function(object) object@selectedGenes)
setMethod("getExecutionTime", "resultSingle1LayerCV", function(object) object@executionTime)
setMethod("getBestOptionValue", "resultSingle1LayerCV", function(object) object@bestOptionValue)

# resultRepeated1LayerCV
#setMethod("getErrorRates", "resultRepeated1LayerCV", function(object) object@summaryErrorRate@cvErrorRate)
setMethod("getOriginal1LayerCV", "resultRepeated1LayerCV", function(object) object@original1LayerCV)
setMethod("getSummaryErrorRate", "resultRepeated1LayerCV", function(object) object@summaryErrorRate)
setMethod("getBestOptionValue", "resultRepeated1LayerCV", function(object) object@bestOptionValue)
setMethod("getFrequencyTopGenes", "resultRepeated1LayerCV", function(object) object@summaryFrequencyTopGenes)
setMethod("getExecutionTime", "resultRepeated1LayerCV", function(object) object@executionTime)




              


            
            
              
