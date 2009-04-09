# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#-------------------------- result2LayerCV-accessors ---------------------------
# Class result2LayerCV Accessors and replacement methods
#
# Object storing the result of a 2 layer external CV
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#---------------------------- Getters ------------------------------------------
# cvErrorRate2ndLayer
setMethod("getFinalErrorRate", "cvErrorRate2ndLayer", function(object) object@finalErrorRate)
setMethod("getSeErrorRate", "cvErrorRate2ndLayer", function(object) object@seFinalErrorRate)
setMethod("getClassErrorRates", "cvErrorRate2ndLayer", function(object) object@classErrorRates)

# errorRate2ndLayerCV
setMethod("getErrorRatePerFold", "errorRate2ndLayerCV", function(object) object@errorRatePerFold)
setMethod("getNoSamplesPerFold", "errorRate2ndLayerCV", function(object) object@noSamplesPerFold)
setMethod("getCvErrorRate", "errorRate2ndLayerCV", function(object) object@cvErrorRate)

# selectedGenes
setMethod("getNoOfGenes", "selectedGenes", function(object) object@noOfGenes)
setMethod("getOptionValue", "selectedGenes", function(object) object@optionValue)
setMethod("getGenesList", "selectedGenes", function(object) object@genesList)

# result2LayerCV
setMethod("getResults1stLayer", "result2LayerCV", function(object) object@results1stLayer)
setMethod("getErrorRates", "result2LayerCV", function(object) object@errorRates)
setMethod("getSelectedGenes", "result2LayerCV", function(object) object@selectedGenes)
setMethod("getExecutionTime", "result2LayerCV", function(object) object@executionTime)
setMethod("getAvgBestOptionValue", "result2LayerCV", function(object) object@avgBestOptionValue)

# resultRepeated2LayerCV
setMethod("getOriginal2LayerCV", "resultRepeated2LayerCV", function(object) object@original2LayerCV)
setMethod("getAvgBestOptionValue", "resultRepeated2LayerCV", function(object) object@avgBestOptionValue)
setMethod("getSummaryErrorRate", "resultRepeated2LayerCV", function(object) object@summaryErrorRate)
setMethod("getExecutionTime", "resultRepeated2LayerCV", function(object) object@executionTime)