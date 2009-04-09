# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#--------------------------- finalClassifier-accessors -------------------------
# Class finalClassifier Accessors and replacement methods
#
# Object storing the final classifier corresponding to an exepriment.
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#---------------------------- Getters ------------------------------------------
setMethod("getGenesFromBestToWorst", "finalClassifier", function(object) object@genesFromBestToWorst)
setMethod("getModels", "finalClassifier", function(object) object@models)

#---------------------------- Setters ------------------------------------------
setReplaceMethod("getGenesFromBestToWorst", "finalClassifier",
    function(object, value) {
        object@genesFromBestToWorst <- value
        object
    }
)

setReplaceMethod("getModels", "finalClassifier",
    function(object, value) {
        object@models <- value
        object
    }
)