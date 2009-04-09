# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#-------------------------- thresholds-accessors ----------------------------------
# Class thresholds Accessors and replacement methods
#
# Object storing the thresholds to be used by NSC
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#---------------------------- Getters ------------------------------------------
setMethod("getNoThresholds", "thresholds", function(object) object@noOfOptions)
setMethod("getOptionValues", "thresholds", function(object) object@optionValues)

#---------------------------- Setters ------------------------------------------
setReplaceMethod("getOptionValues", "thresholds",
    function(object, value) {
        object@optionValues <- value
        object@noOfOptions <- length(object@optionValues)
        object
        }
)