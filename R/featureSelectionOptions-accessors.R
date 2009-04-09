# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#------------------- featureSelectionOptions-accessors -------------------------
# Class featureSelectionOptions Accessors and replacement methods
#
# Object storing the options corresponding to an assessment.
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#---------------------------- Getters ------------------------------------------
setMethod("getOptionValues", "featureSelectionOptions",
                function(object) object@optionValues)
setMethod("getNoOfOptions", "featureSelectionOptions",
                function(object) object@noOfOptions)