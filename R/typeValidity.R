# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#--------------------------- typeValidity --------------------------------------
# Helpers to check the validity of the objects
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Test if this input matrix contains only rate values (>=0 and <=1)
#
# @param    (matrix) Matrix to be checked
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
setMethod("test.rateValue", "matrix", function(attValue, attName, greaterThan0=FALSE) {
    nbRows <- nrow(attValue)
    nbCols <- ncol(attValue)
    for (i in 1:nbRows){
      for (j in 1:nbCols){
        elementRateValue <- test.rateValue( attValue[i,j],
                                            paste(attName, "[", i, ",", j, "]", sep=""),
                                            greaterThan0)
        if (elementRateValue != TRUE)  {
          return(elementRateValue)
        }
      }
    }
    return(TRUE)
  })

#-------------------------------------------------------------------------------
# Test if this input list contains only rate values (>=0 and <=1)
#
# @param    (list) List to be checked
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
setMethod("test.rateValue", "list", function(attValue, attName, greaterThan0=FALSE) {
    lengthList <- length(attValue)
    for (i in 1:lengthList){
      elementRateValue <- test.rateValue( attValue[[i]],
                                          paste(attName, "[[", i, "]]", sep=""),
                                          greaterThan0)
      if (elementRateValue != TRUE)  {
        return(elementRateValue)
      }
    }
    return(TRUE)
  })

#-------------------------------------------------------------------------------
# Test if this input numeric is a rate value (>=0 and <=1)
#
# @param    (numeric) Value to be checked
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
setMethod("test.rateValue", "numeric", function(attValue, attName, greaterThan0=FALSE) {
    nbElements <- length(attValue)
    
    if (nbElements == 1){
      if (greaterThan0){
        condition1 <- attValue > 0
        textCondition1 <- ">0"
      } else {
        condition1 <- attValue >= 0
        textCondition1 <- ">=0"
      }
      if (!condition1 || !attValue <= 1 ){
        return(paste( "'", attName, "', must be ",
                      textCondition1," and <=1 (instead of ",
                      attValue, ")",
                      sep=""))
      } else {
        return(TRUE)
      }
    } else {
      for (i in 1:nbElements){
        if (greaterThan0){
          condition1 <- attValue[i] > 0
          textCondition1 <- ">0"
        } else {
          condition1 <- attValue[i] >= 0
          textCondition1 <- ">=0"
        }
        if (!condition1 || !attValue[i] <= 1 ){
        return(paste( "'", attName, "[", i, "]' ",
                      "must be ", textCondition1, " and <=1 (instead of ",
                      paste(attValue[i], collapse=" "),
                      ")",
                      sep=""))
        }
      }
    }
    return(TRUE)
  })

#-------------------------------------------------------------------------------
# Test if this input matrix contains only integers
#
# @param    (matrix) Matrix to be checked
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
setMethod("test.integer", "matrix", function(attValue, attName) {
    nbRows <- nrow(attValue)
    nbCols <- ncol(attValue)
    for (i in 1:nbRows){
      for (j in 1:nbCols){
        elementRateValue <- test.integer( attValue[i,j],
                                            paste(attName, "[", i, ",", j, "]", sep=""))
        if (elementRateValue != TRUE)  {
          return(elementRateValue)
        }
      }
    }
    return(TRUE)
  })

#-------------------------------------------------------------------------------
# Test if this input list contains only integers
#
# @param    (list) List to be checked
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
setMethod("test.integer", "list", function(attValue, attName) {
    lengthList <- length(attValue)
    for (i in 1:lengthList){
      elementRateValue <- test.integer( attValue[[i]],
                                          paste(attName, "[[", i, "]]", sep=""))
      if (elementRateValue != TRUE)  {
        return(elementRateValue)
      }
    }
    return(TRUE)
  })

#-------------------------------------------------------------------------------
# Test if this input numeric is an integer
#
# @param    (numeric) Value to be checked
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
setMethod("test.integer", "numeric", function(attValue, attName) {
    nbElements <- length(attValue)
    if (nbElements == 1){
      if ( round(attValue) != attValue ){
        return(paste("'", attName, "' must be an integer (instead of ",
            attValue,
            ")",
            sep=""))
      } else {
        return(TRUE)
      }
    } else {
      for (i in 1:nbElements){
        if ( round(attValue[i]) != attValue[i] ){
          return(paste("'", attName, "[",i,"]' ",
            "must be an integer (instead of ",
            attValue[i],
            ")",
            sep=""))
        }
      }
    }
    return(TRUE)
  })

#-------------------------------------------------------------------------------
# Test if this input matrix contains only positive integers
#
# @param    (matrix) Matrix to be checked
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
setMethod("test.positiveInteger", "matrix", function(attValue, attName, greaterThan0=FALSE) {
    nbRows <- nrow(attValue)
    nbCols <- ncol(attValue)
    for (i in 1:nbRows){
      for (j in 1:nbCols){
        elementRateValue <- test.positiveInteger( attValue[i,j],
                                            paste(attName, "[", i, ",", j, "]", sep=""),
                                            greaterThan0)
        if (elementRateValue != TRUE)  {
          return(elementRateValue)
        }
      }
    }
    return(TRUE)
  })

#-------------------------------------------------------------------------------
# Test if this input list contains only positive integers
#
# @param    (matrix) List to be checked
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
setMethod("test.positiveInteger", "list", function(attValue, attName, greaterThan0=FALSE) {
    lengthList <- length(attValue)
    for (i in 1:lengthList){
      elementRateValue <- test.positiveInteger( attValue[[i]],
                                          paste(attName, "[[", i, "]]", sep=""),
                                          greaterThan0)
      if (elementRateValue != TRUE)  {
        return(elementRateValue)
      }
    }
    return(TRUE)
  })

#-------------------------------------------------------------------------------
# Test if this input numeric is a positive integer
#
# @param    (numeric) Value to be checked
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
setMethod("test.positiveInteger", "numeric", function(attValue, attName, greaterThan0=FALSE) {
    nbElements <- length(attValue)
    if (nbElements == 1){
      test <- test.integer(attValue, attName)
      if (test != TRUE){
        return(test)
      }
      if (greaterThan0) {
        condition <- attValue > 0
        conditionText <- ">0"
      } else {
        condition <- attValue >= 0
        conditionText <- ">=0"
      }
      if (! condition){
        return( paste("'", attName, "' ",
                "must be ", conditionText,
                " (instead of ",
                attValue,
                ")",
                sep=""))
      }
      return(TRUE)
    } else {
      for (i in 1:nbElements){
        test <- test.integer(attValue[i], paste(attName, "[", i, "]", sep=""))
        if (test != TRUE){
          return(test)
        }
        if (greaterThan0) {
          condition <- attValue[i] > 0
          conditionText <- ">0"
        } else {
          condition <- attValue[i] >= 0
          conditionText <- ">=0"
        }
        if (! condition){
          return( paste("'", attName, "[", i, "]' ",
                  "must be ", conditionText,
                  " (instead of ",
                  attValue[i],
                  ")",
                  sep=""))
        }
      }
    }
    return(TRUE)
  })

#-------------------------------------------------------------------------------
# Test if this input matrix contains only positive floats
#
# @param    (matrix) Matrix to be checked
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
setMethod("test.positiveFloat", "numeric", function(attValue, attName, greaterThan0=FALSE) {
    nbElements <- length(attValue)
    if (nbElements == 1){
      if (greaterThan0) {
        condition <- attValue > 0
        conditionText <- ">0"
      } else {
        condition <- attValue >= 0
        conditionText <- ">=0"
      }
      if (! condition){
        return( paste("'", attName, "' ",
                "must be ", conditionText,
                " (instead of ",
                attValue,
                ")",
                sep=""))
      }
      return(TRUE)
    } else {
      for (i in 1:nbElements){
        if (greaterThan0) {
          condition <- attValue[i] > 0
          conditionText <- ">0"
        } else {
          condition <- attValue[i] >= 0
          conditionText <- ">=0"
        }
        if (! condition){
          return( paste("'", attName, "[", i, "]' ",
                  "must be ", conditionText,
                  " (instead of ",
                  attValue[i],
                  ")",
                  sep=""))
        }
      }
    }
    return(TRUE)
  })

#-------------------------------------------------------------------------------
# Test if this input list contains only positive intergers in ascending order
#
# @param    (list) List to be checked
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
setMethod("test.orderedPositiveInteger", "list", function(attValue, attName, greaterThan0=FALSE) {
    lengthList <- length(attValue)
    for (i in 1:lengthList){
      elementRateValue <- test.orderedPositiveInteger( attValue[[i]],
                                          paste(attName, "[[", i, "]]", sep=""),
                                          greaterThan0)
      if (elementRateValue != TRUE)  {
        return(elementRateValue)
      }
    }
    return(TRUE)
  })

#-------------------------------------------------------------------------------
# Test if this input numeric contains only positive intergers in ascending order
#
# @param    (numeric) Value or vector to be checked
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
setMethod("test.orderedPositiveInteger", "numeric", function(attValue, attName, greaterThan0=FALSE) {
    nbElements <- length(attValue)
    previousValue <- 0
    for (i in 1:nbElements){
      test <- test.integer(attValue[i], paste(attName, "[", i, "]", sep=""))
      if (test != TRUE){
        return(test)
      }
      if (greaterThan0) {
        condition <- attValue[i] > 0
        conditionText <- ">0"
      } else {
        condition <- attValue[i] >= 0
        conditionText <- ">=0"
      }
      if (! condition){
        return( paste("'", attName, "[", i, "]' ",
                "must be ", conditionText,
                " (instead of ",
                attValue[i],
                ")",
                sep=""))
      }
      if (attValue[i] < previousValue){
        return( paste("'", attName, "' must be in croissant order",
                " (",attName,"[",
                i,
                "])",
                sep=""))
      }
      previousValue <- attValue[i]
    }
    return(TRUE)
  })

#-------------------------------------------------------------------------------
# Test if this input numeric contains only positive floats in ascending order
#
# @param    (numeric) Value or vector to be checked
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
setMethod("test.orderedPositiveFloat", "numeric", function(attValue, attName, greaterThan0=FALSE) {
    nbElements <- length(attValue)
    previousValue <- 0
    for (i in 1:nbElements){
      if (greaterThan0) {
        condition <- attValue[i] > 0
        conditionText <- ">0"
      } else {
        condition <- attValue[i] >= 0
        conditionText <- ">=0"
      }
      if (! condition){
        return( paste("'", attName, "[", i, "]' ",
                "must be ", conditionText,
                " (instead of ",
                attValue[i],
                ")",
                sep=""))
      }
      if (attValue[i] < previousValue){
        return( paste("'", attName, "' must be in croissant order",
                " (",attName,"[",
                i,
                "])",
                sep=""))
      }
      previousValue <- attValue[i]
    }
    return(TRUE)
  })

#-------------------------------------------------------------------------------
# Test if the input <attValue> named <attName> is not of length 0
#
# @param    attValue Value to be checked
#           attName (character) nems of the object
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
test.nonLength0Vector <- function(attValue, attName) {
  if (! is.vector(attValue)){
    return( paste( "'",
                   attName,
                   "' must be a vector (instead of ",
                   class(attValue),
                   ")",
                   sep="" ) )
  }
  if (length(attValue) == 0){
    return( paste( "'",
                   attName,
                   "' must be a vector of length > 0 (instead of ",
                   length(attValue),
                   ")",
                   sep="" ) )
  }
  return(TRUE)
}

#-------------------------------------------------------------------------------
# Test if the input <attValue> named <attName> is not of length 0
#
# @param    attValue Value to be checked
#           attName (character) nems of the object
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
test.nonLength0 <- function(attValue, attName) {
  if (length(attValue) == 0){
    return( paste( "'",
                   attName,
                   "' must be of length > 0 (instead of ",
                   length(attValue),
                   ")",
                   sep="" ) )
  }
  return(TRUE)
}

#-------------------------------------------------------------------------------
# Test if the matrix <attValue> named <attName> does not have a 0 dimension
#
# @param    attValue Matrix to be checked
#           attName (character) nems of the object
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
test.nonDim0NumericMatrix <- function(attValue, attName) {
  if (! is.numeric(attValue) ){
    return( paste( "'",
                   attName,
                   "' must be a numeric matrix (instead of ",
                   mode(attValue),
                   " matrix )",
                   sep="" ) )
  }
  if (any(dim(attValue) == 0)){
    return( paste( "'",
                   attName,
                   "' must be a matrix of non nul dimensions (instead of ",
                   paste(dim(attValue), collapse=" "),
                   ")",
                   sep="" ) )
  }
  return(TRUE)
}

#-------------------------------------------------------------------------------
# Test that the input <attValue> has length 1
#
# @param    attValue object to be checked
#           attName (character) nems of the object
# @return   (logical) TRUE if the object passed the ckecking
#           (character) error message if the object did not pass the checking
#-------------------------------------------------------------------------------
test.oneElement <- function(attValue, attName) {
  if (! is.vector(attValue) ){
    return( paste( "'",
                   attName,
                   "' must be a vector of length 1 (instead of ",
                   class(attValue),
                   ")",
                   sep="" ) )
  }
  if (length(attValue) != 1){
    return( paste( "'",
                   attName,
                   "' must be a vector of length 1 (instead of ",
                   length(attValue),
                   ")",
                   sep="" ) )
  }
  return(TRUE)
}