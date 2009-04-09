# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#--------------------------- show-methods --------------------------------
# Show methods of the objects
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

# geneSubsets
setMethod("show", "geneSubsets", function(object) showWithPrefix(object))
setMethod("showWithPrefix", "geneSubsets", function(object, prefix="", short=FALSE) {
              cat( prefix,
                   "optionValues: ",
                   paste(round(getSubsetsSizes(object), digits=getOption("digits")),
                          collapse=" "),
                   sep="")
              cat( "  (maxSubsetSize: ",
                   round(getMaxSubsetSize(object), digits=getOption("digits")),
                   ", speed:",
                   getSpeed(object),
                   ", noOfOptions:",
                   getNoModels(object),
                   ")\n",
                   sep="")
              })

# thresholds
setMethod("show", "thresholds", function(object) showWithPrefix(object))
setMethod("showWithPrefix", "thresholds", function(object, prefix, short=FALSE) {
              cat( prefix,
                    "optionValues:")
                    if (length(getOptionValues(object)) > 0){
                        cat(round(getOptionValues(object), digits=getOption("digits")),
                                collapse=" ")
                    } else {
                        cat("EMPTY, can be used in an assessment")
                    }
              cat( "( noOfOptions:")
                    if (length(getNoThresholds(object)) > 0){
                        cat(getNoThresholds(object),
                                collapse=" ")
                    } else {
                        cat("EMPTY, can be used in an assessment")
                    }
               cat( ")\n",
                   sep="")
              })

# cvErrorRate
setMethod("show", "cvErrorRate", function(object) showWithPrefix(object) )
setMethod("showWithPrefix", "cvErrorRate", function(object, prefix, short=FALSE) {
    cat( prefix,
         "cvErrorRate:\n",
         prefix,
         "  ",
         paste( round( getCvErrorRate(object),
                       digits=getOption("digits")),
                       collapse=" "),
         "\n",
         prefix,
         "seErrorRate:\n",
         prefix,
         "  ",
         paste( round( getSeErrorRate(object),
                              digits=getOption("digits")),
                              collapse=" "),
        "\n",
        prefix,
        "classErrorRates:\n",
        sep="")
    noRows <- dim(getClassErrorRates(object))[1]
    classes <- rownames(getClassErrorRates(object))
    for (i in 1:noRows){
       cat( prefix,
            "  ",
            classes[i],
            ": ", 
            paste(round(getClassErrorRates(object)[i,], digits=getOption("digits")), collapse=" "),
            "\n",
            sep="")
    }         
  })

# errorRate1stLayerCV
setMethod("show", "errorRate1stLayerCV", function(object) showWithPrefix(object) )
setMethod("showWithPrefix", "errorRate1stLayerCV", function(object, prefix, short=FALSE) {
     cat( prefix,
          "errorRatePerFold: (",
          paste( dim(getErrorRatePerFold(object)),
                 c("folds", "subsets)"),
                 collapse=" and "),
          sep="")
    if (! short) {
        cat("\n")
        show(getErrorRatePerFold(object))
    } else {
        cat(": use 'getErrorRatePerFold(object)'\n")
    }
    cat( prefix,
         "noSamplesPerFold: ", sep="")
    if (! short) {
        cat( paste( getNoSamplesPerFold(object), collapse=", "),
             "\n",
             sep="")
    } else {
        cat("use 'getNoSamplesPerFold(object)'\n")
    }
    cat( prefix,
         "cvErrorRates: \n",
         sep="")
    showWithPrefix(getCvErrorRate(object), paste(prefix, "  ", sep=""))
    })

# frequencyGenes
setMethod("show", "frequencyGenes", function(object) showWithPrefix(object) )
setMethod("showWithPrefix", "frequencyGenes", function(object, prefix, short=FALSE) {
             cat(prefix, "frequency: ", getFrequency(object), "\n", sep="")
             cat(prefix, "genesList: ", paste(getGenesList(object), collapse=" "), "\n", sep="")
          })

# frequencyTopGenePerOneModel
setMethod("show", "frequencyTopGenePerOneModel", function(object) showWithPrefix(object) )
setMethod("showWithPrefix", "frequencyTopGenePerOneModel", function(object, prefix, short=FALSE) {
            if (length(object) >= 1){
                for (i in 1:length(object)){
                    showWithPrefix(object[[i]], paste(prefix, "  "), FALSE)
                }
            }
          })

# selectedGenesPerOneOption
setMethod("show", "selectedGenesPerOneOption", function(object) showWithPrefix(object) )
setMethod("showWithPrefix", "selectedGenesPerOneOption", function(object, prefix, short=FALSE) {
            if (! short){
                showDefault(object)
            } else {
                cat(prefix, "optionValue: ", getOptionValue(object), "\n",
                    prefix, "List of genes selected in ",
                    length(object),
                    " fold(s)\n",
                    sep="")
            }
            #show(getGenesList(object))
          })

# noModels <- length(object)
#                cat(  prefix,
#                      noModels,
#                      " values of option and ",
#                      ifelse(length(noModels)>=1, length(object[[1]]), 0),
#                      " fold(s)\n",
#                      sep="")

# TODO #3 revoir cette fonction
# selectedGenes1stLayerCV
setMethod("show", "selectedGenes1stLayerCV", function(object) showWithPrefix(object) )
setMethod("showWithPrefix", "selectedGenes1stLayerCV", function(object, prefix, short=FALSE) {
            noModels <- length(getSelectedGenesPerFold(object))
            cat( prefix, "selectedGenesPerFold: (",
                paste(  c(noModels,
                        ifelse( noModels>=1,
                                length(getSelectedGenesPerFold(object)[[1]]),
                                0)),
                        c("value(s) of option", "fold(s))"), collapse=" and "),
                sep="")
            if (! short) {
                cat("\n")
                show(getSelectedGenesPerFold(object))
            } else {
                cat(": use 'getSelectedGenesPerFold(object)'\n")
            }
            cat( prefix, "frequencyTopGene: (",
               length(getFrequencyTopGene(object)),
               " value(s) of option)",
               sep="")
            if (! short) {
                cat("\n")
                show(getFrequencyTopGene(object))
            } else {
                cat(": use 'getFrequencyTopGene(object)'\n")
            }
        })



# resultSingle1LayerCV
setMethod("show", "resultSingle1LayerCV", function(object) showWithPrefix(object) )
setMethod("showWithPrefix", "resultSingle1LayerCV", function(object, prefix, short=FALSE) {
              cat( prefix, "errorRates:\n", sep="")
              showWithPrefix(getErrorRates(object), paste(prefix,"  ", sep=""), TRUE)
              cat(prefix, "selectedGenes:\n", sep="")
              showWithPrefix(getSelectedGenes(object), paste(prefix,"  ", sep=""), TRUE)
              cat( prefix, "-> bestOptionValue:",
                   getBestOptionValue(object),
                   "\n", sep="")
              cat( prefix, "executionTime:",
                   round(getExecutionTime(object),
                        digits=getOption("digits")),
                   "s",
                   "\n", sep="")
              })


# resultRepeated1LayerCV
setMethod("show", "resultRepeated1LayerCV", function(object) showWithPrefix(object) )
setMethod("showWithPrefix", "resultRepeated1LayerCV",  function(object, prefix, short=FALSE) {
             cat( prefix,
                   "original1LayerCV: ",
                   length(getOriginal1LayerCV(object)),
                   "  combined 1 layer ",
                   dim(getResultsLayer(getOriginal1LayerCV(object)[[1]], topic='errorRate', errorType='fold'))[1],
                   "-folds CV",
                   sep="")
              if (! short){
                  cat("\n")
                  for (i in 1:length(getOriginal1LayerCV(object))){
                    cat( prefix, "  ",
                         "[[", i, "]] CV number ", i, "\n", sep="")
                    showWithPrefix(getOriginal1LayerCV(object)[[i]], paste(prefix,"  ", sep=""))
                  }
              } else {
                cat(": use 'getOriginal1LayerCV(object)'\n")
              }
              cat(prefix,
                  "summaryFrequencyTopGenes: use 'getFrequencyTopGenes(object)'",
                  sep="")
              cat("\n", prefix, "summaryErrorRate:\n", sep="")
              showWithPrefix(getSummaryErrorRate(object), paste(prefix,"  ", sep=""))
              cat(prefix, "=> bestOptionValue: ",
                getBestOptionValue(object),
                " (Est. Error rate: ",
                min(getResultsLayer(object, topic='errorRate', errorType='cv')),")\n",
                sep="")
              cat( prefix, "executionTime: ",
                getExecutionTime(object),
                "s\n",
                sep="")
              })


# cvErrorRate2ndLayer
setMethod("show", "cvErrorRate2ndLayer", function(object) showWithPrefix(object) )
setMethod("showWithPrefix", "cvErrorRate2ndLayer", function(object, prefix, short=FALSE) {
    cat( prefix,
         "finalErrorRate:\n",
         prefix,
         "  ",
         paste( round( getFinalErrorRate(object),
                       digits=getOption("digits")),
                       collapse=" "),
         "\n",
         prefix,
         "seFinalErrorRate:\n",
         prefix,
         "  ",
         paste( round( getSeErrorRate(object),
                              digits=getOption("digits")),
                              collapse=" "),
        "\n",
        prefix,
        "classErrorRates:\n",
        sep="")
    noClasses <- length(getClassErrorRates(object))[1]
    classes <- names(getClassErrorRates(object))
    for (i in 1:noClasses){
       cat( prefix,
            "  ",
            classes[i],
            ": ",
            paste(round(getClassErrorRates(object)[i], digits=getOption("digits")), collapse=" "),
            "\n",
            sep="")
    }
  })

# errorRate2ndLayerCV
setMethod("show", "errorRate2ndLayerCV", function(object) showWithPrefix(object))
setMethod("showWithPrefix", "errorRate2ndLayerCV", function(object, prefix="", short=FALSE){
        if (! short){
            cat(prefix, "errorRatePerFold:\n", prefix, "  ", sep="")
            for (i in 1:length(getErrorRatePerFold(object))){
              cat(round(getErrorRatePerFold(object)[[i]], digits=getOption("digits")), " ")
            }
            cat("\n", prefix, "noSamplesPerFold: ",
                paste(getNoSamplesPerFold(object), collapse=" "),
                "\n", sep="")
        } else {
            cat(prefix, "errorRatePerFold: use 'getErrorRatePerFold(object)'\n",
                prefix, "noSamplesPerFold: use 'getNoSamplesPerFold(object)'\n",
                sep="")
        }
        cat(prefix, "cvErrorRate: \n", sep="")
        showWithPrefix(getCvErrorRate(object), paste(prefix, "  ", sep=""), TRUE)
      })

# selectedGenes
setMethod("show", "selectedGenes", function(object) showWithPrefix(object))
setMethod("showWithPrefix", "selectedGenes", function(object, prefix="", short=FALSE){
        cat( prefix, "optionValue: ",
             getOptionValue(object),
             "\n", sep="")
        cat( prefix, "noOfGenes: ",
             getNoOfGenes(object),
             "\n", sep="")
        if (! short){
            cat( prefix, "genesList: \n",
                 sep="")
            show(getGenesList(object))
        } else {
            cat( prefix, "genesList: user 'getGenesList(object)'\n", sep="")
        }
        cat("\n")
      })  


# selectedGenes2ndLayerCV
setMethod("show", "selectedGenes2ndLayerCV", function(object) showWithPrefix(object))
setMethod("showWithPrefix", "selectedGenes2ndLayerCV", function(object, prefix="", short=FALSE){
        for (i in 1:length(object)){
          cat( prefix, "  [[",i, "]]  Fold ", i, "\n", sep="")
          showWithPrefix(object[[i]], paste(prefix, "  ", sep=""), short)
        }
      })



# result2LayerCV
setMethod("show", "result2LayerCV", function(object) showWithPrefix(object))
setMethod("showWithPrefix", "result2LayerCV", function(object, prefix="", short=FALSE) {
            cat( prefix, "results1stLayer:",
                length(getResults1stLayer(object)), " inner one-layer CV", sep="" )
            if (! short){
                cat("\n")
                for (i in 1:length(getResults1stLayer(object))){
                    showWithPrefix(getResults1stLayer(object)[[i]], paste(prefix, "  "), TRUE)
                }
            } else {
                cat(": use 'getResults1stLayer(object)'\n")
            }
            cat(prefix, "errorRates:\n", sep="")
                showWithPrefix(getErrorRates(object), paste(prefix, "  "), short)
            if (! short){
                cat(prefix, "selectedGenes:\n", sep="")
                    showWithPrefix(getSelectedGenes(object), paste(prefix, "  "), TRUE)
            } else {
                cat(prefix, "selectedGenes: use 'getSelectedGenes(object)'\n", sep="")
            }
            cat( prefix, "avgBestOptionValue: ",
                getAvgBestOptionValue(object),
                "\n",
                sep="")
            cat( prefix, "executionTime: ",
                getExecutionTime(object),
                "s\n",
                sep="")
              })

# resultRepeated2LayerCV
setMethod("show", "resultRepeated2LayerCV", function(object) showWithPrefix(object) )
setMethod("showWithPrefix", "resultRepeated2LayerCV",  function(object, prefix, short=FALSE) {
             cat( prefix,
                   "original2LayerCV: ",
                   length(getOriginal2LayerCV(object)),
                   "  combined 2 layer ",
                   length(getResultsLayer(getOriginal2LayerCV(object)[[1]], topic='errorRate', errorType='fold')),
                   "-folds CV",
                   sep="")
              if (! short){
                  cat("\n")
                  for (i in 1:length(getOriginal2LayerCV(object))){
                    cat( prefix, "  ",
                         "[[", i, "]] CV number ", i, "\n", sep="")
                    showWithPrefix(getOriginal2LayerCV(object)[[i]], paste(prefix,"  ", sep=""), TRUE)
                  }
              } else {
                cat(": use 'getOriginal2LayerCV(object)'\n")
              }
              cat(prefix, "summaryErrorRate:\n", sep="")
              showWithPrefix(getSummaryErrorRate(object), paste(prefix,"  ", sep=""))
              cat(prefix, "=> avgBestOptionValue: ",
                getAvgBestOptionValue(object),
                " (Est. Error rate: ",
                min(getResultsLayer(object, topic='errorRate', errorType='cv')),")\n",
                sep="")
              cat( prefix, "executionTime: ",
                getExecutionTime(object),
                "s\n",
                sep="")
              })

# finalClassifier
setMethod("show", "finalClassifier", function(object) showWithPrefix(object) )
setMethod("showWithPrefix", "finalClassifier",  function(object, prefix, short=FALSE) {
    cat(prefix, "genesFromBestToWorst:\n",
        prefix,
        "  ",
        paste(getGenesFromBestToWorst(object), collapse=" "),
        "\n",
        sep="")
    if (!short){
        cat(prefix, "models:\n", sep="")
        show(getModels(object))
    } else {
        cat(prefix, "models: use 'getModels(object)'\n", sep="")
    }

})

# assessment
setMethod("show", "assessment", function(object) showWithPrefix(object))
setMethod("showWithPrefix", "assessment", function(object, prefix="", short=FALSE) {
              cat( prefix, "assessment\n", sep="")
              cat(  prefix, "noFolds1stLayer:",
                    getNoFolds1stLayer(object), "\n", sep="")
              cat(  prefix, "noFolds2ndLayer:",
                    getNoFolds2ndLayer(object), "\n", sep="")
              cat(  prefix, "classifierName:", getClassifierName(object),
                    sep="")
              if ( getClassifierName(object) == "svm"){
                cat( " (svmKernel: ",
                     getSvmKernel(object),
                     ")\n",
                     sep="")
              }
              cat(  prefix, "featureSelectionMethod:", 
                    getFeatureSelectionMethod(object),
                    "\n",
                    sep="")
              cat(prefix, "featureSelectionOptions\n", sep="")
              showWithPrefix(getFeatureSelectionOptions(object), paste(prefix, "  "))
              cat(  prefix, "typeFoldCreation:",
                    getTypeFoldCreation(object), "\n", sep="")
              cat(  prefix, "noOfRepeats:", getNoOfRepeats(object), "\n", sep="")
              cat(prefix, "dataset\n", sep="")
              #show(getDataset(object))
              cat(": use 'getDataset(object)'\n")
              if ( ! is.null(getResult1LayerCV(object)) ){
                cat("\n", prefix, "resultRepeated1LayerCV\n", sep="")
                showWithPrefix( getResult1LayerCV(object), paste(prefix,"  "), TRUE)
              } else {
                cat("\n", prefix, "No Results for external CV (1 layer)\n", sep="")
              }
              if ( ! is.null( getResult2LayerCV(object)) ){
                cat("\n", prefix, "resultRepeated2LayerCV\n", sep="")
                showWithPrefix( getResult2LayerCV(object), paste(prefix,"  "), TRUE)
              } else {
                cat("\n", prefix, "No Results for 2 layers external CV\n", sep="")
              }
              if ( ! is.null( getFinalClassifier(object)) ){
                cat("\n", prefix, "finalClassifier\n", sep="")
                showWithPrefix( getFinalClassifier(object), paste(prefix,"  "), TRUE)
              } else {
                cat("\n", prefix, "Final Classifier has not been computed yet\n", sep="")
              }
            } )