# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#-------------------------- assessment-accessors -------------------------------
# Class assessment Accessors and replacement methods
#
# Object storing the options and results of an assessment.
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#---------------------------- Getters ------------------------------------------
# Method to access the results from both 1 Layer CV and 2 Layers of CV
setMethod("getResults", "assessment", function(object, layer, ...){
        deepLayer <- length(layer)
        layerType <- layer[1]
        if ( deepLayer > 4 || deepLayer < 1 ) {
            stop("'layer' must contains at least 1 number and at most 4")
        }
        if ( layerType != 1 && layerType != 2){
            stop("First number in 'layer' must be '1' for the one-layer CV or '2' for 2-layers CV")
        }
        if ( layerType== 1 && is.null(object@resultRepeated1LayerCV)){
            stop("Result of one-layer CV are not available yet, call 'runOneLayerExtCV'")
        }
        if ( layerType==2 && is.null(object@resultRepeated2LayerCV)){
            stop("Result of two-layers CV are not available yet, call 'runTwoLayerExtCV'")
        }
        if ( deepLayer>2 && layerType!=2){
            stop("'layer' can contains more than 2 numbers only for 2-layers CV ('layer' starting with a 2)")
        }
        if (layerType == 1) { # 1-layer CV
            if (deepLayer == 1){ # Summary of the one-layer
                getResultsLayer(object@resultRepeated1LayerCV, ...)
            } else {
                if (layer[2] > length( object@resultRepeated1LayerCV@original1LayerCV)) {
                    stop(paste("In one-layer CV, Repeat ", layer[2], " does not exist", sep=""))
                } else {
                    # Single one-layer CV
                    getResultsLayer(object@resultRepeated1LayerCV@original1LayerCV[[layer[2]]], ...)
                }
            }
        } else { # 2-layers CV
            if (deepLayer == 1){ # Summary of the two-layers
                getResultsLayer(object@resultRepeated2LayerCV, ...)
            } else {
                if (deepLayer == 2) {
                    # Single two-layers CV
                    if (layer[2] > length(object@resultRepeated2LayerCV@original2LayerCV)) {
                        stop(paste("In two-layers CV, Repeat ", layer[2], " does not exist", sep=""))
                    } else {
                        getResultsLayer( object@resultRepeated2LayerCV@original2LayerCV[[layer[2]]], ...)
                    }
                } else {
                    if (deepLayer == 3){
                        if (layer[2] > length(object@resultRepeated2LayerCV@original2LayerCV)) {
                            stop(paste("In two-layers CV, Repeat ", layer[2], " does not exist", sep=""))
                        } else {
                            if(layer[3] > length(object@resultRepeated2LayerCV@original2LayerCV[[layer[2]]]@results1stLayer)) {
                                stop(paste("In two-layers CV Repeat ", layer[2],", Fold ", layer[3], " does not exist", sep=""))
                            } else {
                                getResultsLayer(object@resultRepeated2LayerCV@original2LayerCV[[layer[2]]]@results1stLayer[[layer[3]]], ...)
                            }
                        }
                    } else {
                        if (layer[2] > length(object@resultRepeated2LayerCV@original2LayerCV)) {
                            stop(paste("In two-layers CV, Repeat ", layer[2], "does not exist"))
                        } else {
                            if(layer[3] > length(object@resultRepeated2LayerCV@original2LayerCV[[layer[2]]]@results1stLayer)) {
                                stop(paste("In two-layers CV Repeat ", layer[2],", Fold ", layer[3], " does not exist", sep=""))
                            } else {
                                if(layer[4] > length(object@resultRepeated2LayerCV@original2LayerCV[[layer[2]]]@results1stLayer[[layer[3]]]@original1LayerCV)) {
                                    stop(paste("In two-layers CV Repeat ", layer[2],", Fold ", layer[4], ", Repeat ", layer[3], " does not exist", sep=""))
                                } else {
                                    getResultsLayer(object@resultRepeated2LayerCV@original2LayerCV[[layer[2]]]@results1stLayer[[layer[3]]]@original1LayerCV[[layer[4]]], ...)
                                }
                            }
                        }
                    }
                }
            }
        }
    })

setMethod("getResultsLayer", "resultRepeated1LayerCV", function(object, topic, genesType, errorType){
    switch( topic,
            errorRate =
                if (missing(errorType)) {
                    object@summaryErrorRate
                } else {
                    switch( errorType,
                            all = object@summaryErrorRate,
                            cv = object@summaryErrorRate@cvErrorRate,
                            se = object@summaryErrorRate@seErrorRate,
                            class = object@summaryErrorRate@classErrorRates,
                            stop("In 'getResults' with 'layer'=1, 'errorType' must be 'all', 'cv', 'se', 'class' or missing"))
                } ,
            genesSelected =
                if (missing(genesType) || genesType=='frequ'){
                        object@summaryFrequencyTopGenes
                } else {
                    stop("In 'getResults' with 'layer'=1, 'genesType' must be 'frequ' or missing")
                },
            bestOptionValue = object@bestOptionValue,
            executionTime = object@executionTime,
            stop("In 'getResults' with 'layer'=1, 'topic' must be 'errorRate', 'genesSelected', 'bestOptionValue' or 'executionTime'")
            )
})

setMethod("getResultsLayer", "resultSingle1LayerCV", function(object, topic, genesType, errorType){
    if (missing(genesType)) {
        # By default display the frequency table of the genes
        genesType <- 'frequ'
    }
    switch( topic,
            errorRate =
                if (missing(errorType)) {
                    object@errorRates
                } else {
                    switch( errorType,
                            all = object@errorRates,
                            fold = object@errorRates@errorRatePerFold,
                            noSamplesPerFold = object@errorRates@noSamplesPerFold,
                            cv = object@errorRates@cvErrorRate@cvErrorRate,
                            se = object@errorRates@cvErrorRate@seErrorRate,
                            class = object@errorRates@cvErrorRate@classErrorRates,
                            stop("In 'getResults' with 'layer'=c(1,i), 'errorType' must be 'all', 'fold', 'noSamplesPerFold', 'cv', 'se', 'class' or missing") )
                } ,
            genesSelected =
                switch( genesType,
                        frequ = object@selectedGenes@frequencyTopGene,
                        fold = object@selectedGenes@selectedGenesPerFold,
                        stop("In 'getResults' with 'layer'=c(1,i), 'genesType' must be 'frequ', 'fold' or missing")
                ),
            bestOptionValue = object@bestOptionValue,
            executionTime = object@executionTime,
            stop("In 'getResults' with 'layer'=c(1,i), 'topic' must be 'errorRate', 'genesSelected', 'bestOptionValue' or 'executionTime'")
            )
})

setMethod("getResultsLayer", "resultRepeated2LayerCV", function(object, topic, genesType, errorType){
   switch(  topic,
            errorRate =
                if (missing(errorType)) {
                    object@summaryErrorRate
                } else {
                    switch( errorType,
                            all = object@summaryErrorRate,
                            cv = object@summaryErrorRate@finalErrorRate,
                            se = object@summaryErrorRate@seFinalErrorRate,
                            class = object@summaryErrorRate@classErrorRates,
                            stop("In 'getResultsLayer', 'errorType' must be 'all', 'cv', 'se', 'class' or missing")
                    )
                } ,
            bestOptionValue = object@avgBestOptionValue,
            executionTime = object@executionTime,
            stop("In 'getResultsLayer', 'topic' must be 'errorRate', 'bestOptionValue' or 'executionTime'")
    )
})

setMethod("getResultsLayer", "result2LayerCV", function(object, topic, genesType, errorType){
   switch(  topic,
            errorRate =
                if (missing(errorType)) {
                    object@errorRates
                } else {
                    switch( errorType,
                            all = object@errorRates,
                            fold = object@errorRates@errorRatePerFold,
                            noSamplesPerFold = object@errorRates@noSamplesPerFold,
                            cv = object@errorRates@cvErrorRate@finalErrorRate,
                            se = object@errorRates@cvErrorRate@seFinalErrorRate,
                            class = object@errorRates@cvErrorRate@classErrorRates,
                            stop("In 'getResultsLayer', 'errorType' must be 'all', 'fold', 'cv', 'se', 'class' or missing")
                    )
                } ,
            genesSelected = if (missing(genesType) || genesType=='fold'){
                                object@selectedGenes
                            } else {
                                stop("In 'getResults' with 'layer'=c(1,i), 'genesType' must be 'fold' or missing")
                            },
            bestOptionValue = object@avgBestOptionValue,
            executionTime = object@executionTime,
            stop("In 'getResultsLayer', 'topic' must be 'errorRate', 'genesSelected', 'bestOptionValue' or 'executionTime'")
    )
})


#---------------------------- Setters ------------------------------------------
# Return the informations corresponding to the dataset of the assessment
setMethod("getDataset", "assessment",
    function(object) {
        return(object@dataset)
    }
)

setMethod("getFeatureSelectionOptions", "assessment",
    function(object, topic) {
        if (missing(topic)){
            return(object@featureSelectionOptions)
        }
        if (topic == 'optionValues') return(getOptionValues(object@featureSelectionOptions))
        if (topic == 'noOfOptions') return(getNoOfOptions(object@featureSelectionOptions))
        if (is(object@featureSelectionOptions, "thresholds")){
            switch(topic,
            thresholds = return(getOptionValues(object@featureSelectionOptions)),
            noThresholds = return(getNoThresholds(object@featureSelectionOptions)),
            stop(paste("In 'getFeatureSelectionOptions', with object of class 'thresholds', ",
                    "'topic' must be 'optionValues', 'noOfOptions', 'thresholds', 'noThresholds' or missing ",
                     "(and not", topic, ")", sep="" ))
            )
        } else {
            switch(topic,
            maxSubsetSize = return(getMaxSubsetSize(object@featureSelectionOptions)),
            subsetsSizes = return(getSubsetsSizes(object@featureSelectionOptions)),
            speed = return(getSpeed(object@featureSelectionOptions)),
            noModels = return(getNoModels(object@featureSelectionOptions)),
            stop(paste("In 'getFeatureSelectionOptions', with object of class 'geneSubsets', ",
                    "'topic' must be 'optionValues', 'noOfOptions', 'maxSubsetSize', 'subsetsSizes', ",
                    "'speed' or 'noModels' or missing (and not ", topic, ")", sep="" ))
            )
        }
    }
)


setMethod("getNoFolds1stLayer", "assessment", function(object) object@noFolds1stLayer)
setMethod("getNoFolds2ndLayer", "assessment", function(object) object@noFolds2ndLayer)
setMethod("getClassifierName", "assessment", function(object) object@classifierName)
setMethod("getFeatureSelectionMethod", "assessment", function(object) object@featureSelectionMethod)
setMethod("getSvmKernel", "assessment", function(object) object@svmKernel)
setMethod("getTypeFoldCreation", "assessment", function(object) object@typeFoldCreation)
setMethod("getResult1LayerCV", "assessment", function(object) object@resultRepeated1LayerCV)
setMethod("getResult2LayerCV", "assessment", function(object) object@resultRepeated2LayerCV)
setMethod("getNoOfRepeats", "assessment", function(object) object@noOfRepeats)

setMethod("getFinalClassifier", "assessment",
    function(object, topic) {
        if (missing(topic)){
            return(object@finalClassifier)
        }
        switch(topic,
            genesFromBestToWorst = return(getGenesFromBestToWorst(object@finalClassifier)),
            models = return(getModels(object@finalClassifier)),
            stop("In 'getFinalClassifier', 'topic' must be 'genesFromBestToWorst', 'models' or missing")
        )
    }
)

setReplaceMethod("getFeatureSelectionOptions", "assessment",
    function(object, topic, value) {
        if (! (is.null(getResult1LayerCV(object)) && is.null(getResult2LayerCV(object))
                && is.null(getFinalClassifier(object))) ){
            stop(paste("In 'getFeatureSelectionOptions': featureSelectionOptions",
                        "can't be modified if one-layer CV, two layers CV or the final classifier",
                        "has already been computed"))
        }
        if (missing(topic)){
            object@featureSelectionOptions <- value
        } else {
            if (topic == 'optionValues') {
                getOptionValues(object@featureSelectionOptions) <- value
            } else {
                if (topic == 'noOfOptions') {
                    getNoOfOptions(object@featureSelectionOptions) <- value
                } else {
                    if (is(object@featureSelectionOptions, "thresholds")){
                        switch(topic,
                        thresholds = getOptionValues(object@featureSelectionOptions) <- value,
                        noThresholds = getNoThresholds(object@featureSelectionOptions) <- value,
                        stop(paste("In 'getFeatureSelectionOptions', with object of class 'thresholds', ",
                                "'topic' must be 'optionValues', 'noOfOptions', 'thresholds', 'noThresholds' or missing ",
                                 "(and not", topic, ")", sep="" ))
                        )
                    } else {
                        switch(topic,
                        maxSubsetSize = getMaxSubsetSize(object@featureSelectionOptions) <- value,
                        subsetsSizes = getSubsetsSizes(object@featureSelectionOptions) <- value,
                        speed = getSpeed(object@featureSelectionOptions) <- value,
                        noModels = getNoModels(object@featureSelectionOptions) <- value,
                        stop(paste("In 'getFeatureSelectionOptions', with object of class 'geneSubsets', ",
                                "'topic' must be 'optionValues', 'noOfOptions', 'maxSubsetSize', 'subsetsSizes', ",
                                "'speed' or 'noModels' or missing (and not ", topic, ")", sep="" ))
                        )
                    }
                }
            }
        }
        object
    }
)

setReplaceMethod("getDataset", "assessment",
    function(object, value) {
        if (! (is.null(getResult1LayerCV(object)) && is.null(getResult2LayerCV(object))
                && is.null(getFinalClassifier(object))) ){
            stop(paste("In 'getDataset': dataset",
                        "can't be modified if one-layer CV, two layers CV or the final classifier",
                        "has already been computed"))
        }
        object@dataset <- value
        object
    }
)

setReplaceMethod("getNoFolds1stLayer", "assessment",
    function(object, value) {
        if (! (is.null(getResult1LayerCV(object)) && is.null(getResult2LayerCV(object))
                && is.null(getFinalClassifier(object))) ){
            stop(paste("In 'getNoFolds1stLayer': noFolds1stLayer",
                        "can't be modified if one-layer CV, two layers CV or the final classifier",
                        "has already been computed"))
        }
        object@noFolds1stLayer <- value
        object
    }
)
setReplaceMethod("getNoFolds2ndLayer", "assessment",
    function(object, value) {
        if (! (is.null(getResult1LayerCV(object)) && is.null(getResult2LayerCV(object))
                && is.null(getFinalClassifier(object))) ){
            stop(paste("In 'getNoFolds2ndLayer': noFolds2ndLayer",
                            "can't be modified if one-layer CV, two layers CV or the final classifier",
                            "has already been computed"))
        }
        object@noFolds2ndLayer <- value
        object
    }
)


setReplaceMethod("getClassifierName", "assessment",
    function(object, value) {
        stop(paste("In 'getClassifierName': classifierName",
                        "can't be modified if one-layer CV, two layers CV or the final classifier",
                        "has already been computed"))
        object@classifierName <- value
        object
    }
)

setReplaceMethod("getSvmKernel", "assessment",
    function(object, value) {
        if (! (is.null(getResult1LayerCV(object)) && is.null(getResult2LayerCV(object))
                && is.null(getFinalClassifier(object))) ){
            stop(paste("In 'getSvmKernel': svmKernel",
                        "can't be modified if one-layer CV, two layers CV or the final classifier",
                        "has already been computed"))
        }
        object@svmKernel <- value
        object
    }
)

setReplaceMethod("getTypeFoldCreation", "assessment",
    function(object, value) {
        if (! (is.null(getResult1LayerCV(object)) && is.null(getResult2LayerCV(object))
                && is.null(getFinalClassifier(object))) ){
            stop(paste("In 'getTypeFoldCreation': typeFoldCreation",
                        "can't be modified if one-layer CV, two layers CV or the final classifier",
                        "has already been computed"))
        }
        object@typeFoldCreation <- value
        object
    }
)

setReplaceMethod("getNoOfRepeats", "assessment",
    function(object, value) {
        if (! (is.null(getResult1LayerCV(object)) && is.null(getResult2LayerCV(object))
                && is.null(getFinalClassifier(object))) ){
            stop(paste("In 'getNoOfRepeats': noOfRepeats",
                        "can't be modified if one-layer CV, two layers CV or the final classifier",
                        "has already been computed"))
        }
        object@noOfRepeats <- value
        object
    }
)

setReplaceMethod("getFinalClassifier", "assessment",
    function(object, value) {
        object@finalClassifier <- value
        object
    }
)


setReplaceMethod("getResult1LayerCV", "assessment",  function(object, value) {
                    object@resultRepeated1LayerCV <- value
                    object
                 })
setReplaceMethod("getResult2LayerCV", "assessment",  function(object, value) {
                    object@resultRepeated2LayerCV <- value
                    object
                 })
