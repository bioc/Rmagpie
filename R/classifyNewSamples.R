# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#-------------------------- classifyNewSamples ---------------------------------
# This file contains all the function and methods useful to classify new samples
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Classify one or several new sample(s) stored in the file <newSampleFile>
# by applying the final classifier (genes found by performing an RFE on the
# whole data set or by thresholding the NSC on the whole data set) for a subset
# of <optionValue> genes or a threshold of <optionValue>.
#
# @param object (assessment) Object of class assessment
#        newSamplesFile (character) URL of the file containing the new samples
#           gene expression profiles(in columns)
#        optionValue (numric) Size of subset (or thresholds) to be considered to
#           build the final classifier. If missing, pick the one corresponding
#           to the best error rate by one-layer of cross-validation. If this has
#           not been performed yet returns an error.
#-------------------------------------------------------------------------------
setMethod("classifyNewSamples", "assessment",
    function(object, newSamplesFile, optionValue){
        # If the size of the subset has not been set we pick the best subset size
        # according to the one-layer of cross-validation
        if (missing(optionValue)){
            if (is.null(getResult1LayerCV(object))){
                stop("One external layer of CV must be performed before classifying a new sample without 'optionValue' specified: call 'runOneLayerExtCV'
                        or call 'classifyNewSample' with a non-missing 'optionValue' argument.")
            } else {
                # Select the size corresponding to the best subset
                optionValue <- getResults(object, 1, 'bestOptionValue')
                # Find the position of the subset in geneSubsets
                model <- which(getFeatureSelectionOptions(object, 'optionValues') == optionValue)[1]
            }
        } else {
            # Select the genes corresponding to the size 'subsetSize' if there is
            # no subset corresponding to this size then the closest bigger
            # subset is chosen
            if (optionValue > max(getFeatureSelectionOptions(object, 'optionValues'))){
                # If there if bigger subset the smallest subset is chosen arbitrarily
                model = 1
            } else {
                model <- which(getFeatureSelectionOptions(object, 'optionValues') >= optionValue)[1]
            }
        }

        # If the final classifier has not been find yet
        if (is.null(getFinalClassifier(object))){
            stop('You must run "findFinalClassifier" before classifying new samples')
        }

        # Read the new samples file
        newSamples <- read.table(newSamplesFile, header=TRUE, row.names=1)
        rownames(newSamples) <- rownames(exprs(getDataset(object)))

        if (getClassifierName(object) == 'svm'){
            newSamples <- t(newSamples)
            finalClassifier <- getFinalClassifier(object, "models")[[model]]$model
            relevantGenes <- getFinalClassifier(object, "models")[[model]]$modelFeatures
            predictedClass <- predict(finalClassifier,newSamples[,relevantGenes])
        } else {
            if (getClassifierName(object) == 'nsc')  {
                finalClassifier <- getFinalClassifier(object, "models")[[1]]
                thres <- finalClassifier$threshold[model]
                predictedClass <- pamr.predict( fit=finalClassifier,
                                                newx=newSamples,
                                                threshold=thres)
                names(predictedClass) <- colnames(newSamples)
            }
    # Code useful for lda or naiveBayes, unused in current version ################
    ###############################################################################
        #else {
    #
    #            if (getFeatureSelectionMethod(object) == 'rfe'){
    #                relevantGenes <- getFinalClassifier(object, "models")[[model]]$modelFeatures
    #            } else {
    #                if (getFeatureSelectionMethod(object) == 'nsc'){
    #                    finalClassifier <- getFinalClassifier(object, "models")[[1]]
    #                    thres <- finalClassifier$threshold[length(finalClassifier$threshold) - model + 1]
    #                    mydata <- list( x=exprs(getDataset(object, 'eset')),
    #                                    y=pData(getDataset(object, 'eset'))[[1]],
    #                                    genenames=featureNames(getDataset(object, 'eset')),
    #                                    geneid=featureNames(getDataset(object, 'eset')) )
    #                    relevantGenes <- pamr.listgenes(fit=finalClassifier, data=mydata, threshold=thres)[,1]
    #                } else{
    #                    stop(paste("In 'classifyNewSamples' 'getFeatureSelectionMethod(object)' must be 'nsc' or 'svm' and not ", getFeatureSelectionMethod(object)))
    #                }
    #            }
    #
    #            trainData <- getDataset(object, 'eset')
    #            noSamples <- dim(trainData)[2]
    #
    #            #dataExprSet <- getDataset(object, 'eset')
    #            dftrainData <- (as.data.frame(t(exprs(getDataset(object, 'eset'))[relevantGenes,])))
    #            dftestData <- (as.data.frame(t(newSamples[relevantGenes,])))
    #            row.names(dftestData) <- paste("pred", row.names(dftestData), sep="")
    #
    #            allData <- as.data.frame(rbind(dftrainData, dftestData))
    #
    #            allData$type <- factor(c(as.character(pData(trainData)[[1]]),
    #                as.character(rep(pData(trainData)[[1]][1], dim(newSamples)[2]))))
    #            if (getClassifierName(object) == 'naiveBayes'){
    #                trainTest <- MLearn(type~., allData, naiveBayesI, trainInd=seq(1, noSamples))
    #            } else {
    #                if (getClassifierName(object) == 'lda'){
    #                    print("lda classification")
    #                    trainTest <- MLearn(type~., allData, ldaI, trainInd=seq(1, noSamples))
    #                }
    #            }
    #            predictedClass <- testPredictions(trainTest)
    #            names(predictedClass) <- colnames(newSamples)
    #        }
    # Code useful for lda or naiveBayes, unused in current version ################
    ###############################################################################
        }
        return(predictedClass)
})

#-------------------------------------------------------------------------------
# Train the final classifier corresponding to the assessment <object>.
#
# @param object (assessment) Object of class assessment
#-------------------------------------------------------------------------------
setMethod("findFinalClassifier", "assessment", function(object){
    optionValues <- getFeatureSelectionOptions(object, 'optionValues')
    featureNames <- featureNames(getDataset(object))

    if (getFeatureSelectionMethod(object) == 'rfe'){
        # SVM and RFE
        rfeModels <- rfe.fit(  x=aperm(exprs(getDataset(object))),
                            y=pData(getDataset(object))[[1]],
                            minf=1,
                            noSelectedFeatures=optionValues,
                            kern=getSvmKernel(object),
                            verbose=TRUE)

        svmModels <- alist()
        noModels <- length(rfeModels$TestedModels)
        for (i in 1:noModels){
            svmModels[[noModels-i+1]] <- alist()
            svmModels[[noModels-i+1]]$model <- rfeModels$TestedModels[[i]]$model
            svmModels[[noModels-i+1]]$modelFeatures <-  featureNames[rfeModels$TestedModels[[i]]$modelFeatures]
        }
        getFinalClassifier(object) <- new(  "finalClassifier",
                                        genesFromBestToWorst=featureNames[rfeModels$Flist],
                                        models=svmModels)
    } else {
        if (getFeatureSelectionMethod(object) == 'nsc'){
            # NSC
            thresholds <- getFeatureSelectionOptions(object, 'thresholds')
            mydata <- list(x=exprs(getDataset(object)),
                                                y=pData(getDataset(object))[[1]],
                                                genenames=featureNames,
                                                geneid=featureNames )
            nscModels <- pamr.train(data=mydata,
                                    threshold=thresholds)

            getFinalClassifier(object) <- new(  "finalClassifier",
                                            genesFromBestToWorst=character(0),
                                            models=list(nscModels))
        }
    }
    return(object)
})
