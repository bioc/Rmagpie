# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#--------------------------- two_layer_ext_CV --------------------------------
# Function to perform a two-layer external cross-validation (called by
# runTwoLayerCV to be applied directly on an assessment)
#
# Author: Camille Maumet
# Creation: Feb. 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Run a repeated two-layer CV and return an object of class
# resultRepeated2LayerCV containing all the information of this CV.
#
# @param    eset (ExpressionSet) microarray data
#           noFoldsExtLayer (numeric) number of folds for the CV in the outer
#               layer
#           noFoldsIntLayer (numeric) number of folds for the CV in the inner
#               layer
#           foldCreation (character) type of fold creation in "original" (the
#               one used in Christophe Ambroise's code), "naive" or "balanced"
#           classifierName (character) Name of the classifier "svm" or "nsc"
#           featureSelectionMethod (character) Method of feature selection "rfe"
#                   or "nsc"
#           verbose (logical) display debugging outputs on screen ?
#           kernel (character) Kernel for the SVM in "linear", "radial",
#                   "polynomial"
#           noOfRepeats (numeric) Number of repeats to be excuted
#           optionValues (vector(numeric)) Values of option for feature
#                   selection
#
# @return   object (resultRepeated2LayerCV) storing the results of the
#               repeated two-layer CV
#-------------------------------------------------------------------------------
twoLayerExtCV <- function(  eset,
                            noFoldsExtLayer,
                            noFoldsIntLayer,
                            foldCreation="naive",
                            classifierName,
                            featureSelectionMethod,
                            verbose=getOption('verbose'),
                            kernel="linear",
                            optionValues,
                            noOfRepeats){
    if (verbose){
        print("********** Two-Layers Cross-Validation **********")
    }

    # ------ Variables ------
    beginning <- Sys.time()  # Starting time
    noClasses <- length(unique(pData(eset)[,1]))  # Number of classes
    noFolds2 <- noFoldsExtLayer # No of folds in the 2nd layer
    noFolds1 <- noFoldsIntLayer # No of folds in the 1st layer
    noSamples <- length(eset[[1]])  # No of samples
    noSamplesPerClass <- table(eset[[1]]) # No of samples in each class
    classes <- factor(levels(eset[[1]])) # Output classes
    noFeatures <- length(featureNames(eset)) # Total number of features in the dataset

    # Variables to be updated during the computation
    errorRate <- vector("numeric", noFolds2) # Error rate on the external fold
    result1stLayerCV <- alist() # Result of the inner layer
    result2Layers <- alist() # Result of the repeats
    selectedGenes <- alist() # Genes selected by RFE on each fold of the 2nd external layer
    repAvgNoGenes <- 0 # Average best number of genes in the current repeat
    sumFinalErrorRate <- 0 # Sum of the final error rate from each repeat
    sumClassErrorRates <- 0 # Sum of the classes error rate from each repeat
    avgNoGenes <- 0  # Average best number of genes
    sumNoErrorsPerClass <- vector("numeric", noClasses) # No of errors per class made on the 2nd layer
    # Cross validated error rate for each size of for each fold added over the repeats
    allErrorRatesPerFold <- vector("numeric", noFolds2)

    # Total execution time
    executionTime <- 0

    # ------ Test of the input parameters -----
    if (noFolds2 < 1 || noFolds2 > noSamples){
        stop("Number of external folds (noFolds2) incorrect -> Must have 1 < noFolds2 < noSamples")
    }
    if (noFolds1 < 1 || noFolds1 > (noSamples*(noFolds2-1)/noFolds2)){
        stop("Number of internal folds (noFolds1) incorrect -> Must have 1 < noFolds2 < size external fold")
    }

    # ------ Core of the function ------
    # For each repeat
    for (r in 1:noOfRepeats){
        if (verbose){
            print(paste("---------- Repeat", r, " ----------"))
        }
        startRepeat <- Sys.time()
        repAvgNoGenes <- 0
        # Average the error rate per class over the folds (weigthed mean)
        summaryErrorRatePerClass <- 0
        sumNoErrorsPerClass <- 0

        # Divide in <noFolds2> folds
        extInds <- alist()
        # if <noFolds> is 1, we perform a Leave One Out CV
        if (noFolds2 == 1){
            noFolds2 <- noSamples
        }
        # Creation of the second layer folds
        extInds <- divideInFolds(   noFolds=noFolds2,
                                    typeFold=foldCreation,
                                    train=TRUE,
                                    outputClasses=pData(eset)[[1]])
        outerTestInd <- lapply(extInds, function(y) setdiff(x=1:noSamples, y))

        # Number of samples in each fold
        noSamplesPerFold <- as.numeric(lapply(extInds, length))

        # For each fold in the 2nd layer (outer layer)
        for (i in noFolds2:1 ) {
            if (verbose){
                print(paste("---------- Fold", i, " ----------"))
            }

            # If there is only a difference of one between the number of folds
            # in the outer layer and in the inner layer, the folds of the inner
            # layer are determined by the one of the outer layer
            innerTestInd <- NULL
            if (noFolds1 == noFolds2 - 1){
                innerTestInd <- alist()
#                noSamplesInRemovedFold <- length(outerTestInd[[i]])

                # If it's the last fold just keep the 'noFolds1' others as they are
                if (i==noFolds2){
                    innerTestInd <- lapply(outerTestInd, function(index) {
                        which(! is.na(match(extInds[[i]], index) ))
                    })
                    length(innerTestInd) <- noFolds1
                } else {
                    if (i > 1){
                        for (before in 1:(i-1)){
                            innerTestInd[[before]] <- which(! is.na(match(extInds[[i]], outerTestInd[[before]])))
                        }
                    }
                    for (after in (i+1):noFolds2){
                        innerTestInd[[after - 1]] <- which(! is.na(match(extInds[[i]], outerTestInd[[after]])))
                    }
                }
            }
            # Check that the indices are the same minus one fold
            #   ** sampleNamesIntIndices <- lapply(innerTestInd, function(indices) sampleNames(eset[,extInds[[i]]][, indices]))
            #   ** print(sampleNamesIntIndices)
            #   ** sampleNamesOutIndices <- lapply(outerTestInd, function(indices) sampleNames(eset[,indices]))
            #   ** print(sampleNamesOutIndices)
            # --> Perform a one-layer CV on test data to find the best number of genes
            result1stLayerCV[[i]] <- oneLayerExtCV( eset = eset[,extInds[[i]]],
                                                    noFolds = noFolds1,
                                                    foldCreation = foldCreation,
                                                    testFoldInds = innerTestInd,
                                                    classifierName = classifierName,
                                                    featureSelectionMethod = featureSelectionMethod,
                                                    #noTopGenes = noTopGenes,
                                                    verbose = FALSE,
                                                    kernel = kernel,
                                                    optionValues = optionValues,
                                                    noOfRepeats = 1)

            # Best no of genes obtained by the one-layer CV
            bestOptionValue <- getBestOptionValue(result1stLayerCV[[i]])
            if (verbose){
                print(paste("=> Value of the best option:", bestOptionValue))
            }

            # Average of the best no of genes over the folds
            repAvgNoGenes <- repAvgNoGenes + bestOptionValue

            # The RFE must be performed up to <bestOptionValue> by trying all the
            # number of subsets greater than <bestOptionValue> in <optionValues>
            upToBestSelectedFeature <- optionValues[which(optionValues==bestOptionValue):length(optionValues)]

            # --> Perform the feature selection on the whole training data to find
            #     the <bestNoGenes> best genes
            # --> Train with the <bestNoGenes> best genes and test on test data
            switch( classifierName,
                svm = trainTestFun <- SVMTrainTest,
# Code useful for lda or naiveBayes, unused in current version ################
###############################################################################
#                lda = trainTestFun <- MLITrainTest,
#                naiveBayes = trainTestFun <- MLITrainTest,
# Code useful for lda or naiveBayes, unused in current version ################
###############################################################################
                nsc = trainTestFun <- NSCTrainTest,
                stop(paste("In two_layer_ext_cv, 'classifierName' ", classifierName, " is not supported" ))
            )

            switch( featureSelectionMethod,
                rfe = featureSelectionFun <- RFESVMFeatureSelection,
                nsc = featureSelectionFun <- NSCFeatureSelection,
                stop(paste("In two_layer_ext_cv, 'featureSelectionMethod' ", featureSelectionMethod, " is not supported" ))
            )

            resValid <- performAValidation( data=eset,
                                            #noFolds=noFolds2,
                                            trainIndices=extInds[[i]],
                                            noGenes=upToBestSelectedFeature,
                                            classes=classes,
                                            classifierName=classifierName,
                                            verbose=verbose,
                                            kernel=kernel,
                                            featureSelectionFun=featureSelectionFun,
                                            trainTestFun=trainTestFun )

            # -- Store the genes selected for each fold
            selectedGenes[[i]] <- new("selectedGenes",
                              genesList=resValid$relevantGenes,
                              optionValue=bestOptionValue,
                              noOfGenes=length(resValid$relevantGenes) )

            # -- Store the error rates
            # error rate for each fold
            errorRate[i] <- resValid$errorRate
            # Total number of errors per class
            #resValid$noErrorsPerClass <- 5
            sumNoErrorsPerClass <- sumNoErrorsPerClass + resValid$noErrorsPerClass
            #errorRatePerClass[i,] <- resValid$classErrorRates
        }
        # Average the error rate per class over the folds (weigthed mean)
        summaryErrorRatePerClass <- sumNoErrorsPerClass/noSamplesPerClass
        summaryErrorRatePerClass <- as.numeric(summaryErrorRatePerClass)
        names(summaryErrorRatePerClass) <- classes

        # -- Store the error rates
        finalErrorRate <- weighted.mean(errorRate, w=noSamplesPerFold)
        allErrorRatesPerFold <- allErrorRatesPerFold + errorRate
        cvErrorRates <- new( "cvErrorRate2ndLayer",
                                finalErrorRate=finalErrorRate,
                                seFinalErrorRate=sqrt(var(errorRate)/noFolds2),
                                classErrorRates=summaryErrorRatePerClass)
        errorRate2ndLayerCV <- new("errorRate2ndLayerCV",
                                    errorRatePerFold=errorRate,
                                    noSamplesPerFold=noSamplesPerFold,
                                    cvErrorRate=cvErrorRates)

        # -- Store the genes selected
        selectedGenes2ndLayerCV <- new("selectedGenes2ndLayerCV", selectedGenes)

        avgBestOptionValue <- repAvgNoGenes/noFolds2

        result2Layers[[r]] <- new("result2LayerCV",
                        results1stLayer=result1stLayerCV,
                        errorRates=errorRate2ndLayerCV,
                        selectedGenes=selectedGenes2ndLayerCV,
                        executionTime=as.numeric(difftime(Sys.time(), startRepeat, units="secs")),
                        avgBestOptionValue=avgBestOptionValue)

        avgNoGenes <- avgNoGenes + avgBestOptionValue
        sumFinalErrorRate <- sumFinalErrorRate + finalErrorRate
        sumClassErrorRates <- sumClassErrorRates + summaryErrorRatePerClass
    }

    result <- new("resultRepeated2LayerCV",
                        original2LayerCV=result2Layers,
                        summaryErrorRate=new("cvErrorRate2ndLayer",
                                    finalErrorRate=sumFinalErrorRate/noOfRepeats,
                                    seFinalErrorRate=sqrt(var(allErrorRatesPerFold/noOfRepeats)/(noFolds2)),
                                    classErrorRates=sumClassErrorRates/noOfRepeats),
                        avgBestOptionValue=avgNoGenes/noOfRepeats,
                        executionTime=as.numeric( difftime( Sys.time(),
                                                    beginning,
                                                    units="secs")))

    if (verbose){
        print("********** END Two-Layers Cross-Validation **********")
    }

    return(result)
}