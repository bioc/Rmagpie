# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#--------------------------- one_layer_ext_CV --------------------------------
# Function to perform a one-layer external cross-validation (called by
# runOneLayerCV to be applied directly on an assessment)
#
# Author: Camille Maumet
# Creation: Feb. 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Run a repeated one-layer CV and return an object of class
# resultRepeatedOneLayerCV containing all the information of this CV.
#
# @param    eset (ExpressionSet) microarray data
#           testFoldInds (list(numeric)) List of the fold indices, if NULL
#                   the folds are automatically generated
#           noFolds (numeric) number of folds for the CV
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
# @return   object (resultRepeatedOneLayerCV) storing the results of the
#               repeated one-layer CV
#-------------------------------------------------------------------------------
oneLayerExtCV <- function(eset, testFoldInds=NULL, noFolds=10, foldCreation,
                classifierName, featureSelectionMethod,
                verbose=getOption('verbose'), kernel="linear", noOfRepeats=1,
                optionValues) {
    if (verbose){
        print("********** One-Layer Cross-Validation **********")
    }

    # ------ Variables ------

    # Starting time
    beginning <- Sys.time()
    # Number of classes
    noClasses <- length(unique(pData(eset)$type))
    # Number of samples
    noSamples <- length(eset$type)
    # Total number of features (genes) in the dataset
    noFeatures <- length(featureNames(eset))
    # Sizes of subsets to be tried by the forward feature selection algorithm
    noModels <- length(optionValues)

    # If the number of folds is set to 1, we perform a leave one out CV
    if (noFolds == 1){
        noFolds <- noSamples
    }

    # -- Variables to be updated during the computation
    # Cross validated error rate for each size of subset
    cvErrorRate <- vector()
    # Standard error for each size of subset
    seErrorRate <- vector()
    # Error rate for each fold and each size of subsets in a given repeat
    errorRatePerFold <- matrix()
    # Labels of the output classes
    classes <- factor(levels(eset$type))
    # List of single one layer CV corresponding to a given repeated one-layer CV
    resultOriginal1LayerCV <- alist()

    # Cross validated error rate for each size of subset and each repeat
    allCvErrorRates <- matrix(nrow=noOfRepeats, ncol=noModels)
    # Cross validated error rate for each size of for each fold and each model
    # added over the repeats
    allErrorRatesPerFold <- matrix(0,ncol=noModels, nrow=noFolds)
    # Sum of the CV class error rate for each model from each fold and each repeat
    sumClassesErrorRate <-  matrix(0, ncol=noModels, nrow=noClasses)
    # Genes selected for a given size of subset, a given repeat and a given fold
    selectedGenesPerModel <- alist()
    selectedGenesPerModel[[noModels]] <- character(0)

    # ------ Core of the function ------
    for (i in 1:noOfRepeats){
        # Starting time of the current repeat
        startRepeat <- Sys.time()

        if (verbose){
            print(paste("---------- Repeat", i, " ----------"))
        }

        # Results of the cv perform for each size of subset
        resCvPerModel <- alist()

        # ------ Divide in <noFolds> folds ------
        # If the folds have been passed by the outer layer we keep them for the 1st repeat
        if (is.null(testFoldInds) || i>1){
            testFoldInds <- divideInFolds(  noFolds=noFolds,
                                            typeFold=foldCreation,
                                            train=FALSE,
                                            outputClasses=pData(eset)$type)
        }
        if (classifierName == 'svm') {
                if (featureSelectionMethod == 'rfe'){
                    resCvPerModel <- performACvPerModelWithSVMAndRFE(
                                        dataset = eset,
                                        testFoldInds = testFoldInds,
                                        foldCreation = foldCreation,
                                        classes = classes,
                                        #noTopGenes = noTopGenes,
                                        verbose = verbose,
                                        kernel = kernel,
                                        optionValues = optionValues)
                } else {
                    stop("In 'oneLayerExtCV' if 'classifierName' is 'svm' then 'featureSelectionMethod' must be 'rfe'")
                }
        } else {
            if (classifierName == 'nsc') {
                if (featureSelectionMethod == 'nsc') {
                    resCvPerModel <- performACvPerModelWithNSCAndNSC(
                                        dataset = eset,
                                        testFoldInds = testFoldInds,
                                        foldCreation = foldCreation,
                                        classes = classes,
                                        #noTopGenes = noTopGenes,
                                        verbose = verbose,
                                        kernel = kernel,
                                        optionValues = optionValues)
                } else {
                    stop("In 'oneLayerExtCV' if 'classifierName' is 'nsc' then 'featureSelectionMethod' must be 'nsc' too")
                }
            }
# Code useful for lda or naiveBayes, unused in current version ################
###############################################################################
                #else {
#                switch(featureSelectionMethod,
#                    rfe = featureSelFun <- RFESVMFeatureSelection,
#                    nsc = featureSelFun <- NSCFeatureSelection,
#                    stop(paste("In one_layer_ext_cv, 'featureSelectionMethod' ", featureSelectionMethod, " is not supported")))
#                switch (classifierName,
#                    lda = trainTestFun <- MLITrainTest,
#                    naiveBayes = trainTestFun <- MLITrainTest,
#                    stop(paste("In one_layer_ext_cv, 'classifierName' ", classifierName, " is not supported" ))
#                )
#                    # -->  Perform a CV: For each fold
#                    #      and for each size of subset :
#                    #       * Perform an RFE on the training data to find the best subset of
#                    #         genes of a given size
#                    #       * Train with the best subset of genes of a given size and test
#                    #         on the test data
#                    resCvPerModel <- performACvPerModel(dataset = eset,
#                                                        testFoldInds = testFoldInds,
#                                                        foldCreation = foldCreation,
#                                                        classifierName = classifierName,
#                                                        classes = classes,
#                                                        #noTopGenes = noTopGenes,
#                                                        verbose = verbose,
#                                                        kernel = kernel,
#                                                        optionValues = optionValues,
#                                                        featureSelectionFun = featureSelFun,
#                                                        trainTestFun = trainTestFun)
#            }
# Code useful for lda or naiveBayes, unused in current version ################
###############################################################################
        }
        # Storage of the error rates
        # CV error rate
        cvErrorRate <- resCvPerModel$cvErrorRates
        # Keep the value of the CV error rate for each model and each repeat
        # Will be useful to determine the standard error and avg CV error rate
        noPracticalThres <- length(cvErrorRate)
        if (i==1 &&  noPracticalThres!=noModels){
            # Cross validated error rate for each size of subset and each repeat
            allCvErrorRates <- matrix(nrow=noOfRepeats, ncol=noPracticalThres)
            # Cross validated error rate for each size of for each fold of each repeat
            allErrorRatesPerFold <- matrix(0, ncol=noPracticalThres, nrow=noFolds)
            # Sum of the CV class error rate for each model from each fold and each repeat
            sumClassesErrorRate <-  matrix(0, ncol=noPracticalThres, nrow=noClasses)
        }

        allCvErrorRates[i,] <- cvErrorRate
        # class error rates
        classErrorRates <- resCvPerModel$classErrorRates

        sumClassesErrorRate <- sumClassesErrorRate + classErrorRates # sum of class errors
        # error rate per fold
        # Keep the value of the error rate per fold for each fold, each size of
        # subset and each repeat
        errorRatePerFold <- resCvPerModel$errorRatesPerFold
        allErrorRatesPerFold <- allErrorRatesPerFold + errorRatePerFold

        # Store the error rates corresponding to this repeats of one-layer CV
        errorRate1stLayerCV <- new( "errorRate1stLayerCV",
                    errorRatePerFold = errorRatePerFold,
                    noSamplesPerFold = resCvPerModel$noSamplesPerFold, # Number of samples per fold
                    cvErrorRate = new("cvErrorRate",
                                        cvErrorRate = cvErrorRate,
                                        seErrorRate = resCvPerModel$seErrorRates, # standard error on CV error
                                        classErrorRates = classErrorRates))


        # Best number of genes (size of subset corresponding to the smallest)
        # CV error rate
        bestOptionValue <- resCvPerModel$bestOptionValue
        if (verbose){
            print(paste("=> value of the best option:", paste(bestOptionValue, collapse=",")))
        }

        # Determine the frequency of selection of the genes among the folds
        # List storing the frequencies and corresponding list of genes
        listDfTopGenes <- alist()
        if (noPracticalThres > 1) {
            for (j in 1:(noPracticalThres-1)){
                listDfTopGenes[[j]] <- matrix()
                listDfTopGenes[[j]] <- new("frequencyTopGenePerOneModel",
                        frequencyTopGenes(resCvPerModel$selectedGenes[[j]]))
                # Store the genes selected for each size of subset and each fold
                selectedGenesPerModel[[j]] <- c( selectedGenesPerModel[[j]],
                                                 resCvPerModel$selectedGenes[[j]])
            }
            if (featureSelectionMethod == 'rfe') {
                # In the last subset where all the genes are selected, they all have
                # a frequency of 1
                listDfTopGenes[[noPracticalThres]] <- new( "frequencyTopGenePerOneModel",
                   list( new(  "frequencyGenes",
                         frequ=1,
                         genesList=featureNames(eset))))

                selectedGenesPerModel[[noPracticalThres]] <- c( selectedGenesPerModel[[noPracticalThres]],
                                                        resCvPerModel$selectedGenes[[noPracticalThres]])
            } else {
                listDfTopGenes[[noPracticalThres]] <- matrix()
                listDfTopGenes[[noPracticalThres]] <- new("frequencyTopGenePerOneModel",
                        frequencyTopGenes(resCvPerModel$selectedGenes[[noPracticalThres]]))
                # Store the genes selected for each size of subset and each fold
                selectedGenesPerModel[[noPracticalThres]] <- c( selectedGenesPerModel[[noPracticalThres]],
                                                 resCvPerModel$selectedGenes[[noPracticalThres]])
            }
        }

        # Store the selected genes corresponding to this repeats of one-layer CV
        selectedGenes1stLayerCV <- new("selectedGenes1stLayerCV",
              selectedGenesPerFold = resCvPerModel$selectedGenes,
              frequencyTopGene = listDfTopGenes)

        # Store the information of the current repeat
        resultOriginal1LayerCV[[i]] <- new("resultSingle1LayerCV",
              errorRates = errorRate1stLayerCV,
              selectedGenes = selectedGenes1stLayerCV,
              executionTime = as.numeric( difftime( Sys.time(),
                                                    startRepeat,
                                                    units="secs")),
              bestOptionValue = bestOptionValue)
    }

    # Create a summary of the error rates from the repeats
    summaryErrorRate1stLayerCV <- new("cvErrorRate",
        # This works but is too long, it's shorter to avg the cvErrorRate
        # cvErrorRate = apply(allErrorRatesPerFold, 2, weighted.mean, pds),
        cvErrorRate = apply(allCvErrorRates, 2, sum)/noOfRepeats,
        # Wrong formula
        #seErrorRate = sqrt(apply(allErrorRatesPerFold, 2, var)/(noFolds*noOfRepeats)),
        seErrorRate = sqrt(apply(allErrorRatesPerFold/noOfRepeats, 2, var)/(noFolds)),
        classErrorRates = sumClassesErrorRate/noOfRepeats )


    # Create a summary of the frequency of selection of the genes
    summaryFrequencyTopGenes <-alist()
    for (j in 1:noPracticalThres){
        summaryFrequencyTopGenes[[j]] <- array()
        summaryFrequencyTopGenes[[j]] <- new("frequencyTopGenePerOneModel",
                                  frequencyTopGenes(selectedGenesPerModel[[j]]))
    }
    overallBestNoOfGenes <- optionValues[which.min(apply(allCvErrorRates, 2, mean, na.rm=TRUE))]

    if (verbose){
        print(paste("=> OVERALL value of the best option:", overallBestNoOfGenes))
    }

    # Finally, store all the results of this repeated one-layer CV
    result <- new("resultRepeated1LayerCV",
                        original1LayerCV = resultOriginal1LayerCV,
                        summaryErrorRate = summaryErrorRate1stLayerCV,
                        bestOptionValue = overallBestNoOfGenes,
                        summaryFrequencyTopGenes = summaryFrequencyTopGenes,
                        executionTime=as.numeric(difftime(Sys.time(), beginning, units="secs"))
                      )
    if (verbose){
        print("********** END of One-Layer Cross-Validation **********")
    }
    return(result)
}