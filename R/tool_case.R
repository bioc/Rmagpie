# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#----------------------------- tool_case ---------------------------------------
# This toolcase contain all the methods required both by the double layer of
# external cv and the one layer of external cv
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Perform a cv per model define by <optionValues> using the
# NSC on <data> with the folds <testFoldInds> and return
# the corresponding error rates for each model and the best no of genes:
# -->  Perform a CV: For each fold
#      and for each size of subset :
#       * Perform an RFE on the training data to find the best subset of
#         genes of a given size
#       * Train with the best subset of genes of a given size and test
#         on the test data
#
# @param  dataset (ExpressionSet)Data on which the classification is done
#         testFoldInds (list(numeric)) List of fold indices
#         optionValues (numeric)  Vector containing the thresholds to be
#                                 considered during feature selection
#         classes (vector(character)) Name of the classes
#         verbose (logical)       Debug variable
#         kernel (character) Ignored (just here for the similarity with other
#                                       funcitons)
# @return resCvPerModel (list) Result of cv per model:
#               resCvPerModel$bestOptionValue: Value of the best option
#               resCvPerModel$cvErrorRates: cv error rate for each threshold
#               resCvPerModel$classErrorRates: class error rate for each thresh
#               resCvPerModel$noSamplesPerFold: No of samples in each fold
#               resCvPerModel$errorRatesPerFold: Error rate in each fold
#               resCvPerModel$seErrorRates: standard error on cv error rate
#               resCvPerModel$selectedGenes: genes selected with each thershold
#-------------------------------------------------------------------------------
performACvPerModelWithNSCAndNSC  <- function( dataset,
                                testFoldInds,
                                foldCreation,
                                optionValues,
                                classes,
                                verbose=FALSE,
                                kernel="linear") {
    # --- Variables
    start <- Sys.time() # starting time
    noModels <- length(optionValues) # Number of different thresholds
    noFolds <- length(testFoldInds) # Number of folds
    noSamplesPerFolds <- sapply(testFoldInds, length)
    noFeatures <- length(featureNames(dataset)) # Total number of genes

    # Variable to be updated during the computation
    resCvPerModel <- alist()

    # --- Core of the method ---
    x <- exprs(dataset) # Levels of gene expression
    y <- pData(dataset)$type # class labels

    # Format the data to be usable by the pamr package
    mydata <- list( x=x,
                    y=factor(y),
                    geneid=featureNames(dataset),
                    genenames=featureNames(dataset))
                    
    # Number of samples from each class
    noSamplesPerClass <- tapply(classes,1:2, function(x) sum(y==x) )                    

    # /!\ thresholds must be in ascending order
    thresholdedTrain <-  pamr.train(data=mydata, threshold=optionValues)
    # Real number of threshold used (Some of them are not used if 0 genes
    # are selected)
    noModels <- length(thresholdedTrain$threshold)
    resultCV <- pamr.cv(fit=thresholdedTrain, data=mydata, folds=testFoldInds)
    # Best size of subset
  	#no.genes <- resultCV$size[which.min(resultCV$error)]

    # -- Store the best no of genes
    # reverse 2 times to get the bigger threshold corresponding to the smallest
    # cv error rate
    resCvPerModel$bestOptionValue <- rev(resultCV$threshold)[which.min(rev(resultCV$error))]

    # -- Store the error rates
    resCvPerModel$cvErrorRates <- resultCV$error # CV error rate

    # Number of error for each model in each class
    noErrorPerClass <- tapply(classes, 1:2,
        function(class) apply(resultCV$yhat, 2, function(pred) sum(pred[which(y==class)]!=y[which(y==class)])))
    # Tranform into matrix
    noErrorPerClass <- sapply(noErrorPerClass, function(x) x)
    # Find the corrsponding error rates per class
    resCvPerModel$classErrorRates <- apply(noErrorPerClass,1,'/', noSamplesPerClass)
    rownames(resCvPerModel$classErrorRates) <- classes

    # no samples per folds
    resCvPerModel$noSamplesPerFold <- noSamplesPerFolds

    # Error rate in each fold
    resCvPerModel$errorRatesPerFold <- t( sapply(resultCV$folds,
        function(f){
            rev(apply( resultCV$yhat[f,], 2,
                function(pred) sum(pred!=y[f])/length(f)) )
        } ))

    # Standard error
    resCvPerModel$seErrorRates <- as.numeric(sqrt(apply(resCvPerModel$errorRatesPerFold, 2, var)/noFolds))

    # Genes selected for each size of subset and each fold
    resCvPerModel$selectedGenes <- alist()
    # Values used to rank the genes for each gene in each fold for
    # each subset (RFE scrore)

    for (i in 1:noFolds){
        # /!\ Necessary to be able to get the names of the genes out of pamr.listgenes
        resultCV$cv.objects[[i]]$gene.subset <- 1:noFeatures
    }
    #  -- Store the selected genes
    for (i in 1:noModels){
        genesPerOneModel <- lapply(resultCV$cv.objects,
                function(train) {
                    # /!\ The nsc algo stops at the first nonzero = 0 and
                    # does not try the following thresholds
                    if (length(train$nonzero)>=i && train$nonzero[i] > 0){
                        genes <- pamr.listgenes(train, mydata, threshold=resultCV$threshold[i])
                        # Keep the maximum value for any class
                        genesWithoutNames <-  matrix(as.numeric(genes[,2:dim(genes)[2]]),
                                    dim<-c( nrow(genes),
                                            ncol(genes)-1))
                        scores <- t(genesWithoutNames)
                        # Copy the names of the genes from the first column
                        colnames(scores) <- genes[,1]
                        rownames(scores) <- colnames(genes)[2:length(colnames(genes))]
                    } else {
                        scores <- numeric(0)
                    }
                    return(scores)
                })

        resCvPerModel$selectedGenes[[i]] <- new("selectedGenesPerOneOption",
                                        genesPerOneModel,
                                        optionValue=optionValues[i])
    }
    currIndex <- 1

    return(resCvPerModel)
}

#-------------------------------------------------------------------------------
# Perform a cv per model defined by <optionValues> using the
# RFE-SVM on <dataset> with the folds <testFoldInds> and return
# the corresponding error rates for each model and the best no of genes:
# -->  Perform a CV: For each fold
#      and for each size of subset :
#       * Perform an RFE on the training data to find the best subset of
#         genes of a given size
#       * Train with the best subset of genes of a given size and test
#         on the test data
#
# @param  dataset (ExpressionSet)Data on which the classification is done
#         testFoldInds (list(numeric)) List of fold indices
#         optionValues (numeric)  Vector containing the thresholds to be
#                                 considered during feature selection
#         classes (vector(character)) Name of the classes
#         verbose (logical)       Debug variable
#         kernel (character) Type of kernel for the SVM
#
# @return resCvPerModel (list) Result of cv per model:
#               resCvPerModel$bestOptionValue: Value of the best option
#               resCvPerModel$cvErrorRates: cv error rate for each size of
#                   subset
#               resCvPerModel$classErrorRates: class error rate for each size of
#                   subset
#               resCvPerModel$noSamplesPerFold: No of samples in each fold
#               resCvPerModel$errorRatesPerFold: Error rate in each fold
#               resCvPerModel$seErrorRates: standard error on cv error rate
#               resCvPerModel$selectedGenes: genes selected for each size of
#                   subset
#-------------------------------------------------------------------------------
performACvPerModelWithSVMAndRFE  <- function( dataset,
                                testFoldInds,
                                foldCreation,
                                optionValues,
                                classes,
                                verbose=FALSE,
                                kernel="linear") {
    noFolds <- length(testFoldInds)
    resCvPerModel <- alist()
    # Genes selected per model and per fold and scores of feature selection
    resCvPerModel$selectedGenes <- alist()

    # Number of samples
    noSamples <- length(sampleNames(dataset))

    # Cross-validated recursive feature elimination:
    # -->  Perform a CV: For each fold
    #      and for each size of subset :
    #       * Perform an RFE on the training data to find the best subset of
    #         genes of a given size
    #       * Train with the best subset of genes of a given size and test
    #         on the test data
    model <- rfe.cv(  x = aperm(exprs(dataset)),
                      y = pData(dataset)$type,
                      #nfold = noFolds,
                        foldInds = testFoldInds,
                      noSelectedFeatures = optionValues,
                      foldCreation = foldCreation,
                      kern = kernel,
                      verbose=verbose)

    # Best size of subset
  	no.genes <- model$noFeatures[which.min(model$error.cv)]

    # -- Store the best no of genes
    resCvPerModel$bestOptionValue <- no.genes

    # -- Store the error rates
    resCvPerModel$cvErrorRates <- model$error.cv # CV error rate
    resCvPerModel$seErrorRates <- model$error.se # standard error on CV error rate
    resCvPerModel$classErrorRates <- model$error.ind # class error rate
    resCvPerModel$noSamplesPerFold <- model$noSample.perfold # no samples per folds

    # The models are ordered from the biggest to the smallest subset,
    # we reverse the order
    if (dim(model$error.rate.perfold)[2] > 1){
        # Reverse the order of the models (gene subsets)
        resCvPerModel$errorRatesPerFold <-  model$error.rate.perfold[,dim(model$error.rate.perfold)[2]:1]
    } else {
        resCvPerModel$errorRatesPerFold <-  model$error.rate.perfold
    }

    # Genes selected for each size of subset and each fold
    selectedGenes <- alist()
    # Values used to rank the genes for each gene in each fold for
    # each subset (RFE scrore)
    selectedGenes <- rev(model$rankingCriterion)

    selectedGenes <- lapply(selectedGenes,
        function(listFolds){
            namedListFolds <- lapply(listFolds,
                function(listGenes) {
                    matGenes <- matrix(listGenes, dim<-c(1, length(listGenes)))
                    rownames(matGenes) <- "rfe-scores"
                    colnames(matGenes) <- featureNames(dataset)[as.numeric(names(listGenes))]
                    return(matGenes)
                })
            return(namedListFolds)
        })

    #  -- Store the selected genes
    for (i in 1:length(selectedGenes)){
        resCvPerModel$selectedGenes[[i]] <- new("selectedGenesPerOneOption",
                selectedGenes[[i]],
                optionValue=optionValues[i])
    }
    return(resCvPerModel)
}

# Code useful for lda or naiveBayes, unused in current version ################
###############################################################################
#performACvPerModel <- function( dataset,
#                                #noFolds,
#                                testFoldInds,
#                                foldCreation,
#                                optionValues,
#                                classifierName,
#                                classes,
#                                #noTopGenes,
#                                verbose=FALSE,
#                                kernel="linear",
#                                featureSelectionFun,
#                                trainTestFun) {
#    noFolds <- length(testFoldInds)
#    resCvPerModel <- alist()
#    # Genes selected per model and per fold and scores of feature selection
#    resCvPerModel$selectedGenes <- alist()
#    noSamples <- length(sampleNames(dataset)) # Number of samples
#    noModels <- length(optionValues)# Number of sizes of genes subsets
#
#    # Labels of the output classes
#    presentClasses <- factor(levels(dataset$type))
#    noClasses <- length(presentClasses)
#
#    # -- Variables to be updated during the computation
#    # Store the standard error on CV error rate
#    resCvPerModel$seErrorRates <- vector("numeric", noModels)
#    # Store the error rate for each fold and each size of subset
#    resCvPerModel$errorRatesPerFold <-  matrix(nrow=noFolds, ncol=noModels)
#    # Indices of the genes selected
#    selectedGenes <- alist()
#    # Score of the genes selected during feature selection
#    orderingCriterion <- alist()
#    resCvPerModel$classErrorRates <- alist()
#    for ( i in 1:noClasses){
#        resCvPerModel$classErrorRates[[i]] <- matrix(nrow=noFolds, ncol=noModels)
#    }
#
#    # Indices of the folds
#    # Get the training indices from the test indices
#    foldInds <- lapply(testFoldInds, setdiff, x=seq(1:noSamples))
#
#    # For each fold
#    for (j in 1:noFolds) {
#        if (verbose){
#            print(paste("--- Fold ", j, " ---"))
#        }
#        # ------ Perform an RFE on the training dataset to find the best subset of
#        #         genes of a given size ------
#        # For each fold and each number of genes to be selected, compute a RFE
#        resultFeatSelection <- featureSelectionFun(
#                                            dataset=dataset[,foldInds[[j]]],
#                                            subsetsSizes=optionValues,
#                                            kernel=kernel,
#                                            verbose=verbose)
#
#        # Can be < to noModels with nsc when the threshold which leads to 0 gene
#        # is found before the end
#        practicalNoOfOptions <- length(resultFeatSelection$relevantGenes)
#
#        # For each size of subset
#        for (k in 1:practicalNoOfOptions){
#            if (j==1){
#                orderingCriterion[[k]] <- alist()
#            }
#
#            # Get the genes chosen by feature selection for the current size
#            # of subset
#            relevantGenes <- resultFeatSelection$relevantGenes[[k]]
#            if (length(relevantGenes) > 0){
#                # ----- Train with the best subset of genes of a given size and test
#                #         on the test data -----
#                errors <- trainTestFun( classifierName=classifierName,
#                                        dataset=dataset,
#                                        relevantGenes=relevantGenes,
#                                        trainInd=foldInds[[j]],
#                                        classes=presentClasses,
#                                        kernel=kernel,
#                                        bestThreshold=resultFeatSelection$bestThreshold)
#                # Error rate for the current fold and model
#                resCvPerModel$errorRatesPerFold[j,k] <- errors$errorRate
#                for (p in 1:noClasses){
#                    resCvPerModel$classErrorRates[[p]][j,k] <- errors$classErrorRates[p]
#                }
#                # -- Store the genes selected and their score during feture selection
#                #orderingCriterion[[k]][[j]] <- (rev(resultFeatSelection$TestedModels)[[k]])$orderingCriterion
#                # TODO #5
#                orderingCriterion[[k]][[j]] <- resultFeatSelection$scores[[k]]
#
#                ## Here
#                if (j == noFolds){
#                    resCvPerModel$selectedGenes[[k]] <- new("selectedGenesPerOneOption",
#                                        orderingCriterion[[k]],
#                                        optionValue=optionValues[k])
#                }
#            }
#        }
#    }
#
#    # -- Store the error rates
#    # CV error rate
#    #resCvPerModel$cvErrorRates <- as.numeric(lapply(intConfuMatrices, findErrorRate, classes))
#    noSamplesPerFold <- noSamples - as.numeric(lapply(foldInds, length))
#
#    resCvPerModel$cvErrorRates <- apply(resCvPerModel$errorRatesPerFold, 2, weighted.mean, w=noSamplesPerFold/noSamples)
#    resCvPerModel$seErrorRates <- sqrt(apply(resCvPerModel$errorRatesPerFold, 2, var)/noFolds)
#    # Class error rate
#    resCvPerModel$classErrorRates <- t(matrix(unlist(lapply(resCvPerModel$classErrorRates, function(x) apply(x,2, weighted.mean, w=noSamplesPerFold/noSamples))),ncol=noClasses))
#    rownames(resCvPerModel$classErrorRates) <- classes
#
#    # -- Store the best no of genes
#    resCvPerModel$bestOptionValue <- optionValues[which.min(resCvPerModel$cvErrorRates)]
#    resCvPerModel$noSamplesPerFold <- noSamplesPerFold
#  return(resCvPerModel)
#}
#
# Code useful for lda or naiveBayes, unused in current version ################
###############################################################################

#-------------------------------------------------------------------------------
# Feature Selection with RFE-SVM
#
# @param    dataset (ExpressionSet) microarray data
#           subsetsSizes (vector(numeric)) Sizes of subsets to be tried
#           kernel (character) Kernel for the SVM
#           verbose (logical) debug variable
#-------------------------------------------------------------------------------
RFESVMFeatureSelection <- function(dataset, subsetsSizes, kernel, verbose) {
    # Perform a RFE to find the best genes for each size of subset
    resultRFE <- rfe.fit( x=aperm(exprs(dataset)),
                            y=pData(dataset)[[1]],
                            #minf=min(noSelectedFeatures[1]),
                            # possible bug
                            minf=subsetsSizes[1],
                            noSelectedFeatures=subsetsSizes,
                            kern=kernel,
                            verbose=verbose)
    scores <- lapply(rev(resultRFE$TestedModels), '[[', "orderingCriterion")
    scores <- lapply(scores, function(scorePerOneModel){
                x <- matrix(scorePerOneModel, dim<-c(1, length(scorePerOneModel)))
                rownames(x) <- "rfe-score"
                colnames(x) <- featureNames(dataset[as.numeric(names(scorePerOneModel)),])
                return(x)
            })
    relevantGenes <- lapply(rev(resultRFE$TestedModels), '[[', "modelFeatures")
    return(list(relevantGenes=relevantGenes,
                scores=scores))
}

#-------------------------------------------------------------------------------
# Feature Selection with NSC
#
# @param    dataset (ExpressionSet) microarray data
#           subsetsSizes (vector(numeric)) thresholds to be tried
#           kernel (character) Kernel for the SVM
#           verbose (logical) debug variable
#-------------------------------------------------------------------------------
NSCFeatureSelection <- function(dataset, subsetsSizes, kernel, verbose) {
    # Train the NSC to find the best genes for each size of subset
    mydata <- list( x=exprs(dataset),
                    y=pData(dataset)$type,
                    geneid=featureNames(dataset),
                    genenames=featureNames(dataset))

    trainedNSC <- pamr.train(data=mydata, threshold=subsetsSizes)

    #  -- Store the selected genes
    genesPerOneModel <-tapply(trainedNSC$threshold, 1:length(trainedNSC$threshold),
            function(thres) {
                if (trainedNSC$nonzero[which(trainedNSC$threshold==thres)] > 0){
                    genes <- pamr.listgenes(trainedNSC, mydata, threshold=thres)
                    if (is.matrix(genes)){
                        noColumns <- ncol(genes) - 1
                        noRows <- nrow(genes)
                    } else {
                        # Case in which we have one gene only
                        noColumns <- length(genes) - 1
                        noRows <- 1
                    }
                    genesWithoutNames <-  matrix(as.numeric(genes[,2:dim(genes)[2]]),
                                    dim<-c( noRows, noColumns))
                    scores <- t(genesWithoutNames)
                    # Copy the names of the genes from the first column
                    colnames(scores) <- genes[,1]
                    rownames(scores) <- colnames(genes)[2:length(colnames(genes))]
                } else {
                    scores <- numeric(0)
                }
                return(scores)
            } )
    return(list(relevantGenes=lapply(genesPerOneModel, colnames),
                scores=genesPerOneModel))
}

###############################################################################
# Code useful to compute lda or NaiveBayes, not used in the current version of
# our package
###############################################################################
#MLITrainTest <- function(classifierName, dataset, relevantGenes, trainInd,  classes, kernel, bestThreshold){
#    # --- Variables ---
#    noClasses <- length(classes)
#    noErrorsPerClass <- vector("numeric", noClasses)
#    # --- Core of the function ---
#    switch( classifierName,
#        lda = intCvInfo <- MLearn(type~., dataset[relevantGenes,], ldaI, trainInd=trainInd),
#        naiveBayes = intCvInfo <- MLearn(type~., dataset[relevantGenes,], naiveBayesI, trainInd=trainInd),
#        stop(paste("In 'performACvPerModel' classifier ", classifierName, " is not supported"))
#    )
#    confuMatrix <- confuMat(intCvInfo)
#    for (i in 1:noClasses){
#        noErrorsPerClass[i] <- sum(confuMatrix[i,upper.tri(confuMatrix)[i,]], confuMatrix[i,lower.tri(confuMatrix)[i,]], na.rm=TRUE)
#    }
#
#    # Error rates calculation
#    return(errorRate=findErrorRate(confuMatrix, classes),
#            classErrorRates=findClassErrorRates(confuMatrix, classes),
#            noErrorsPerClass=noErrorsPerClass)
#}

#-------------------------------------------------------------------------------
# Training and testing with SVM-RFE
#
# @param    classifierName (character) Name of the classifier (ignored)
#           dataset (ExpressionSet) microarray data set
#           relevantGenes (vector(character)): Genes previously selected by RFE
#           trainInd (list(numeric)) Indices of the training fold
#           classes (vector(character)) name of the classes
#           kernel (character) kernel of the SVM
#           bestThreshold (numeric) Best size of subset (ignored)
#
# @return list:
#           $errorRate: error rate
#           $classErrorRates: class error rates
#           $noErrorsPerClass: No of errors in each class
#-------------------------------------------------------------------------------
SVMTrainTest <- function(classifierName, dataset, relevantGenes,
            trainInd, classes, kernel, bestThreshold){
    # --- Variables ---
    noClasses <- length(classes) # no of classes

    # Variables to be updated during the computation
    noSamplesPerClass <- vector("numeric") # no of samples per class
    noErrorsPerClass <- vector("numeric") # no of errors per class

    # --- Core of the function ---
    # Train the support vector machine on the training data considering only the
    # relevant genes
    fmodel <- svm(  aperm(exprs(dataset[relevantGenes,trainInd])),
                      pData(dataset[,trainInd])[[1]],
                      kernel=kernel)

    # Test the support vector machine on the test data considering only the
    # relevant genes
    if (is.matrix(aperm(exprs(dataset))[-trainInd,])==TRUE){
      out<-predict(fmodel,aperm(exprs(dataset))[-trainInd,relevantGenes])
    } else{
      out<-predict(fmodel,t(aperm(exprs(dataset))[-trainInd,relevantGenes]))
    }
    # Error rates calculation                     
    realOutput <- pData(dataset)[[1]][-trainInd]
    comparison <- out!= realOutput
    for (i in 1:noClasses)  {
        noErrorsPerClass[i] <- sum(comparison[realOutput==classes[i]], na.rm=TRUE)
        noSamplesPerClass[i] <- sum(realOutput==classes[i])
    }
    noTestSamples <- length(sampleNames(dataset)) - length(trainInd)
    return( list(errorRate=sum(noErrorsPerClass)/noTestSamples,
            classErrorRates=noErrorsPerClass/noSamplesPerClass,
            noErrorsPerClass=noErrorsPerClass))
}

#-------------------------------------------------------------------------------
# Training and testing with SVM-RFE
#
# @param    classifierName (character) Name of the classifier (ignored)
#           dataset (ExpressionSet) microarray data set
#           relevantGenes (vector(character)): Genes previously selected by NSc
#           trainInd (list(numeric)) Indices of the training fold
#           classes (vector(character)) name of the classes
#           kernel (character) kernel of the SVM (ignored)
#           bestThreshold (numeric) Best threshold
#
# @return list:
#           $errorRate: error rate
#           $classErrorRates: class error rates
#           $noErrorsPerClass: No of errors in each class
#-------------------------------------------------------------------------------
NSCTrainTest <- function(classifierName, dataset, relevantGenes, trainInd,
                        classes, kernel, bestThreshold){
    # --- Variables ---
    noClasses <- length(classes) # no of classes
    # Variables to be updated during the computation
    noSamplesPerClass <- vector("numeric") # no of samples per class
    noErrorsPerClass <- vector("numeric") # no of errors per class

    # --- Core of the function ---
    # Train the nearest shrunken centroid on the training data with the
    # best threshold
    trainData <- list( x=exprs(dataset[,trainInd]),
                    y=pData(dataset)$type[trainInd],
                    geneid=featureNames(dataset),
                    genenames=featureNames(dataset))  

    trainedNSC <- pamr.train(data=trainData, threshold=bestThreshold)
    # Test the nsc on the test data considering only the
    # relevant genes
    # threshold of zero, we want to keep all relevant genes
    prediction <- pamr.predict(fit=trainedNSC, exprs(dataset[,-trainInd]), threshold=bestThreshold)

    # Error rates calculation
    realOutput <- pData(dataset)[[1]][-trainInd]
    comparison <- prediction!= realOutput
    for (i in 1:noClasses)  {                                   
        noErrorsPerClass[i] <- sum(comparison[realOutput==classes[i]], na.rm=TRUE)
        noSamplesPerClass[i] <- sum(realOutput==classes[i])
    }
    noTestSamples <- length(sampleNames(dataset)) - length(trainInd)

    return( list(errorRate=sum(noErrorsPerClass)/noTestSamples,
            classErrorRates=noErrorsPerClass/noSamplesPerClass,
            noErrorsPerClass=noErrorsPerClass))
}

#-------------------------------------------------------------------------------
# Perform a cv the classifier <classifierName> on <dataset> with <trainIndices>
# as the indices of the training data.
#
# @param    dataset (ExpressionSet) microarray data set
#           trainIndices (list(numeric)) Indices of the training fold
#           noGenes: size of subset or threshold condidered
#           classes (vector(character)) name of the classes
#           classifierName (character) name of the classifier
#           kernel (character) kernel of the SVM (ignored)
#           featureSelectionFun (function) function to be used for feature
#               selection
#           trainTestFun (function) function to be used for training/testing
#
# @return list:
#           $errorRate: error rate
#           $classErrorRates: class error rates
#           $noErrorsPerClass: No of errors in each class
#           $relevantGenes: relevant genes sected by RFE or NSC
#-------------------------------------------------------------------------------
performAValidation <- function( dataset, trainIndices, noGenes, classes,
                                classifierName, verbose=FALSE, kernel="linear",
                                featureSelectionFun, trainTestFun) {
    # Number of classes
    noClasses <- length(classes)

    # Number of errors for each class
    noErrorsPerClass <- vector("numeric", noClasses)

    resFeatSelection <- featureSelectionFun(dataset=dataset[,trainIndices],
                                            subsetsSizes=noGenes,
                                            kernel=kernel,
                                            verbose=verbose)                                            
    genes <- resFeatSelection$relevantGenes[[1]]

    # rankingCriterion carry the information given by genes (names of
    # the genes) and the score of each of this gene during feature selection
    #   --> print(paste("genes", paste(genes, collapse=" ")))
    #   --> print(rankingCriterion)
    rankingCriterion <- resFeatSelection$scores[[1]]

    errors <- trainTestFun( classifierName=classifierName,
                            dataset=dataset,
                            relevantGenes=genes,
                            trainInd=trainIndices,
                            kernel=kernel,
                            classes=classes,
                            bestThreshold=noGenes[1])


     return(list(errorRate=errors$errorRate, classErrorRates=errors$classErrorRates,
                noErrorsPerClass=errors$noErrorsPerClass,
                relevantGenes=rankingCriterion))
}



#-------------------------------------------------------------------------------
# Return the frequency of all genes in <selectedGenesPerFold>
#
# @param  selectedGenesPerFold (list(vector("character")) List of genes selected
#                   for each fold. selectedGenesPerFold[[i]]
#                   corresponds to the genes selected in the ith fold
#                                 
# @return topGenes (frequencyGenes) 
#                   topGenes[[i]]@frequ contains the ith frequency
#                   topGenes[[i]]@geneList contains the corresponding genes
#-------------------------------------------------------------------------------
frequencyTopGenes <- function(selectedGenesPerFold) {
    topGenes <- alist()
    currListInd <- 1

    noGenes <- length(rownames(selectedGenesPerFold[[1]]))
    noTrials <- length(selectedGenesPerFold)

    genes <- character()
    noFolds <- length(selectedGenesPerFold)
    for (i in 1:noFolds){
        genes <- c(genes, colnames(selectedGenesPerFold[[i]]))
    }

    tableGenes <- table(genes)
    for (i in (noTrials*noGenes):1){
        if (length(labels(which(tableGenes==i))) != 0){
            topGenes[[currListInd]] <- new("frequencyGenes",
            frequ=i/noTrials,
            genesList=labels(which(tableGenes==i)))
            currListInd <- currListInd + 1
        }
    }
    return(topGenes)
}

#-------------------------------------------------------------------------------
# Return the indices of <noFolds> balanced testing folds according to the
# output <outputClasses>
#
# @param  outputClasses (factor) Vector of outputs
#         noFolds (numeric) Number of folds to be created
#-------------------------------------------------------------------------------
divideTestFolds <- function(outputClasses, noFolds) {
  # Order the outputCalsses to create balanced bins
  orderedIndices <- order(outputClasses)

  noSamples <- length(outputClasses)
  noBins <- floor(noSamples/noFolds)  # Number of bins to be created
  noSamplesPerBin <- noFolds*noBins   # Min number of samples per bin

  # Random assignment of the element of each bin to each fold
  foldallocs <- vector()
  for (j in 1:noBins) {
  	foldallocs <- c(foldallocs,sample(1:noFolds, replace=F))
  }

  # Consider the remaining indices
  numLeft = noSamples-noSamplesPerBin
  foldallocs <- c(foldallocs,sample(1:noFolds,numLeft, replace=F))

  # The indices of the folds correspond to the ordered indices (orderedIndices),
  # splitted according to their assignment by foldallocs
  foldsIndices <- split(orderedIndices,foldallocs)

  return(foldsIndices)
}

#-------------------------------------------------------------------------------
# Return the indices of <noFolds> naive testing folds according to the
# output <outputClasses>
#
# @param  noSamples (numeric) Number of samples
#         noFolds (numeric) Number of folds to be created
#-------------------------------------------------------------------------------
divideNaiveTestFolds <- function(noSamples, noFolds){
   # --- Test of the parameters ----
   if ( ! is.numeric(noSamples) || noSamples<=0 ){
      stop("Invalid parameter: noSamples must be a positive integer")
   }
   if ( ! is.numeric(noFolds) || noFolds<=0 ){
      stop("Invalid parameter: noFolds must be a positive integer")
   }
   if (noFolds != trunc(noFolds)){
      warning("Value truncated: noFolds truncated")
   }
   if (noSamples != trunc(noSamples)){
      warning("Value truncated: noFolds truncated")
   }
   if ( noFolds<2 || noFolds>noSamples ){
      stop("Invalid value: noFolds must greater than 2 and smaller than noSamples")
   }

   foldInds <- alist()
   groups <- alist()
   sizeFold <- floor(noSamples/noFolds)
   rest <- vector()
   for (i in 1:sizeFold){
     groups[[i]] <- seq(((i-1)*noFolds+1), i*noFolds)
   }

   for (i in 1:noFolds){
      foldInds[[i]] <- vector()
      for (j in 1:sizeFold){
        foldInds[[i]] <- c(foldInds[[i]], groups[[j]][i])
      }
   }
   if (sizeFold!=noSamples/noFolds){

      rest <- seq(sizeFold*noFolds+1,noSamples)
      for (i in 1:length(rest)){
         foldInds[[i]] <- c(foldInds[[i]], rest[i])
      }
   }

   return (foldInds)
}

#-------------------------------------------------------------------------------
# Return the indices of <noFolds> balanced training folds according to the
# output <outputClasses>
#
# @param  outputClasses (factor) Vector of outputs
#         noFolds (numeric) Number of folds to be created
# @return trainFoldInds (list) List of folds indices trainFoldInds[[i]]
#                       contains the indices for folds i
#-------------------------------------------------------------------------------
divideTrainFolds <- function(outputClasses, noFolds)
{
  # --- Test of the parameters ----
  if ( ! is.factor(outputClasses)){
    stop("Invalid parameter: outputClasses must be a vector of factors")
  }
  if ( ! is.numeric(noFolds) || noFolds<=0 ){
    stop("Invalid parameter: noFolds must be a positive integer")
  }
  if (noFolds != trunc(noFolds)){
    warning("Value truncated: noFolds truncated")
  }
  noSamples <- length(outputClasses)

  if ( noFolds<2 || noFolds>noSamples ){
    stop("Invalid value: noFolds must greater than 2 and smaller than the number of samples in outputClasses")
  }
  trainFoldInds <- alist()
  testFoldInds <- divideTestFolds(outputClasses, noFolds)
  for (i in 1:length(testFoldInds))
  {
    trainFoldInds[[i]] <- setdiff(1:noSamples, testFoldInds[[i]])
  }
  return (trainFoldInds)
}

#-------------------------------------------------------------------------------
# Return the indices of <noFolds> neive testing folds according to the
# output <outputClasses>
#
# @param  noSamples (numeric) Number of samples
#         noFolds (numeric) Number of folds to be created
#-------------------------------------------------------------------------------
divideNaiveTrainFolds <- function(noSamples, noFolds)
{
  if ( ! is.numeric(noFolds) || noFolds<=0 ){
    stop("Invalid parameter: noFolds must be a positive integer")
  }
  if (noFolds != trunc(noFolds)){
    warning("Value truncated: noFolds truncated")
  }

  if ( noFolds<2 || noFolds>noSamples ){
    stop("Invalid value: noFolds must greater than 2 and smaller than the number of samples in outputClasses")
  }
  trainFoldInds <- alist()
  testFoldInds <- divideNaiveTestFolds(noSamples, noFolds)
  for (i in 1:length(testFoldInds))
  {
    trainFoldInds[[i]] <- setdiff(1:noSamples, testFoldInds[[i]])
  }
  return (trainFoldInds)
}

# ------------------------------------------------------------------------------
# Divide a data set into <noFolds> according to the <typeFold> method.
# If <train> is TRUE the indices of the training folds are returned,
# if FALSE, the indices of the testing folds are returned
#
# @param    noFolds (numeric): Number of folds to be created
#           typeFold (character): Type of fold creation (balanced, naive
#                                   original)
#           train (logical): TRUE if the function must return training indices
#                            FALSE for testing indices
#           outputClasses (character): vector of class labels corresponding
#                                   to the samples
# ------------------------------------------------------------------------------
divideInFolds <- function(noFolds, typeFold, train, outputClasses) {
    # Indices of the different folds
    foldInds <- alist()
    # Number of samples
    noSamples <- length(outputClasses)
    if (typeFold=="original"){
        leave.out <- trunc(noSamples/noFolds)
      	o <- sample(1:noSamples)
      	foldTestInds <- vector("list", noFolds)
      	for (j in 1:(noFolds - 1)) {
      	    jj <- (1 + (j - 1) * leave.out)
      	    foldTestInds[[j]] <- (o[jj:(jj + leave.out - 1)])
      	}
      	foldTestInds[[noFolds]] <- o[(1 + (noFolds - 1) * leave.out):noSamples]
    } else{
      if (typeFold=="naive"){
        foldTestInds <- divideNaiveTestFolds(noSamples, noFolds)
      }
      else{
        foldTestInds <- divideTestFolds(outputClasses, noFolds)
      }
    }

    if (train) {
        # If we want the training indices
        foldInds <- lapply(foldTestInds, setdiff, x=seq(1,noSamples))
    } else {
        foldInds <- foldTestInds
    }
    return(foldInds)
}

# Code useful for lda or naiveBayes, unused in current version ################
###############################################################################
##-------------------------------------------------------------------------------
## Merge the confusion matrices stored in <listConfuMatrices> to create a new
## confusion matrix which summarise their results
##
## @param  listConfuMatrices (list) List of confusion matrices to be merged
## @return newConfuMat () New confusion matrix
##-------------------------------------------------------------------------------
#findConfuMat <- function(listConfuMatrices, classes){
#    if (! is.list(listConfuMatrices)){
#        extConfuMat <- listConfuMatrices
#        listConfuMatrices <- alist()
#        listConfuMatrices[[1]] <- extConfuMat
#    }
#    for ( i in 1:length(listConfuMatrices) ){
#        if (is.null(dimnames(listConfuMatrices[[i]]))){
#          dimnames(listConfuMatrices[[i]]) <- list(given=classes, predicted=classes)
#        }
#    }
#
#    newConfuMat <- matrix(0, nrow=nlevels(classes), ncol=nlevels(classes))
#    dimnames(newConfuMat) <- list(given=classes, predicted=classes)
#
#    for (i in 1:length(listConfuMatrices)){
#        # The resulting confusion matrix corresponds to the sum of the other matrices
#        # Handle the case when one or more classes are not present in the predicted classes
#        newConfuMat[match(dimnames(listConfuMatrices[[i]])$given,classes),match(dimnames(listConfuMatrices[[i]])$predicted,classes)] <- newConfuMat[match(dimnames(listConfuMatrices[[i]])$given,classes),match(dimnames(listConfuMatrices[[i]])$predicted,classes)] + listConfuMatrices[[i]]
#    }
#
#    return(newConfuMat)
#}
#
##-------------------------------------------------------------------------------
## Find the error rate corresponding to the confusion matrix <confuMatrix>
##
## @param  confuMatrix () Confusion matrix to be analysed
## @return (numeric) Error rate
##-------------------------------------------------------------------------------
#findErrorRate <- function(confuMatrix, classes){
#    # If it's not a square matrix, convert it to a square one
#    if (! all(dim(confuMatrix) == length(classes))){
#        confuMatrix <- findConfuMat(confuMatrix, classes)
#    }
#    noErrors <- sum(confuMatrix[upper.tri(confuMatrix)], confuMatrix[lower.tri(confuMatrix)])
#    return( noErrors / sum(confuMatrix) )
#}
#
##-------------------------------------------------------------------------------
## Find the individual classes error rates corresponding to the confusion matrix <confuMatrix>
##
## @param  confuMatrix () Confusion matrix to be analysed
## @return (numeric) Error rate
##-------------------------------------------------------------------------------
#findClassErrorRates <- function(confuMatrix, classes){
#    noClasses <- length(classes)
#    classErrorRates <- vector("numeric", noClasses)
#    # If it's not a square matrix, convert it to a square one
#    if (! all(dim(confuMatrix) == noClasses)){
#        confuMatrix <- findConfuMat(confuMatrix, classes)
#    }
#    for (i in 1:noClasses){
#        noErrors <- sum(confuMatrix[i,upper.tri(confuMatrix)[i,]], confuMatrix[i,lower.tri(confuMatrix)[i,]])
#        noSamplesInClass <- sum(confuMatrix[i,])
#        if (noSamplesInClass != 0){
#            classErrorRates[i] <- noErrors/sum(confuMatrix[i,])
#        } else {
#            classErrorRates[i] <- 0
#        }
#    }
#    return( classErrorRates )
#}
# Code useful for lda or naiveBayes, unused in current version ################
###############################################################################