# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#------------------------------------- rfe -------------------------------------
# Functions to perform an RFE.
#
# Author: Christophe Ambroise
#   Further modified: Justin Zhu, Camille Maumet
# Last Modified: 18 Jan. 2008
#-------------------------------------------------------------------------------

############################################################
# Compute the Standart deviation, the mean of the error, and the individual error
############################################################
crossval.errors<- function(y, cv.fit, nfold, groups)
{
  nomodels<-ncol(cv.fit)
  nc<- length(table(y))
  n<-length(y)          # number of observations
  ni <- rep(NA, nfold)  # number of observations in each fold

  # Fold error for each model: a column per model, a line per fold
  err <- matrix(NA, ncol = nomodels, nrow = nfold)
  temp <- matrix(y, ncol = nomodels, nrow = n)
  for (i in 1:nfold){
    ii<-groups[[i]]
    ni[i]<-length(ii)
    if (is.null(dim(temp[ii,]))){
      err[i, ] <- (temp[ii, ] != cv.fit[ii, ])/ni[i]
    } else {
      err[i, ] <- apply(temp[ii, ] != cv.fit[ii, ], 2, sum)/ni[i]
    }
  }

  error.se <- sqrt(apply(err, 2, var)/nfold)
  error.cv <- apply(err,2,weighted.mean,w=(ni/n))

  # Individual error per class: a column per model, a line per class
  error.ind<-matrix(NA, ncol = nomodels, nrow = nc)
  for (i in 1:nomodels){
    s <- table(y,cv.fit[,i])
    diag(s) <- 0
    error.ind[, i] <- apply(s, 1, sum)/table(y)
  }
  rownames(error.ind)<-dimnames(table(y))[[1]]

  return(list( error.se=error.se,
               error.cv=error.cv,
               error.ind=error.ind,
               error.rate.perfold=err,
               noSample.perfold=ni) )
}

#################################################################
# Calculates the primal variables w stored in warray
#################################################################
svm.weight<-function(model)
{
  # Calculates the primal variables w stored in warray
  # warray[k,l,] is the weight vector for the binary pb class k against class l
  nclass<-length(model$labels)
  classk<- rep(1:nclass,model$nSV)
  p<-dim(model$SV)[2]
  # array of the weight vectors
  warray<-array(0,dim<-c(nclass,nclass,p))
  # loop to use the coefs
  for (s in 1:dim(model$SV)[1])
  for (co in 1:(nclass-1)){
    # find the two class problem
    k<- classk[s]
    l <- ((1:nclass)[-k])[co]
    warray[k,l,]<-  warray[k,l,] + model$coefs[s,co]*model$SV[s,]
    warray[l,k,]<-  warray[l,k,] + model$coefs[s,co]*model$SV[s,]
  }
  return(warray)
}


#################################################################
# Return twice the sum of the absolute value of primal variables w
#################################################################
linear.feasel<-function(model)
{
  # return twice the sum of the absolute value of primal variables w
	sumabsw<-apply(abs(svm.weight(model)),3,sum)
}

noFvec<- function(p,speed="high")
{
  if (speed=="high"){
    nomodels<-ceiling(log2(p))+1
    noFeatures<-2^(0:(nomodels-1))
    noFeatures[nomodels]<-p
    return(noFeatures)
  } else {
    noFeatures<-(1:p)
  }
  return(noFeatures)

}
#####################################################
# orderFeatures runs a SVM and order the features
#####################################################
orderFeatures<- function(x, y, FeaturesFromBest2Worst, noFeatures, speed, noSelectedFeatures, kern)
{
    # Reorganize the <noFeatures> first features of the list <featuresFromBest2Worst>
    p <- dim(x)[2]
    if (noFeatures == p) {
        FeaturesToTest <- (1:p)
    } else {
        FeaturesToTest <- FeaturesFromBest2Worst[1:noFeatures]
    }

    model <- svm (x[,FeaturesToTest], y , kernel = kern)
    sumabsw <- findDJ(model=model, kernel=kern)
    FeatureOrder <- order(sumabsw,decreasing = TRUE)
    FeaturesFromBest2Worst[1:noFeatures] <- FeaturesToTest[FeatureOrder]

    # Keep only the weights of the genes that are kept for the next step
    #orderingCriterion <- sort(sumabsw, decreasing=TRUE)[1:noFeatures]
    # Keep all the weights
    orderingCriterion <- sort(sumabsw, decreasing=TRUE)
    names(orderingCriterion) <- FeaturesFromBest2Worst[1:noFeatures]

  if ( noSelectedFeatures[1] == 0){
    if (speed=="low"){
      # Remove Feature by Feature
      noFeatures <- noFeatures-1
    } else {
      # Remove Half of the Features at each iteration
      if (noFeatures==2^(floor(log2(noFeatures)))) { # if the actual number of features is a power of 2
        noFeatures <- max(1,noFeatures/2)
      } else {
          noFeatures<- 2^(floor(log2(noFeatures)))
      }
    }
  } else {
    index <- max(1, which(noSelectedFeatures == noFeatures)-1)
    noFeatures <- noSelectedFeatures[index]
  }


  return(list(  noFeatures=noFeatures,
                FeaturesFromBest2Worst=FeaturesFromBest2Worst,
                model=model,
                modelFeatures=FeaturesToTest,
                orderingCriterion=orderingCriterion))
}


######################################################
# RFE for multi-class SVM
#
# rfe functions for cross-validation
#######################################################
rfe.predict<-function(fit,x) {
    nomodels<-length(fit$TestedModels)

    if (! is.null(dim(x))){
      yhat<-data.frame(matrix(0,dim(x)[1],nomodels))
      for (j in 1:nomodels){
        yhat[,j]<-predict((fit$TestedModels[[j]])$model,x[,(fit$TestedModels[[j]])$modelFeatures])
      }
    } else {
      yhat<-data.frame(matrix(0,1,nomodels))
      for (j in 1:nomodels){
        yhat[,j]<-predict((fit$TestedModels[[j]])$model,t(matrix(x[(fit$TestedModels[[j]])$modelFeatures])))
      }
    }
    return(yhat)
}

########################################################
rfe.cv<-function(x, y, foldInds, speed="high", noSelectedFeatures, foldCreation, kern, verbose )
{
   p <- dim(x)[2]
   nfold <- length(foldInds)
    # Find the number of models to be tested
    if ( noSelectedFeatures[1] == 0) {
        if (speed=="high"){
            nomodels<- ceiling(log2(p))+1
        }
        else{
          nomodels<- p
        }
    } else {
        nomodels <- length(noSelectedFeatures)
    }

    # -->  Perform a CV: For each fold
    #      and for each size of subset :
    #       * Perform an RFE on the training data to find the best subset of
    #         genes of a given size
    #       * Train with the best subset of genes of a given size and test
    #         on the test data
    n <- length(y)
    model.fit <- vector("list", nfold)

    cv.fit <- data.frame(matrix(y, n, nomodels))
    # For each fold
    for (j in 1:nfold) {
        if (verbose){
            print(paste("--- Fold ", j, " ---"))
        }
        #  For each size of subset :
        #       * Perform an RFE on the training data to find the best subset of
        #         genes of a given size       
        model.fit[[j]] <- rfe.fit(x=x[-foldInds[[j]], ], y=y[-foldInds[[j]]], kern=kern, verbose=verbose)
        #       * Test on the test data
        cv.fit[foldInds[[j]],] <- rfe.predict(model.fit[[j]], x[foldInds[[j]], ])
    }

    errors <- crossval.errors(y,cv.fit, nfold, foldInds)
    # Present the results
    if (noSelectedFeatures[1] == 0){
        noFeatures <- noFvec(p,speed)
    } else{
        noFeatures <- noSelectedFeatures
    }
    Flist=matrix(0,p,nfold)       # matrix of the selected features for each fold

    selFeaturesPerModel <- alist()
    # Criterion use to order the feature during RFE
    rankingCriterion <- alist()
    selFeaturesPerFold <- alist()

    for (j in 1:nomodels){
        selFeaturesPerModel[[j]] <- matrix()
        selFeaturesPerFold[[j]] <- alist()
        rankingCriterion[[j]] <- alist()
        for (i in 1:nfold){
            if (i==1){
                length(selFeaturesPerModel[[j]]) <- nfold*length(model.fit[[i]]$TestedModels[[j]]$modelFeatures)
                dim(selFeaturesPerModel[[j]]) <- c(length(model.fit[[i]]$TestedModels[[j]]$modelFeatures), nfold)
            }
            Flist[,i]<- model.fit[[i]]$Flist
            selFeaturesPerModel[[j]][,i] <- model.fit[[i]]$TestedModels[[j]]$modelFeatures
            selFeaturesPerFold[[j]][[i]] <- model.fit[[i]]$TestedModels[[j]]$modelFeatures
            rankingCriterion[[j]][[i]] <- model.fit[[i]]$TestedModels[[j]]$orderingCriterion[as.character(model.fit[[i]]$TestedModels[[j]]$modelFeatures)]
        }
    }

    return(list(  error.cv = rev(errors$error.cv),
                error.se = rev(errors$error.se),
                error.ind = errors$error.ind[,nomodels:1],
                error.rate.perfold = errors$error.rate.perfold,
                noSample.perfold = errors$noSample.perfold,
                noFeatures = noFeatures,
                Flist = Flist,
                selectedFeaturesPerModel = selFeaturesPerModel,
                selectedFeaturesPerFold = selFeaturesPerFold,
                rankingCriterion=rankingCriterion))
}
########################################################
rfe.ae<-function(x,y,speed="high")
{
  # RFE apparent error rate estimation for different number of
  # selected features (from 1 to p)
  fit<-rfe.fit(x,y,speed=speed)
  yhat<-rfe.predict(fit,x)
  nomodels<-dim(yhat)[2]
  error.ae<- rev( colMeans(yhat != data.frame(matrix(rep(y,nomodels),length(y),nomodels))))
  noFeatures<- noFvec(dim(x)[2],speed)

  return(list(error.ae=error.ae,noFeatures=noFeatures,FList=fit$Flist))
}

# Recursive Feature Elimination for multiclass SVM until minf features are
# selected Returns a list with all the  selected models
rfe.fit<-function(x, y, minf=1, speed="high", noSelectedFeatures=0, kern, verbose) 
{
  if (all(y == y[1])){
    stop("samples from x can't all come from the same class")
  }
  p <- dim(x)[2]
  FeaturesFromBest2Worst <- seq(1,p)
  TestedModels <- list()

  if (noSelectedFeatures[1] == 0){
    noFeatures <- p
  } else {
    noFeatures <- noSelectedFeatures[length(noSelectedFeatures)]
    minf <- min(noSelectedFeatures)
  }

  # If the greatest number of features to be tested is not equal to the total number
  # of features we must compute a first step to get the <noFeatures> best genes
  if (noFeatures != p) {
    tmp <- orderFeatures( x,
                          y,
                          FeaturesFromBest2Worst,
                          noFeatures,
                          speed,
                          noSelectedFeatures=noSelectedFeatures,
                          kern)
    FeaturesFromBest2Worst <-tmp$FeaturesFromBest2Worst
  }

  i <- 1

  while (noFeatures > minf){
    if (verbose){
        print(paste("no. features=", noFeatures))
    }
    # Fit an SVM and remove the genes that are less relevant
    tmp <- orderFeatures( x,
                          y,
                          FeaturesFromBest2Worst,
                          noFeatures,
                          speed,
                          noSelectedFeatures=noSelectedFeatures,
                          kern)

    FeaturesFromBest2Worst <- tmp$FeaturesFromBest2Worst
    orderingCriterion <- tmp$orderingCriterion

    # keep track of the models
    TestedModels[[i]] <- list(  model=tmp$model,
                                modelFeatures=tmp$modelFeatures,
                                orderingCriterion=orderingCriterion[as.character(tmp$modelFeatures)])

    i <- i+1
    noFeatures <-  tmp$noFeatures
  }


    lastSelectedGenes <- FeaturesFromBest2Worst[1:noFeatures]
    noSelectedGenes <- length(lastSelectedGenes)

    # For the last subset size (smallest) we must compute the ranking criterion
    # by applying another time the svm
    model <- svm (x[,lastSelectedGenes], y , kernel = kern)
    sumabsw <- findDJ(model=model, kernel=kern)
    # Keep the weights
    orderingCriterion <- sort(sumabsw, decreasing=TRUE)
    names(orderingCriterion) <- FeaturesFromBest2Worst[1:noFeatures]

    TestedModels[[i]]<- list( model=svm(  x[,FeaturesFromBest2Worst[1:noFeatures]],
                                        y,
                                        kernel=kern),
                            modelFeatures=lastSelectedGenes,
                            orderingCriterion=orderingCriterion[1:minf])
    return(list(Flist=FeaturesFromBest2Worst,TestedModels=TestedModels))
}