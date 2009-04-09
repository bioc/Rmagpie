# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#--------------------------- non_linear_svm --------------------------------
# Function to hamdle non-linear kernel with Support Vector Machines.
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Calculate and return the criterion that will be used to order the genes in the
# RFE using the kernel <kernel>.
#
# @param    model (SVM model) svm model fitted to the features
#           kernel (character) type of kernel in "linear", "radial"
#           or "polynomial"
#
# @return   vector of numeric containing the mark of each feature
#-------------------------------------------------------------------------------
findDJ <- function(model, kernel){
    # Linear kernel
    if (kernel=="linear"){
        # Much faster !! than using K = vanilladot (?)
        DJ = 0.5*getW(model)^2
        return(DJ)
    } else {
        # Radial kernel
        if (kernel == "radial"){
            K <- rbfdot(sigma = model$gamma)
        }
        else {
            # Polynomial kernel
            if (kernel == "polynomial"){
                K <- polydot(degree = model$degree)
            }
        }
    }
    matK <- alist()
    nbSV <- dim(model$SV)[1]
    nbFeatures <- dim(model$SV)[2]
    DJ <- vector("numeric", nbFeatures)

    for (i in 1:nbFeatures){
    matK[[i]] <- (kernelMatrix(K, as.matrix(model$SV[,i])))
    DJ[i] <- t(matrix(model$coefs))%*%matK[[i]]%*%matrix(model$coefs)
    }

    # It works but it's too long
    #  for (k in 1:nbSV){
    #    for (l in 1:nbSV)
    #      for (i in 1:nbFeatures){
    #        DJ[i] <- DJ[i] + 0.5*model$coefs[k]*model$coefs[l]*K(gamm=model$gamma, u=model$SV[k,i], v=model$SV[l,i])
    #      }
    #  }
    return(DJ)
}

#-------------------------------------------------------------------------------
# Return w as stated in the original RFE paper from the SVM model.
#
# @param    model (SVM model) svm model fitted to the features
# @return   numeric
#-------------------------------------------------------------------------------
getW <- function(model){
  nbFeatures <- dim(model$SV)[2]
	w <- vector("numeric", nbFeatures)
  nbSV <- dim(model$SV)[1]
  for (k in 1:nbSV){
    w <- w + model$coefs[k]*model$SV[k,]
  }
  return(w)
}

#-------------------------------------------------------------------------------
# Return alpha as stated in the original RFE paper from the SVM model.
#
# @param    model (SVM model) svm model fitted to the features
# @return   numeric
#-------------------------------------------------------------------------------
getAlpha <- function(model){
  alpha <- abs(model$coefs)
  return(alpha)
}

#-------------------------------------------------------------------------------
# Return y as stated in the original RFE paper from the SVM model.
#
# @param    model (SVM model) svm model fitted to the features
# @return   numeric
#-------------------------------------------------------------------------------
getY <- function(model){
  y <- sign(model$coefs)
  return(y)
}
