# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#----------------------------- user_interface ----------------------------------
# This file provides the function and methods needed to plot one-layer and
# two-layer CV
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot the results of two-layer CV
#
# @param    object (assessment) assessment of interest
#-------------------------------------------------------------------------------
setMethod("plotErrorsFoldTwoLayerCV", "assessment", function(object){
    foldErrorRates <- matrix(nrow=getNoOfRepeats(object), ncol=getNoFolds2ndLayer(object))
    bestOptionValue <- matrix(nrow=getNoOfRepeats(object), ncol=getNoFolds2ndLayer(object))
    for (i in 1:getNoOfRepeats(object)){
        foldErrorRates[i,] <- getResults(object, c(2,i), topic='errorRate', errorType='fold')
        for (j in 1:getNoFolds2ndLayer(object)){
            bestOptionValue[i,j] <- getResults(object, c(2,i,j), topic='bestOptionValue')
        }
    }
    if (is(getFeatureSelectionOptions(object), "geneSubsets")) {
        msgTitle <- "Error rate per fold in the 2nd Layer of CV  vs subset sizes"
        if (getFeatureSelectionOptions(object, 'speed') == 'high'){
            optionValues <- log2(getFeatureSelectionOptions(object, 'optionValues'))
            bestOptionValue <- log2(bestOptionValue)
            msgX <- "log2(Number of genes)"
        } else {
            optionValues <- getFeatureSelectionOptions(object, 'optionValues')
            msgX <- "Number of genes"
        }
    } else {
        msgTitle <- "Error rate per fold in the 2nd Layer of CV vs value of thresholds"
        optionValues <- rev(getFeatureSelectionOptions(object, 'optionValues'))
        msgX <- "Thresholds"
    }
    msgY <- "Error rate per fold in the second layer of CV"
    plotErrorRatePerFold( foldErrorRates, bestOptionValue, optionValues,
                            noRepeats=getNoOfRepeats(object), msgX, msgY, msgTitle)
})

#-------------------------------------------------------------------------------
# Plot the error rate for each fold in two-layer CV
#
# @param    foldErrorRates (numeric) Error rate in each fold
#           bestOptionValue (numeric) best option value
#           optionValues (numeric) value of options
#           noRepeats (numeric) No of repeats
#           msgX (character) Message for x axis
#           msgY (character) Message for y axis
#           msgTitle (character) Title for the graph
#-------------------------------------------------------------------------------
plotErrorRatePerFold <- function(foldErrorRates, bestOptionValue, optionValues,
                                noRepeats, msgX="", msgY="", msgTitle="") {

    minError <- min(foldErrorRates)
    maxError <- max(foldErrorRates)# + seErrorRate)

    legendWidth <- 0.45
    legendHeight <- 0.11*noRepeats

    posXLegend <- min(optionValues)#max(foldId) - (max(foldId) - min(foldId))/2
    positions <- findPositionsLegendWithFold(bestOptionValue, foldErrorRates, optionValues, maxError, legendWidth, legendHeight)

    posYLegend <- positions$posYLegend
    maxY <- positions$maxY

    yRange <- c(0, maxY)
    xRange <- c(min(optionValues), max(optionValues))
    plot(bestOptionValue,foldErrorRates,
        main=msgTitle,
        type="p", col=3, pch=19,
        ylim=yRange, xlim=xRange,
        lty=1,
        xlab=msgX, ylab=msgY)

    # Plot a new line for each repeat
    legendText <- paste("repeat", 1)
    legendCol <- c(3)
    if (noRepeats >= 2){
        for (i in 2:noRepeats){
            points(bestOptionValue[i,], foldErrorRates[i,],
            type="p", col=i+2, pch=19)
            legendText <- c(legendText, paste("repeat", i))
            legendCol <- c(legendCol, i+2)
        }
    }

    title(msgTitle)

    # Legend
    legend(x=posXLegend, y=posYLegend,
          legendText,
          merge=TRUE, pch=c(19),
          col=legendCol, lty=c(1))
}

#-------------------------------------------------------------------------------
# Find the position of the legend
#
# @param    foldErrorRates (numeric) Error rate in each fold
#           bestOptionValue (numeric) best option value
#           optionValues (numeric) value of options
#           maxError (numeric) Maximum value of the y axis
#           legendWidth (numeric) Approximative width of the legend
#           legendHeight (numeric) Approximative height of the legend
#-------------------------------------------------------------------------------
findPositionsLegendWithFold <- function(bestOptionValue, foldErrorRates,
            optionValues, maxError, legendWidth, legendHeight){
    # Last feature that could be covered by the legend
    lastPosition <- max(optionValues)*legendWidth
    # Position of the points that could be covered by the legend
    lastPosition <- which(bestOptionValue<=lastPosition)

    # Minimum error value of the points that could be covered by the legend
    minCovered <- min(foldErrorRates[lastPosition])
    # Maximum error value of the points that could be covered by the legend
    maxCovered <- max(foldErrorRates[lastPosition])
    # If there is enough room the legend goes at the bottom of the graph
    necessarySpace <- maxError*legendHeight
    if (minCovered >= necessarySpace){
        # Legend at the bottom
        posYLegend <- necessarySpace
        maxY <- maxError
    } else {
        # Legend on the top
        topPosition <- max(maxError, max(foldErrorRates[lastPosition])+necessarySpace)
        posYLegend <- topPosition
        maxY <- max(topPosition, maxError)
    }
    return(list(maxY=maxY, posYLegend=posYLegend))
}

#-------------------------------------------------------------------------------
# Find the position of the legend
#
# @param    nbSelFeatures (numeric) Values on the x axis
#           cvErrorRate (numeric) Values on the y axis
#           maxError (numeric) Maximum value of the y axis
#           legendWidth (numeric) Approximative width of the legend
#           legendHeight (numeric) Approximative height of the legend
#-------------------------------------------------------------------------------
findPositionsLegend <- function(nbSelFeatures,  cvErrorRate,
            seErrorRate=0, maxError, legendWidth, legendHeight){
    # Last feature that could be covered by the legend
    lastPosition <- max(nbSelFeatures)*legendWidth
    # Position of the last feature that could be covered by the legend
    lastPosition <- which(nbSelFeatures>lastPosition)[1]

    # Minimum error value of the points that could be covered by the legend
    minCovered <- min(cvErrorRate[1:lastPosition]-seErrorRate[1:lastPosition])
    # Maximum error value of the points that could be covered by the legend
    maxCovered <- max(cvErrorRate[1:lastPosition]+seErrorRate[1:lastPosition])
    # If there is enough room the legend goes at the bottom of the graph
    necessarySpace <- maxError*legendHeight
    if (minCovered >= necessarySpace){
        # Legend at the bottom
        posYLegend <- necessarySpace
        maxY <- maxError
    } else {
        # Legend at the top
        topPosition <- max(maxError, max(cvErrorRate[1:lastPosition]+seErrorRate[1:lastPosition])+necessarySpace)
        posYLegend <- topPosition
        maxY <- max(topPosition, maxError)
    }
    return(list(maxY=maxY, posYLegend=posYLegend))
}

#-------------------------------------------------------------------------------
# Plot the results of one-layer CV (summary)
#
# @param    object (assessment) assessment of interest
#-------------------------------------------------------------------------------
setMethod("plotErrorsSummaryOneLayerCV", "assessment", function(object){
    if (is(getFeatureSelectionOptions(object), "geneSubsets")) {
        msgTitle <- "Ext. CV Error rate vs number of selected genes"
        if (getFeatureSelectionOptions(object, 'speed') == 'high'){
            optionValues <- log2(getFeatureSelectionOptions(object, 'optionValues'))
            msgX <- "log2(Number of genes)"
            cvErrorRate <- getResults(object, 1, 'errorRate', errorType='cv')
            seErrorRate <- getResults(object, 1, 'errorRate', errorType='se')
        } else {
            optionValues <- getFeatureSelectionOptions(object, 'optionValues')
            msgX <- "Number of genes"
            cvErrorRate <- getResults(object, 1, 'errorRate', errorType='cv')
            seErrorRate <- getResults(object, 1, 'errorRate', errorType='se')
        }
    } else {
        msgTitle <- "Ext. CV Error rate vs value of thresholds"
        optionValues <- rev(getFeatureSelectionOptions(object, 'optionValues'))
        msgX <- "Thresholds"
        cvErrorRate <- rev(getResults(object, 1, 'errorRate', errorType='cv'))
        seErrorRate <- rev(getResults(object, 1, 'errorRate', errorType='se'))
    }
    msgY <- "Error rate with external cross validation"
    plotErrors( cvErrorRate=cvErrorRate,
                nbSelFeatures=optionValues,
                seErrorRate=seErrorRate,
                msgX=msgX,
                msgY=msgY,
                msgTitle=msgTitle) })

#-------------------------------------------------------------------------------
# Plot the results of one-layer CV (summary)
#
# @param    cvErrorRate (numeric) Error rate for each option value
#           seErrorRate (numeric) standard error
#           nbSelFeatures (numeric) value of options
#           msgX (character) Message for x axis
#           msgY (character) Message for y axis
#           msgTitle (character) Title for the graph
#-------------------------------------------------------------------------------
plotErrors <- function(cvErrorRate, nbSelFeatures=seq(along=cvErrorRate),
    seErrorRate=0, msgX="", msgY="", msgTitle=""){
    minError <- min(cvErrorRate)
    maxError <- max(cvErrorRate + seErrorRate)

    # Legend covers 45% of the graph in width and 25% in height
    legendWidth <- .45
    legendHeight <- .20

    posXLegend <- min(nbSelFeatures)

    positions <- findPositionsLegend(nbSelFeatures, cvErrorRate, seErrorRate, maxError, legendWidth, legendHeight)

    posYLegend <- positions$posYLegend
    maxY <- positions$maxY

    yRange <- c(0, maxY)
    xRange <- c(min(nbSelFeatures), max(nbSelFeatures))

    plot(nbSelFeatures, cvErrorRate,
        type="o", col="black", pch=19,
        ylim=yRange, xlim=xRange,
        lty=1,
        xlab=msgX, ylab=msgY)
    if (!all(seErrorRate==0))    {
        superpose.eb(nbSelFeatures, cvErrorRate, seErrorRate);
    }

    # Minumum in red
    # Warning, do not use which.min instead, here we want all the minima
    xMins <- nbSelFeatures[which(cvErrorRate==min(cvErrorRate))]
    yMins <- rep(min(cvErrorRate), times=length(xMins))
    points(x=xMins, y=yMins, pch=15, col="red")

    arrowLength <- .1
    arrowTextLength <- .13

    # Scale of the x axis
    scaleLengthX <- xRange[2] - xRange[1]
    # Scale of the y axis
    scaleLengthY <- yRange[2] - yRange[1]

    # Find the position of the arrow pointing at the first minimum
    positionArrows <- findPositionArrow(xMins[1], yMins[1], nbSelFeatures, cvErrorRate,
                                        seErrorRate, arrowLength, arrowTextLength,
                                        scaleLengthX, scaleLengthY)
    arrows(positionArrows$xStartArrows, positionArrows$yStartArrows,
            positionArrows$xEndArrows, positionArrows$yEndArrows, length=.1, col='red')
    text(x=positionArrows$xStartTexts, y=positionArrows$yStartTexts, "Min. error (biased)", col="red", cex=0.8)

    # Title of the graph
    title(msgTitle)

    # Legend
    legend(x=posXLegend, y=posYLegend,
          c("Ext. cv error"),
          merge=TRUE, pch=c(19),
          col=c("black"), lty=c(1))
}
#------------------------------------------------------------------------------
# Find the possition where the arrow pointing at the first min of error
# must be located
#
# @param    xMins (numeric) Value along the x axis corresponding to min along y
#           yMins (numeric) Min along the y axis
#           nbSelFeatures (numeric) Values of options
#           cvErrorRate (numeric) CV error rate
#           seErrorRate (numeric) standard error
#           arrowLength (numeric) length of the arrow
#           arrowTextLength (numeric) length of the text going with the arrow
#           scaleLengthX (numeric) length of the x axis
#           scaleLengthY (numeric) length of the y axis
#------------------------------------------------------------------------------
findPositionArrow<- function(xMins, yMins, nbSelFeatures, cvErrorRate, seErrorRate,
                            arrowLength, arrowTextLength, scaleLengthX, scaleLengthY){
    goRight <- FALSE
    # end of the arrow points to the minimum of the error rate
    xEndArrow <- xMins
    yEndArrow <- yMins
    # By default we draw an arrow that goes from bottom to top and left to right
    goRight <- TRUE
    goTop <- TRUE

    # Conversion from % to unit of the x axis

    # width of the arrow and the text
    widthTextAndArrow <- max(nbSelFeatures)*(arrowLength + arrowTextLength)
    # width of the arrow alone
    widthArrow <- max(nbSelFeatures)*arrowLength

    # --- First trial: Just make an horizontal arrow ---
    # Where the text should start
    xStartAll <-  xEndArrow - widthTextAndArrow
    # Where the arrow should start
    xStartArrow <- xEndArrow - widthArrow
    yStartArrow <- yEndArrow # Horizontal arrow

    # - Test if the arrow is out of the plot on the left -
    # If the arrow is out, we draw it from right to left instead
    if (xStartAll < 0 ) {
        goRight <- FALSE
        xStartArrow <- xEndArrow + widthArrow
        xStartAll <-  xEndArrow + widthTextAndArrow
    }

    # - Test if this arrow cover part of the plot -

    # Position of the features that might be covered by the arrow or its legend
    if (goRight == TRUE){
        positions <- intersect(which(nbSelFeatures >= xStartAll),
                            which(nbSelFeatures < xEndArrow))
    } else {
        positions <- intersect(which(nbSelFeatures <= xStartAll),
                            which(nbSelFeatures > xEndArrow))
    }
    # Minimum of the error rate (considering the error bars) of the dots that
    # could be covered
    minThatCanBeCovered <- cvErrorRate[positions] - seErrorRate[positions]
    # Maximum of the error rate (considering the error bars) of the dots that
    # could be covered
    maxThatCanBeCovered <- cvErrorRate[positions] + seErrorRate[positions]

    # Which positions, the arrow is going through
    criticalPosition <- positions[which(minThatCanBeCovered <= yEndArrow &  maxThatCanBeCovered >= yEndArrow)]
    # If we have at least one dot of the plot that is covered by the arrow
    # --- Second trial: diagonal arrow from bottom to top ---
    if (length(criticalPosition) >= 1){
        # We let 5% of blank space between the graph and the arrow
        blankSpace <- .05
        # Conversion from % to unit of the y axis
        blankSpace <- max(cvErrorRate[positions] + seErrorRate[positions])*blankSpace
        yStartArrow <- min(cvErrorRate[criticalPosition] - seErrorRate[criticalPosition]) - blankSpace
        # - Test if the arrow is out of the plot at the bottom -
        # --- Third trial: diagonal arrow from top to bottom ---
        if (yStartArrow < 0){
            yStartArrow <- max(cvErrorRate[criticalPosition] + seErrorRate[criticalPosition]) + blankSpace
        }
    } else {
        # No problem we keep the horizontal arrow
    }

    # Position of the text
    if (goRight == TRUE){
        xStartText <- xStartArrow - scaleLengthX*(arrowTextLength)
    } else {
        xStartText <- xStartArrow + scaleLengthX*(arrowTextLength)
    }
    return(list(xStartArrows=xStartArrow, yStartArrows=yStartArrow,
                xEndArrows=xEndArrow, yEndArrows=yEndArrow,
                xStartTexts=xStartText, yStartTexts=yStartArrow))
}

#-------------------------------------------------------------------------------
# Plot the results of one-layer CV (all repeats)
#
# @param    object (assessment) assessment of interest
#-------------------------------------------------------------------------------
setMethod("plotErrorsRepeatedOneLayerCV", "assessment", function(object)
{
    if (is(getFeatureSelectionOptions(object), "geneSubsets")) {
        msgTitle <- "Ext. CV Error rate vs number of selected genes"
        if (getFeatureSelectionOptions(object, 'speed') == 'high'){
            optionValues <- log2(getFeatureSelectionOptions(object, 'optionValues'))
            msgX <- "log2(Number of genes)"
        } else {
            optionValues <- getFeatureSelectionOptions(object, 'optionValues')
            msgX <- "Number of genes"
        }
    } else {
        msgTitle <- "Ext. CV Error rate vs value of thresholds"
        optionValues <- getFeatureSelectionOptions(object, 'optionValues')
        msgX <- "Thresholds"
    }
    msgY <- "Error rate with external cross validation"

    originalCVs <- getResult1LayerCV(object)@original1LayerCV

    # Error rate retained after all the repeats
    summaryCVErr <- getResults(object, 1, topic='errorRate', errorType='cv')
    summarySEErr <- getResults(object, 1, topic='errorRate', errorType='se')
    noOriginalCV <- length(originalCVs)

    # Store the error rate of each original 1 Layer CV
    originalCVErr <- alist()

    for (i in 1:noOriginalCV){
        originalCVErr[[i]] <- as.numeric(originalCVs[[i]]@errorRates@cvErrorRate@cvErrorRate)
    }
    noOptionValues <- length(summaryCVErr)
    minError <- min(summaryCVErr, min(unlist(originalCVErr)))
    maxError <- max(summaryCVErr + summarySEErr, max(unlist(originalCVErr)))

    # Legend covers 45% of the graph in width and 25% in height
    legendWidth <- .45
    legendHeight <- .1+.1*noOriginalCV
    posXLegend <- min(optionValues)

    # Tip to get the real min and max of CV error over the repeats
    maximumDecalage <- vector("numeric", noOptionValues)
    for (i in 1:noOptionValues){
        maximumDecalage[i] <- max (max(as.numeric(lapply(originalCVErr, "[[", i))) - summaryCVErr[i],
                                   summaryCVErr[i] -  min(as.numeric(lapply(originalCVErr, "[[", i))),
                                   summarySEErr[i])
    }

    positions <- findPositionsLegend(optionValues, summaryCVErr, maximumDecalage, maxError, legendWidth, legendHeight)

    posYLegend <- positions$posYLegend
    maxY <- positions$maxY

    yRange <- c(0, maxY)
    xRange <- c(min(optionValues), max(optionValues))

    plot(optionValues, summaryCVErr,
        type="o", col="black", pch=19,
        ylim=yRange, xlim=xRange,
        lty=1,lwd=2,
        xlab=msgX, ylab=msgY)
    if (!all(summarySEErr==0))    {
        superpose.eb(optionValues, summaryCVErr, summarySEErr);    
    }

    # Plot a new line for each repeat
    legendText <- c("Ext. cv error")
    legendCol <- c("black")
    for (i in 1:noOriginalCV){
        points(optionValues, originalCVErr[[i]],
        type="o", col=i+2, pch=19)
        legendText <- c(legendText, paste("repeat", i))
        legendCol <- c(legendCol, i+2)
    }
    points(optionValues, summaryCVErr,
        type="o", col="black", pch=19,
        lty=1,lwd=2)
    if (!all(summarySEErr==0))    {
        superpose.eb(optionValues, summaryCVErr, summarySEErr);
    }

    # Minumum in red
    # Warning, do not use which.min instead, here we get all the minima
    xMins <- optionValues[which(summaryCVErr==min(summaryCVErr))]
    yMins <- rep(min(summaryCVErr), times=length(xMins))
    points(x=xMins, y=yMins, pch=15, col="red")

    arrowLength <- .1
    arrowTextLength <- .13

    scaleLengthX <- xRange[2] - xRange[1]
    scaleLengthY <- yRange[2] - yRange[1]
    positionArrows <- findPositionArrow(xMins[1], yMins[1], optionValues, summaryCVErr,
                                        maximumDecalage, arrowLength, arrowTextLength,
                                        scaleLengthX, scaleLengthY)
    arrows(positionArrows$xStartArrows[1], positionArrows$yStartArrows[1],
            positionArrows$xEndArrows[1], positionArrows$yEndArrows[1], length=.1, col='red')

    text(x=positionArrows$xStartTexts, y=positionArrows$yStartTexts, "Min. error (biased)", col="red", cex=0.8)
    title(msgTitle)

    # Legend
    legend(x=posXLegend, y=posYLegend,
          legendText,
          merge=TRUE, pch=c(19, rep(19, noOriginalCV)),
          col=legendCol, lty=c(1, rep(1, noOriginalCV)))
})

#-------------------------------------------------------------------------------
# Plot error bars from R blog by Raoul Grasman available at
# http://users.fmg.uva.nl/rgrasman/rpages/2005/09/error-bars-in-plots.html
#
# @param    x : X coordinate of point where we want to add an error bar
# @param    y : Y coordinate of point where we want to add an error bar
# @param    ebl : Half height of the bar from bootom to point coordinate
# @param    ebu : Half height of the bar from point to top (default = abl)
# @param    length : Length of the horizontal bars   
#-------------------------------------------------------------------------------
superpose.eb <-
function (x, y, ebl, ebu = ebl, length = 0.08, ...){
    arrows(x, y + ebu, x, y - ebl, angle = 90, code = 3,
    length = length, ...)    
}