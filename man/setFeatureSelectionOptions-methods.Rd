\name{setFeatureSelectionOptions-methods}
\docType{methods}
\alias{getFeatureSelectionOptions<-,assessment-method}
\alias{getFeatureSelectionOptions<--methods}
\alias{getFeatureSelectionOptions<-}
\title{getFeatureSelectionOptions<- Method to modify the attributes of a featureSelectionOptions from an assessment}
\description{
This method provides an easy interface to modify the attributes of the object of class
featureSelectionOptions related to a particular assessment, directly from this object assessment.
The argument \code{topic} specifies which part of the featureSelectionOptions is of interest.
This method is only available none of the one-layer
CV or two-layers CV have been performed and the final classifier has not been determined yet.
}

\section{Methods}{
\describe{
\item{object = "assessment"}{The method is only applicable on objects of class
        assessment.}
}}

\arguments{
  \item{object}{\code{Object of class assessment}. Object assessment of interest}
  \item{topic}{\code{character}. Optional argument that specifies which attribute of
    the featureSelectionOptions must be replaced, the possible values are:
    \code{"optionValues"} (slot \code{optionValues} of the featureSelectionOptions),
    \code{"noOfOptions"} (slot \code{noOfOptions} of the featureSelectionOptions),
    if the featureSelectionOptions object is an object of class \code{geneSubsets}, then
    the following values are also available for the argument \code{topic} 
    \code{"subsetsSizes"} (slot \code{optionValues} of the geneSubsets),
    \code{"noModels"} (slot \code{noOfOptions} of the geneSubsets),
    \code{"maxSubsetSize"} (slot \code{maxSelectedFeatures} of the geneSubsets),
    \code{"speed"} (slot \code{speed} of the featureSelectionOptions),
    if the featureSelectionOptions object is an object of class \code{thresholds}, then
    the following values are also available for the argument \code{topic} 
    \code{"thresholds"} (slot \code{optionValues} of the object thresholds),
    \code{"noThresholds"} (slot \code{noOfOptions} of the object thresholds)

    if the \code{topic} is missing then the whole featureSelectionOptions object is replaced.}
}

\value{
    The methods modifies the object of class assessment and returned the slot modified
    accordingly to the request provided by \code{topic}.

    If \code{topic} is missing
    \code{object of class featureSelectionOptions} featureSelectionOptions corresponding to the assessment is replaced by \code{value}.

    If \code{topic} is \code{"optionValues"}
    \code{numeric} Slot \code{optionValues} of the featureSelectionOptions is replaced by \code{value}.

    If \code{topic} is "noOfOptions"
    \code{numeric} Slot \code{noOfOptions} of the featureSelectionOptions is replaced by \code{value}.

    If \code{object} is of class \code{geneSubsets} and \code{topic} is "maxSubsetSize"
    \code{numeric} Slot \code{maxSubsetSize} of the geneSubsets is replaced by \code{value}.
    
    If \code{object} is of class \code{geneSubsets} and \code{topic} is "subsetsSizes"
    \code{numeric} Slot \code{optionValues} of the geneSubsets is replaced by \code{value}.
    
    If \code{object} is of class \code{geneSubsets} and \code{topic} is "noModels"
    \code{numeric}Slot \code{noOfOptions}  of the geneSubsets is replaced by \code{value}.
    
    If \code{object} is of class \code{geneSubsets} and \code{topic} is "speed"
    \code{numeric} Slot \code{speed}  of the geneSubsets is replaced by \code{value}.

    If \code{object} is of class \code{thresholds} and \code{topic} is "thresholds"
    \code{numeric} Slot \code{optionValues} of the object of class thresholds is replaced by \code{value}.
    
    If \code{object} is of class \code{thresholds} and \code{topic} is "noThresholds"
    \code{numeric} Slot \code{noOfOptions} of the object of class thresholds is replaced by \code{value}.

}
\author{ Camille Maumet }

\seealso{
    \code{\linkS4class{featureSelectionOptions}}, \code{\linkS4class{assessment}}
}
\examples{

# With an assessment using RFE
data('vV70genesDataset')

mySubsets <- new("geneSubsets", optionValues=c(1,2,3,4,5,6))
myExpe <- new("assessment", dataset=vV70genes,
                                   noFolds1stLayer=10,
                                   noFolds2ndLayer=9,
                                   classifierName="svm",
                                   typeFoldCreation="original",
                                   svmKernel="linear",
                                   noOfRepeat=2,
                                   featureSelectionOptions=mySubsets)

# Modify the size of the biggest subset
getFeatureSelectionOptions(myExpe, topic='maxSubsetSize') <- 70
getFeatureSelectionOptions(myExpe, topic='maxSubsetSize')
# Modify all the sizes of subsets
getFeatureSelectionOptions(myExpe, topic='subsetsSizes') <- c(1,5,10,25,30)
getFeatureSelectionOptions(myExpe, topic='subsetsSizes')
# Modify the speed
getFeatureSelectionOptions(myExpe, topic='speed') <- 'slow'
getFeatureSelectionOptions(myExpe, topic='speed')
# Modify the entire geneSubsets
getFeatureSelectionOptions(myExpe) <- mySubsets
getFeatureSelectionOptions(myExpe, topic='maxSubsetSize')
getFeatureSelectionOptions(myExpe, topic='subsetsSizes')
getFeatureSelectionOptions(myExpe, topic='speed')
getFeatureSelectionOptions(myExpe, topic='noModels')


# With an assessment using NSC as a feature selection method
myThresholds <- new("thresholds", optionValues=c(0.1,0.2,0.3))
myExpe2 <- new("assessment", dataset=vV70genes,
                                   noFolds1stLayer=10,
                                   noFolds2ndLayer=9,
                                   classifierName="nsc",
                                   featureSelectionMethod='nsc',
                                   typeFoldCreation="original",
                                   svmKernel="linear",
                                   noOfRepeat=2,
                                   featureSelectionOptions=myThresholds)

otherThresholds <- new("thresholds", optionValues=c(0,0.5,1,1.5,2,2.5,3))
# Modify the whole object 'featureSelectionOptions' (an object of class thresholds)
getFeatureSelectionOptions(myExpe2) <- otherThresholds
getFeatureSelectionOptions(myExpe2, topic='thresholds')
getFeatureSelectionOptions(myExpe2, topic='noThresholds')
}

\keyword{methods}