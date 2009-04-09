\name{getFeatureSelectionOptions-methods}
\docType{methods}
\alias{getFeatureSelectionOptions,assessment-method}
\alias{getFeatureSelectionOptions-methods}
\alias{getFeatureSelectionOptions}
\title{getFeatureSelectionOptions Method to access the attributes of a featureSelectionOptions from an assessment}
\description{
This method provides an easy interface to access the attributes of the object of class
featureSelectionOptions related to a particular assessment, directly from this object assessment.
The argument \code{topic} specifies which part of the featureSelectionOptions is of interest.
}

\section{Methods}{
\describe{
\item{object = "assessment"}{The method is only applicable on objects of class
        assessment.}
}}

\arguments{
  \item{object}{\code{Object of class assessment}. Object assessment of interest}
  \item{topic}{\code{character}. Optional argument that specifies which attribute of
    the featureSelectionOptions is requested, the possible values are:
    \code{"optionValues"} (Access the slot \code{optionValues} of the featureSelectionOptions),
    \code{"noOfOptions"} (Access the slot \code{noOfOptions} of the featureSelectionOptions),
    if the featureSelectionOptions object is an object of class \code{geneSubsets}, then
    the following values are also available for the argument \code{topic} 
    \code{"subsetsSizes"} (Access the slot \code{optionValues} of the geneSubsets), 
    \code{"noModels"} (Access the slot \code{noOfOptions} of the geneSubsets),       
    \code{"maxSubsetSize"} (Access the slot \code{maxSelectedFeatures} of the geneSubsets),
    \code{"speed"} (Access the slot \code{speed} of the featureSelectionOptions),    
    if the featureSelectionOptions object is an object of class \code{thresholds}, then
    the following values are also available for the argument \code{topic} 
    \code{"thresholds"} (Access the slot \code{optionValues} of the object thresholds), 
    \code{"noThresholds"} (Access the slot \code{noOfOptions} of the object thresholds)

    if the \code{topic} is missing then the whole featureSelectionOptions object is returned.}
}

\value{
    The value returned by the method changes accordingly to the 'topic' argument.

    If \code{topic} is missing
    \code{object of class featureSelectionOptions} the featureSelectionOptions corresponding to the assessment of interest

    If \code{topic} is "optionValues"
    \code{numeric} corresponding to the \code{optionValues} of the featureSelectionOptions

    If \code{topic} is "noOfOptions"
    \code{numeric} corresponding to the \code{noOfOptions} of the featureSelectionOptions

    If \code{object} is of class \code{geneSubsets} and \code{topic} is "maxSubsetSize"
    \code{numeric} corresponding to the \code{maxSubsetSize}  of the geneSubsets
    
    If \code{object} is of class \code{geneSubsets} and \code{topic} is "subsetsSizes"
    \code{numeric} corresponding to the \code{optionValues} of the geneSubsets  
    
    If \code{object} is of class \code{geneSubsets} and \code{topic} is "noModels"
    \code{numeric} corresponding to the \code{noOfOptions}  of the geneSubsets
    
    If \code{object} is of class \code{geneSubsets} and \code{topic} is "speed"
    \code{numeric} corresponding to the \code{speed}  of the geneSubsets  

    If \code{object} is of class \code{thresholds} and \code{topic} is "thresholds"
    \code{numeric} corresponding to the \code{optionValues} of the object of class thresholds  
    
    If \code{object} is of class \code{thresholds} and \code{topic} is "noThresholds"
    \code{numeric} corresponding to the \code{noOfOptions} of the object of class thresholds 

}
\author{ Camille Maumet }

\seealso{
    \code{\linkS4class{featureSelectionOptions}}, \code{\linkS4class{assessment}}
}
\examples{
# With an assessment using RFE
#dataPath <- file.path("C:", "Documents and Settings", "c.maumet", "My Documents", "Programmation", "Sources", "SVN", "R package", "data")
#aDataset <- new("dataset", dataId="vantVeer_70", dataPath=dataPath)
#aDataset <- loadData(aDataset)

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

# Return the whole object 'featureSelectionOptions' (an object of class geneSusbsets)
getFeatureSelectionOptions(myExpe)
# Size of the biggest subset
getFeatureSelectionOptions(myExpe, topic='maxSubsetSize')
# All sizes of subsets
getFeatureSelectionOptions(myExpe, topic='subsetsSizes')
# Speed
getFeatureSelectionOptions(myExpe, topic='speed')
# Number of subsets
getFeatureSelectionOptions(myExpe, topic='noModels') == getNoModels(mySubsets)

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

# Return the whole object 'featureSelectionOptions' (an object of class geneSusbsets)
getFeatureSelectionOptions(myExpe2)
# vector of thresholds
getFeatureSelectionOptions(myExpe2, topic='thresholds')
# Number of thresholds
getFeatureSelectionOptions(myExpe2, topic='noThresholds')
}

\keyword{methods}