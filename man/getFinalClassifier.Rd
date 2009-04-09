\name{getFinalClassifier-methods}
\docType{methods}
\alias{getFinalClassifier,assessment-method}
\alias{getFinalClassifier-methods}
\alias{getFinalClassifier}
\title{getFinalClassifier Method to access the attributes of a finalClassifier from an assessment}
\description{
This method provides an easy interface to access the attributes of the object of class
finalClassifier related to a particular assessment, directly from this object assessment.
The argument \code{topic} specifies which part of the finalClassifier is of interest.
}

\section{Methods}{
\describe{
\item{object = "assessment"}{The method is only applicable on objects of class
        assessment.}
}}

\arguments{
  \item{object}{\code{Object of class assessment}. Object assessment of interest}
  \item{topic}{\code{character}. Optional argument that specifies which attribute of
    the finalClassifier is requested, the possible values are
    \code{genesFromBestToWorst} (slot \code{genesFromBestToWorst} of the finalClassifier),
    \code{models} (slot \code{models} of the finalClassifier),
    if the \code{topic} is missing then the whole finalClassifier object is returned.}
}

\value{
    The value returned by the method changes accordingly to the  \code{topic} argument.

    If  \code{topic} is missing
    \code{object of class finalClassifier} the finalClassifier corresponding to the assessment of interest

    If  \code{topic} is  \code{"genesFromBestToWorst"}
    \code{numeric} corresponding to the \code{genesFromBestToWorst} of the finalClassifier

    If  \code{topic} is  \code{"models"}
    \code{numeric} corresponding to the \code{models} of the finalClassifier

}
\author{ Camille Maumet }

\seealso{
    \code{\linkS4class{finalClassifier}}, \code{\linkS4class{assessment}}
}
\examples{
#dataPath <- file.path("C:", "Documents and Settings", "c.maumet", "My Documents", "Programmation", "Sources", "SVN", "R package", "data")
#aDataset <- new("dataset", dataId="vantVeer_70", dataPath=dataPath)
#aDataset <- loadData(aDataset)

mySubsets <- new("geneSubsets", optionValues=c(1,2,3,4,5,6))
data('vV70genesDataset')

# assessment with RFE and SVM
expeOfInterest <- new("assessment", dataset=vV70genes,
                               noFolds1stLayer=10,
                               noFolds2ndLayer=9,
                               classifierName="svm",
                               typeFoldCreation="original",
                               svmKernel="linear",
                               noOfRepeat=2,
                               featureSelectionOptions=mySubsets)

expeOfInterest <- findFinalClassifier(expeOfInterest)

# Return the whole object of class finalClassifier
getFinalClassifier(expeOfInterest)
getFinalClassifier(expeOfInterest, 'genesFromBestToWorst')
getFinalClassifier(expeOfInterest, 'models')

# assessment with NSC
expeOfInterest <- new("assessment", dataset=vV70genes,
                               noFolds1stLayer=10,
                               noFolds2ndLayer=9,
                               featureSelectionMethod='nsc',
                               classifierName="nsc",
                               typeFoldCreation="original",
                               svmKernel="linear",
                               noOfRepeat=2,
                               featureSelectionOptions=new("thresholds"))

expeOfInterest <- findFinalClassifier(expeOfInterest)

# Return the whole object of class finalClassifier
getFinalClassifier(expeOfInterest)
getFinalClassifier(expeOfInterest, 'genesFromBestToWorst')
getFinalClassifier(expeOfInterest, 'models')
}

\keyword{methods}