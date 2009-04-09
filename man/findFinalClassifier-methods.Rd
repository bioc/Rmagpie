\name{findFinalClassifier-methods}
\docType{methods}
\alias{findFinalClassifier}
\alias{findFinalClassifier-methods}
\alias{findFinalClassifier,assessment-method}
\title{findFinalClassifier Method to train and build the final classifier based on an assessment}
\description{
This method generates and stores the final classifier corresponding
to an assessment. This classifier can then be used to classify new samples by
calling \code{classifyNewSamples}. The final classifier is build according to the
classifier selected for a given assessment, applied on the whole data considering
only the genes selected by the feature selction method selected.
}
\section{Methods}{
\describe{

\item{object = "assessment"}{This method is only applicable on objects of class
        assessment.}
}}

\value{
    The methods returns an object of class assessment which finalClassifier has been build.
}

\seealso{
    \code{\linkS4class{finalClassifier}}, \code{\linkS4class{assessment}}
}
\examples{
#dataPath <- file.path("C:", "Documents and Settings", "c.maumet", "My Documents", "Programmation", "Sources", "SVN", "R package", "data")
#aDataset <- new("dataset", dataId="vantVeer_70", dataPath=dataPath)
#aDataset <- loadData(aDataset)
data('vV70genesDataset')

# With the RFE-SVM as feature selection method
expeOfInterest <- new("assessment", dataset=vV70genes,
                                   noFolds1stLayer=10,
                                   noFolds2ndLayer=9,
                                   classifierName="svm",
                                   typeFoldCreation="original",
                                   svmKernel="linear",
                                   noOfRepeat=2,
                                   featureSelectionOptions=new("geneSubsets", optionValues=c(1,2,4,8,16,32,64,70)))

# Build the final classifier
expeOfInterest <- findFinalClassifier(expeOfInterest)

# With the NSC as feature selection method
expeOfInterest <- new("assessment", dataset=vV70genes,
                                   noFolds1stLayer=10,
                                   noFolds2ndLayer=9,
                                   featureSelectionMethod="nsc",
                                   classifierName="nsc",
                                   typeFoldCreation="original",
                                   svmKernel="linear",
                                   noOfRepeat=2,
                                   featureSelectionOptions=new("thresholds"))

# Build the final classifier
expeOfInterest <- findFinalClassifier(expeOfInterest)
}

\keyword{methods}
