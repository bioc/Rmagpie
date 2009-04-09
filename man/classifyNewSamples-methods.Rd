\name{classifyNewSamples-methods}
\docType{methods}
\alias{classifyNewSamples}
\alias{classifyNewSamples-methods}
\alias{classifyNewSamples,assessment-method}
\title{ classifyNewSamples Method to classify new samples for a given assessment}
\description{
This method classify one or several new samples provided in the file 'newSamplesFile'
using the final classifier build by 'findFinalClassifier'.

}
\section{Methods}{
\describe{

\item{object = "assessment"}{This method is only applicable on objects of class
        assessment.}
}}

\arguments{
  \item{object}{\code{object of class assessment}. Object assessment of interest}
  \item{newSamplesFile}{\code{character}. URL of the file containing the
    gene expressions of the samples to be classified. The first line of the file
    must corresponds to the sample names and the first column to the names of the genes.}
  \item{optionValue}{\code{numeric}. Size of subset (for RFE-SVM) or threshold (for NSC)
    to be considered, the option value  must be available in the slot
    \code{featureSelectionOptions} of the assessment. If not, the smallest
    value bigger than 'optionValue' is selected. If this argument is missing the
    best option value according to one-layer cross-validation is used.}
}
\examples{

data('vV70genesDataset')

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

\dontrun{
classifyNewSamples(expeOfInterest, "pathToFile/testSamples_geneExpr.txt", 4)
}
expeOfInterest <- runOneLayerExtCV(expeOfInterest)
\dontrun{
classifyNewSamples(expeOfInterest, "pathToFile/testSamples_geneExpr.txt")
}
}


\keyword{methods}
