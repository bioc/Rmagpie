\name{plotErrorsFoldTwoLayerCV-methods}
\docType{methods}
\alias{plotErrorsFoldTwoLayerCV}
\alias{plotErrorsFoldTwoLayerCV-methods}
\alias{plotErrorsFoldTwoLayerCV,assessment-method}
\title{ plotErrorsFoldTwoLayerCV  Method to plot the error rate of a two-layer Cross-validation}
\description{
 This method creates a plot that reprenset the error rate in each fold of each repeat of the
second layer of cross-validation of the two-layer cross-validation of the assessment at stake. The plot
represents the error rate versus the size of gene subsets (for SVM-RFE) or the threshold values (for NSC).
}
\section{Methods}{
\describe{
\item{object = "assessment"}{The method is only applicable on objects of class
        assessment.}
}}

\seealso{
    \code{\link{plotErrorsSummaryOneLayerCV-methods}}, \code{\link{plotErrorsRepeatedOneLayerCV-methods}}
}

\examples{
data('vV70genesDataset')

expeOfInterest <- new("assessment", dataset=vV70genes,
                                   noFolds1stLayer=3,
                                   noFolds2ndLayer=2,
                                   classifierName="svm",
                                   typeFoldCreation="original",
                                   svmKernel="linear",
                                   noOfRepeat=10,
                                   featureSelectionOptions=new("geneSubsets", optionValues=c(1,2,3,4,5,6)))

expeOfInterest <- runTwoLayerExtCV(expeOfInterest)

plotErrorsFoldTwoLayerCV(expeOfInterest)

}

\keyword{methods}
