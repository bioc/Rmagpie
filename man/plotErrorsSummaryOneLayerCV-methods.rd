\name{plotErrorsSummaryOneLayerCV-methods}
\docType{methods}
\alias{plotErrorsSummaryOneLayerCV}
\alias{plotErrorsSummaryOneLayerCV-methods}
\alias{plotErrorsSummaryOneLayerCV,assessment-method}
\title{ plotErrorsSummaryOneLayerCV  Method to plot the summary estimated error rates of a one-layer Cross-validation}
\description{
 This method creates a plot that represent the summary estimated error rate of the
one-layer cross-validation of the assessment at stake. The plot
represents the summary estimate of the error rate (averaged over the repeats) versus the size
of gene subsets (for SVM-RFE) or the threshold values (for NSC).
}
\section{Methods}{
\describe{
\item{object = "assessment"}{The method is only applicable on objects of class
        assessment.}
}}

\seealso{
    \code{\link{plotErrorsFoldTwoLayerCV-methods}}, \code{\link{plotErrorsRepeatedOneLayerCV-methods}}
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

expeOfInterest <- runOneLayerExtCV(expeOfInterest)

plotErrorsSummaryOneLayerCV(expeOfInterest)

}

\keyword{methods}
