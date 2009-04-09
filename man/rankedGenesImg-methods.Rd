\name{rankedGenesImg-methods}
\docType{methods}
\alias{rankedGenesImg}
\alias{rankedGenesImg-methods}
\alias{rankedGenesImg,assessment-method}
\title{rankedGenesImg Method to plot the genes according to their frequency in a microarray like image}
\description{
 Generate an image per value of option representing the features (on dot per feature). The
color of the dot depends on the frequency of the feature in for the given value of option
(number of genes or threshold).
}
\section{Methods}{
\describe{
\item{object = "assessment"}{The method is only applicable on objects of class
        assessment.}
}}


\arguments{
  \item{object}{\code{Object of class assessment}. Object assessment of interest.}
  \item{storagePath}{\code{character}. URL where the image must be stored.}
}

\examples{

\dontrun{
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

rankedGenesImg(expeOfInterest, storagePath='myPath')
}

}

\keyword{methods}

