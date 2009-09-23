\name{getDataset-methods}
\docType{methods}
\alias{getDataset-methods}
\alias{getDataset}
\alias{getDataset,assessment-method}
\title{getDataset Method to access the attributes of a dataset from an assessment}
\description{
This method provides an easy interface to access the attributes of a dataset
directly from an object assessment. The argument \code{topic} specifies which part
of the dataset is of interest.
}

\section{Methods}{
\describe{
\item{object = "assessment"}{The method is only applicable on objects of class
        assessment.}
}}

\arguments{
  \item{object}{\code{Object of class assessment}. Object assessment of interest}
  \item{topic}{\code{character}. Optional argument that specifies which attribute of
    the dataset is requested, the possible values are
    \code{"dataId"} (slot dataId of the dataset),
    \code{"dataPath"} (slot dataPath of the dataset),
    \code{"geneExprFile"} (slot geneExprFile of the dataset),
    \code{"classesFile"} (slot classesFile of the dataset),
    \code{"eset"} (slot eset of the dataset)
    if the \code{"topic"} is missing then the whole dataset object is returned.}
}

\value{
    The value returned by the method changes accordingly to the \code{"topic"} argument.

    If \code{"topic"} is missing
    \code{object of class dataset} the dataset corresponding to the assessment of interest

    If \code{"topic"} is \code{"dataId"}
    \code{object of class character} corresponding to the \code{dataId} of the dataset

    If \code{"topic"} is \code{"dataPath"}
    \code{object of class character} corresponding to the \code{dataPath} of the dataset

    If \code{"topic"} is \code{"geneExprFile"}
    \code{object of class character} corresponding to the \code{geneExprFile} of the dataset

    If \code{"topic"} is \code{"classesFile"}
    \code{object of class character} corresponding to the \code{classesFile} of the dataset

     If \code{"topic"} is \code{"eset"}
    \code{object of class ExpressionSetOrNull} corresponding to the \code{eset} of the dataset
}
\author{ Camille Maumet }

\seealso{
    \code{\linkS4class{assessment}}
}
\examples{
#dataPath <- file.path("C:", "Documents and Settings", "c.maumet", "My Documents", "Programmation", "Sources", "SVN", "R package", "data")
#aDataset <- new("dataset", dataId="vantVeer_70", dataPath=dataPath)
#aDataset <- loadData(aDataset)

data('vV70genesDataset')

expeOfInterest <- new("assessment", dataset=vV70genes,
                                   noFolds1stLayer=10,
                                   noFolds2ndLayer=9,
                                   classifierName="svm",
                                   typeFoldCreation="original",
                                   svmKernel="linear",
                                   noOfRepeat=2,
                                   featureSelectionOptions=new("geneSubsets", optionValues=c(1,2,3,4,5,6)))

getDataset(expeOfInterest)

}


\keyword{methods}
