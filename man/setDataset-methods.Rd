\name{setDataset-methods}
\docType{methods}
\alias{getDataset<-}
\alias{getDataset<--methods}
\alias{getDataset<-,assessment-method}
\title{getDataset<- Method to modify the attributes of a dataset from an assessment}
\description{
This method provides an easy interface to modify the attributes of a dataset
directly from an object assessment. The argument topic specifies which part
of the dataset should be modified. This method is only available none of the one-layer
CV or two-layers CV have been performed and the final classifier has not been determined yet.
}

\section{Methods}{
\describe{
\item{object = "assessment"}{This method is only applicable on objects of class
        assessment.}
}}

\arguments{
  \item{object}{class assessment. Object assessment of interest}
  \item{topic}{character. Optional argument that specifies which attribute of
    the dataset must be changed, the possible values are
    \code{dataId} (slot dataId of the dataset),
    \code{dataPath} (slot dataPath of the dataset),
    \code{geneExprFile} (slot geneExprFile of the dataset),
    \code{classesFile} (slot classesFile of the dataset),
    if the \code{topic} is missing then the whole dataset object is replaced.}
  \item{value}{The replacement value.}
}

\value{
    The methods modifies the object of class assessment and returned the slot modified
    accordingly to the request provided by \code{topic}.

    If 'topic' is missing
    \code{object of class dataset} the dataset corresponding to the assessment is replaced by 'value'.

    If 'topic' is "dataId"
    \code{object of class character} the 'dataId' of the dataset is replaced by 'value'

    If 'topic' is "dataPath"
    \code{object of class character} the 'dataPath' of the dataset is replaced by 'value'

    If 'topic' is "geneExprFile"
    \code{object of class character} the 'geneExprFile' of the dataset is replaced by 'value'

    If 'topic' is "classesFile"
    \code{object of class character} the 'classesFile' of the dataset is replaced by 'value'
}
\author{ Camille Maumet }

\seealso{
    \code{\linkS4class{assessment}},  \code{\link{getDataset-methods}}
}
\examples{

\dontrun{
aDataset <- new("dataset", dataId="vantVeer_70", dataPath="pathToFile")
aDataset <- loadData(aDataset)

expeOfInterest <- new("assessment", dataset=aDataset,
                                   noFolds1stLayer=10,
                                   noFolds2ndLayer=9,
                                   classifierName="svm",
                                   typeFoldCreation="original",
                                   svmKernel="linear",
                                   noOfRepeat=2,
                                   featureSelectionOptions=new("geneSubsets", optionValues=c(1,2,3,4,5,6)))

# Modify the dataId
getDataset(expeOfInterest, topic='dataId') <- "khan"
getDataset(expeOfInterest, 'dataId')

# Replace the dataset
getDataset(expeOfInterest) <- aDataset
getDataset(expeOfInterest, 'dataId')
}
}


\keyword{methods}

