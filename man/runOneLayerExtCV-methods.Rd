\name{runOneLayerExtCV-methods}
\docType{methods}
\alias{runOneLayerExtCV}
\alias{runOneLayerExtCV-methods}
\alias{runOneLayerExtCV,assessment-method}
\title{runOneLayerExtCV: Method to run an external one-layer cross-validation}
\description{
This method run an external one-layer cross-validation according to the options
stored in an object of class assessment. The concept of external cross-validation has
been introduced by G.J. McLachlan and C. Ambroise in 'Selection bias in gene extraction
on the basis of microarray gene-expression data' (cf. section References). This technique
of cross-validation is used to determine an unbiased estimate of the error rate when
feature selection is involved.
}
\section{Methods}{
\describe{
\item{object = "assessment"}{This method is only applicable on objects of class
        assessment.}
}}

\arguments{
  \item{object}{\code{Object of class assessment}. Object assessment of interest}
}

\value{
    \code{object of class assessment} in which the one-layer external cross-validation
    has been computed, therfore, the slot \code{resultRepeated1LayerCV} is no more NULL.
    This methods print out the key results of the assessment, to access the full detail
    of the results, the user must call the method \code{getResults}.
}

\seealso{
    \code{\linkS4class{assessment}}, \code{\link{getResults}}, \code{\link{runTwoLayerExtCV-methods}}
}

\references{C. Amboise and G.J. McLachlan 2002. selection bias in gene
extraction on the basis of microarray gene-expression data. PNAS,
99(10):6562-6566}

\examples{
data('vV70genesDataset')

# assessment with RFE and SVM
myExpe <- new("assessment", dataset=vV70genes,
                   noFolds1stLayer=9,
                   noFolds2ndLayer=10,
                   classifierName="svm",
                   typeFoldCreation="original",
                   svmKernel="linear",
                   noOfRepeat=2,
                   featureSelectionOptions=new("geneSubsets", optionValues=c(1,2,3,4,5,6)))

myExpe <- runOneLayerExtCV(myExpe)
}

\keyword{methods}
