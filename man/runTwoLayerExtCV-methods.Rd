\name{runTwoLayerExtCV-methods}
\docType{methods}
\alias{runTwoLayerExtCV-methods}
\alias{runTwoLayerExtCV}
\alias{runTwoLayerExtCV,assessment-method}
\title{runTwoLayerExtCV: Method to run an external two-layers cross-validation}
\description{
This method run an external two-layers cross-validation according to the options
stored in an object of class assessment. The concept of two-layers cross-validation has
been introduced by J.X. Zhu,G.J. McLachlan, L. Ben-Tovim Jonesa, I.A.Wood in 'On
selection biases with prediction rules formed from gene expression data' and by
I. A. Wood, P. M. Visscher, and K. L. Mengersen in 'Classification based upon gene
expression data: bias and precision of error rates' (cf. section References). This technique
of cross-validation is used to determine an unbiased estimate of the best error rate
(using the best size of subset for RFE-SVM, of the best threshold for NSC) when
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
    has been computed, therfore, the slot \code{resultRepeated2LayerCV} is no more NULL.
    This methods print out the key results of the assessment, to access the full detail
    of the results, the user must call the method \code{getResults}.
}

\seealso{
    \code{\linkS4class{assessment}}, \code{\link{getResults}}, \code{\link{runOneLayerExtCV-methods}}
}

\references{J.X. Zhu,  G.J. McLachlan, L. Ben-Tovim, I.A. Wood (2008),
"On selection biases with prediction rules formed from gene expression data",
Journal of Statistical Planning and Inference, 38:374-386.

I.A. Wood, P.M. Visscher, and K.L. Mengersen
"Classification based upon gene expression data: bias and precision of error rates"
Bioinformatics, June 1, 2007; 23(11): 1363 - 1370. }

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

myExpe <- runTwoLayerExtCV(myExpe)
}

\keyword{methods}
