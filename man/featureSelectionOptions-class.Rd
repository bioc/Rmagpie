\name{featureSelectionOptions-class}
\docType{class}
\alias{featureSelectionOptions-class}
\alias{getOptionValues,featureSelectionOptions-method}

\title{"featureSelectionOptions": A virtual class to store the options of a feature selection}
\description{This virtual class has two descendants: geneSubsets and thresholds.
As a virtual class, you can't create an object of class featureSelectionOptions.}

\section{Slots}{
	\describe{
    \item{\code{optionValues}:}{\code{numeric} (vector). Value of the possible options}
    \item{\code{noOfOptions}:}{\code{numeric}. Total number of options}
  }
}
\section{Methods}{
  \describe{
    \item{\code{getOptionValues(featureSelectionOptions)}}{
    Retreive the value of options (slot optionValues)}
    \item{\code{getNoOfOptions(featureSelectionOptions)}}{Retreive
    the number of options (slot featureSelectionOptions)}
    }
}

\author{Camille Maumet}

\seealso{
	\code{\linkS4class{geneSubsets}}, \code{\linkS4class{thresholds}}
}
\keyword{classes}
