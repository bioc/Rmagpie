\name{thresholds-class}
\docType{class}
\alias{thresholds-class}
\alias{getNoThresholds,thresholds-method}
\alias{getNoThresholds<-,thresholds-method}
\alias{getNoThresholds}
\alias{getNoThresholds<-}
\alias{getOptionValues,thresholds-method}
\alias{getOptionValues<-,thresholds-method}
\alias{getOptionValues}
\alias{getOptionValues<-}

\title{thresholds: A class to handle the thresholds to be tested during training of the Nearest
    Shrunken Centroid}
\description{The Nearest Shrunken Centroid is computed using a threshold. This
threshold is usually determined by finding the best threshold value over a set of
values by finding the threshold leading to the best error rate assessed by cross-validation.
This class stores the values of thresholds to be tried. If the user wants to use default
values it's also possible.}

\section{Creating objects}{
    \code{new("thresholds")}

    Create an empty thresholds. The default thresholds values will be computed
    and this object updated as soon as it is linked in an assessment.

    \code{new("thresholds", optionValues)}

    Create a thresholds, containing the thresholds values defined by
    \code{optionValues}. The slot \code{noOfOptions} is automatically updated.

}

\section{Slots}{
	 \describe{
    \item{\code{optionValues}:}{\code{numeric} Values of the thresholds, if
    \code{optionValues} has length zero then the default thresholds values must be used. }
    \item{\code{noOfOptions}:}{\code{numeric} Number of thresholds. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{featureSelectionOptions}"}, directly.
}
\section{Methods}{
  \describe{
    \item{\code{getNoThresholds(thresholds)}}{Retreive
    the number of the thresholds (slot noOfOptions)}
    \item{\code{getOptionValues(thresholds)},  \code{getOptionValues(thresholds)<-}}{Retreive
    and modify the values of the thresholds (slot optionValues)}
	 }
}

\author{Camille Maumet}

\seealso{
    \code{\linkS4class{geneSubsets}}, \code{\linkS4class{assessment}}
}
\examples{
# Empty thresholds, the default values will be used when added to an assessment
emptThresholds <- new("thresholds")
getOptionValues(emptThresholds)
getNoThresholds(emptThresholds)

# Another thresholds
thresholds <- new("thresholds", optionValues=c(0,0.1,0.2,1,2))
getOptionValues(thresholds)
getNoThresholds(thresholds)

# Set the thresholds
newThresholds <- c(0.1,0.2,0.5,0.6,1)
getOptionValues(thresholds) <- newThresholds
getOptionValues(thresholds)
getNoThresholds(thresholds)
}
\keyword{classes}
