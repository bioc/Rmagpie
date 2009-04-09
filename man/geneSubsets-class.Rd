\name{geneSubsets-class}
\docType{class}
\alias{geneSubsets-class}
\alias{getMaxSubsetSize,geneSubsets-method}
\alias{getMaxSubsetSize<-,geneSubsets-method}
\alias{getMaxSubsetSize}
\alias{getMaxSubsetSize<-}
\alias{getNoModels,geneSubsets-method}
\alias{getNoModels}
\alias{getSpeed,geneSubsets-method}
\alias{getSpeed<-,geneSubsets-method}
\alias{getSpeed}
\alias{getSpeed<-}
\alias{getSubsetsSizes,geneSubsets-method}
\alias{getSubsetsSizes<-,geneSubsets-method}
\alias{getSubsetsSizes}
\alias{getSubsetsSizes<-}

\title{geneSubsets: A class to handle the sizes of gene susbets to be tested during forward gene selection}
\description{Forward gene selection is usually a computationally expensive
task. To reduce the computation expense one may want to do not consider one gene
at a time but chunks of genes. This class store the sizes of gene susbets to be
tested during forward gene selection.}
\section{Creating objects}{
    \code{new("geneSubsets", optionValues)}

    Create a geneSubsets, the sizes of the different subsets are determined by
    \code{optionValues}. The size of the biggest subset \code{maxSubsetSize} and
    the number of subsets to be tried \code{noOfOptions} are automatically deducted.
    The speed is set to \code{high} is there are less models than the size of the
    biggest subset and 'slow' if not.

    \code{new("geneSubsets", maxSubsetSize, speed="high")}

    Create a geneSubsets, with a biggest subset of size \code{maxSubsetSize}. If
    the \code{speed} is \code{high} the sizes of the subsets are increased by a power of 2
    from 1 to the biggest power of 2 smaller than \code{maxSubsetSize}. If the
    speed is \code{slow} the sizes of the subsets are increased by 1
    from 1 to the \code{maxSubsetSize}.

}
\section{Slots}{
	 \describe{
    \item{\code{maxSubsetSize}:}{\code{numeric}. Size of the biggest subset}
    \item{\code{optionValues}:}{\code{numeric} (vector). Sizes of the subsets in acsending order}
    \item{\code{noOfOptions}:}{\code{numeric}. Total number of subsets
    to be tried during backward gene selection}
    \item{\code{speed}:}{\code{character}. Speed of the backward feature selection.
    \code{high} if the number of models is smaller than the size of the biggest subset,
    \code{slow} if not.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{getMaxSubsetSize(geneSubsets), getMaxSubsetSize(geneSubsets)<-}}{
    Retreive and modify the size of the biggest subset (slot maxSubsetSize)}
    \item{\code{getOptionValues(geneSubsets), getOptionValues(geneSubsets)<-}}{Retreive
    and modify the sizes of the subsets of features (slot optionValues)}
    \item{\code{getNoOfOptions(geneSubsets)}}{Retreive the total number of subsets
    to be tried during backward gene selection (slot noModels)}
    \item{\code{getSpeed(geneSubsets), getSpeed(geneSubsets)<-}}{Retreive and modify the speed of
    the backward feature selection. (slot speed)}
    }
}
\author{ Camille Maumet }


\seealso{
	\code{\linkS4class{thresholds}},\code{\linkS4class{assessment}}
}
\examples{
geneSubset235 <- new("geneSubsets", optionValues=c(2,3,5))
geneSubset235
getSubsetsSizes(geneSubset235)
getSpeed(geneSubset235)
getMaxSubsetSize(geneSubset235)

geneSubsetMax60 <- new("geneSubsets", maxSubsetSize=60, speed="slow")
geneSubsetMax60

geneSubsetSlow <- new("geneSubsets", maxSubsetSize=70, speed="slow")
geneSubsetSlow

getMaxSubsetSize(geneSubsetMax60) <- 70
geneSubsetMax60

newSizes <- c(1,2,3,4,5)
getSubsetsSizes(geneSubsetMax60) <- newSizes
geneSubsetMax60

getSpeed(geneSubset235) <- 'slow'
geneSubset235
}
\keyword{classes}
