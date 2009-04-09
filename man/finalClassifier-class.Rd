\name{finalClassifier-class}
\docType{class}
\alias{finalClassifier-class}
\alias{getGenesFromBestToWorst,finalClassifier-method}
\alias{getGenesFromBestToWorst}
\alias{getModels,finalClassifier-method}
\alias{getModels}


\title{finalClassifier: A class to store the final classifier corresponding to an assessment}
\description{This class stores the properties of the final classifiers associated to
a given assessment. A classifier is usually available for each option value defined in
the slot featureSelectionOptions. This final classifier is obtained by running the feature
selction method on the whole dataset to find the relevant genes and then train the
classifier on the whole data considering only the relevant genes.}
\section{Creating objects}{
To generate the final classifier, call the method 'findFinalClassifier'
on an object of class assessment (\code{\link{findFinalClassifier-methods}}).
}
\section{Slots}{
	 \describe{
    \item{\code{genesFromBestToWorst}:}{\code{character}. If the feature selection
    method is RFE: the genes ordered by the weights obtained with the smallest
    subset size during RFE. If the method of featuure selection is the Nearest
    Shrunken Centroid, this slot is empty.}
    \item{\code{models}:}{\code{list of object of class svm}.If the feature selection
    method is RFE: svm models trained on the whole dataset for each size of subset
    (2 attributes: 'model', the classifier model and
    'modelFeatures' the features selected for each subset).
    If the feature selection method is NSC: the object created by pamr.train
    on the whole dataset.}
    }
}
\section{Methods}{
  \describe{
    \item{\code{getGenesFromBestToWorst(finalClassifier)}}{
    Retreive the genes ordered by their weights obtained with the smallest subset during RFE (slot genesFromBestToWorst)}
    \item{\code{getModels(finalClassifier)}}{Retreive
    the svm models for each size of subset (slot models)}
    }
}

\author{ Camille Maumet }


\seealso{
	\code{\linkS4class{finalClassifier}},\code{\linkS4class{assessment}}, \code{\link{getFinalClassifier-methods}}
}
\examples{
#dataPath <- file.path("C:", "Documents and Settings", "c.maumet", "My Documents", "Programmation", "Sources", "SVN", "R package", "data")
#aDataset <- new("dataset", dataId="vantVeer_70", dataPath=dataPath)
#aDataset <- loadData(aDataset)
data('vV70genesDataset')

mySubsets <- new("geneSubsets", optionValues=c(1,2,4,8,16,32,64,70))
expeOfInterest <- new("assessment", dataset=vV70genes,
                               noFolds1stLayer=10,
                               noFolds2ndLayer=9,
                               classifierName="svm",
                               typeFoldCreation="original",
                               svmKernel="linear",
                               noOfRepeat=2,
                               featureSelectionOptions=mySubsets)

expeOfInterest <- findFinalClassifier(expeOfInterest)

# Return the whole object of class finalClassifier
finalClassifier <- getFinalClassifier(expeOfInterest)

# Svm model corresponding to a subset of size 4 (3rd size of subset)
getModels(finalClassifier)[[3]]$model
# Relevant genes for a subset of size 4 (3rd size of subset)
getModels(finalClassifier)[[3]]$modelFeatures

# Genes ordered according to their weight after performing the RFE up to 1 gene
getGenesFromBestToWorst(finalClassifier)

}
\keyword{classes}
