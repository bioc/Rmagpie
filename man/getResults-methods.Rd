\name{getResults-methods}
\docType{methods}
\alias{getResults-methods}
\alias{getResults}
\alias{getResults,assessment-method}
\title{getResults Method to access the result of one-layer and two-layers cross-validation from an assessment}
\description{
This method provides an easy interface to access the results of one-layer
and two-layers of cross-validation directly from an object assessment.
}

\section{Methods}{
\describe{
\item{object = "assessment"}{The method is only applicable on objects of class
        assessment.}
}}

\arguments{
  \item{object}{\code{Object of class assessment}. Object assessment of interest}
  \item{layer}{\code{numeric}. Indice that states which layer of cross-validation must be
    accessed.
    Set to \code{1} to acces the one-layer cross-validation,
    Set to \code{c(1,i)} to acces the ith repeat of the one-layer cross-validation,
    Set to \code{2} to acces to the two-layers cross-validation,
    Set to \code{c(2,i)} to access the ith repeat of the two-layers cross-validation,
    Set to \code{c(2,i,j)} to access the jth inner layer of ith repeat of the two-layers cross-validation,
    Set to \code{c(2,i,j,k)} to access the kth repeat of the jth inner layer of ith repeat
    of the two-layers cross-validation}
  \item{topic}{character. Argument that specifies which kind of
    result is requested, the possible values are
    \code{"errorRate"}: Access to cross-validation error rate, standard error on cross-validated
        error rate, error rate per fold, number of samples per fold and error rate per class,
    \code{"selectedGenes"}: Access to the genes selected for each fold or their frequency of
        selection among the folds and the repeats,
    \code{"bestOptionValue"}: For one-layer of cross-validation, access to the best
    option value (size of gene subset for SVM-RFE or thresholds for NSC) corresponding
    to the best value of the cross-validated error rate. For the two-layers of
    cross-validation, access the average best option value (over the repeats and folds).
    \code{"executionTime"}: Time used to run the selected layer in seconds.}
  \item{errorType}{character. Optional, ignored if topic is not \code{"errorRate"}.
    Specify the type of error rate requested, the possible values are:
    \code{missing} or \code{"all"} to access all the following error rates
    \code{"cv"} to access the cross-validated error rate,
    \code{"se"} to access the standard error on the cross validated error rate,
    \code{"fold"} to access the error rate per fold (not available in certain cases see
    section value for more details),
    \code{"noSamplesPerFold"} to access the number of samples in each fols (not available in
    certain cases see section value for more details),
    \code{"class"} to acces the error rate per class}
\item{genesType}{character. Optional, ignored if topic is not \code{"selectedGenes"}.
    Specify the type of display of genes selected, the possible values are:
    missing
    \code{"fold"} to access the genes selected for each fold (not available in certain case see
    section value for more details),
    \code{"frequ"} to access the genes order by their frequency among the folds(not available in certain case see
    section value for more details)}

}

\value{
    if there is no error, the value returned by the method depends on the arguments namely,
    \code{layer}, \code{topic}, \code{errorType} and \code{genesType}.

    If \code{layer} is 1
        \item{General}{Get the results of the repeated one-layer cross-validation corresponding to
        the \code{object} of class assessment. If the one-layer cross-validation has not been performed and the user try to access it
        then the function return an error indicating that he must call \code{runOneLayerExtCV} first.}
        \item{if topic is \code{"errorRate"}}{}
            \item{If errorType=\code{"all"} or is \code{missing}}{All the following error rates}
            \item{If errorType=\code{"cv"}}{\code{numeric}. Cross-validated error-rate for each value of option tried
            obtained by one-layer of cross-validation (1 value per value of option).}
            \item{If errorType=\code{"se"}}{\code{numeric}. Standard error on cross-validated error-rate for each
            value of option tried obtained by one-layer of cross-validation (1 value
            per value of option).}
            \item{If errorType=\code{"class"}}{numeric. Class cross-validated error rate error for each value of option
             tried obtained by one-layer of cross-validation (1 value per class and value of option).}
            \item{Else}{Error signaling that the topic is not appropriate.}
        \item{if topic is \code{"genesSelected"}}{}
            \item{If genesType=\code{"freq"} or is missing}{
            \code{list}. Each elelement of the list corresponds to the genes selected for each model
            ordered by frequency.}
            \item{Else}{
            Error signaling that the topic is not appropriate.}
        \item{if topic is \code{"bestOptionValue"}}{
        Size of subset (for RFE-SVM) or threshold (for NSC) corresponding to the minimum cross-validated error rate.}
        \item{if topic is \code{"executionTime"}}{
        Time in second to perform this one-layer cross-validation.}
    If \code{layer} is c(1,i)
        \item{General}{Get the results of the ith repeat of the one-layer cross-validation corresponding to
        the \code{object} of class assessment. If the one-layer cross-validation has not been performed and the user try to access it
        then the function return an error indicating that he must call \code{runOneLayerExtCV} first.}
        \item{if topic is \code{"errorRate"}}{   }
            \item{If errorType=\code{"all"} or is \code{missing}}{
            All the following error rates}
            \item{If errorType=\code{"cv"}}{
            numeric. Cross-validated error-rate for each value of option tried
            obtained by one-layer of cross-validation on the ith repeat(1 value
            per subset).}
            \item{If errorType=\code{"se"}}{
            numeric. Standard error on cross-validated error-rate for each value of option
             tried obtained by one-layer of cross-validation on the ith repeat
            (1 value per value of option).}
            \item{If errorType=\code{"class"}}{
            numeric. Class cross-validated error rate error for each value of option tried
            obtained by one-layer of cross-validation on the ith repeat
            (1 value per class and value of option).}
            \item{If errorType=\code{"fold"}}{
            numeric. Class cross-validated error rate error for each fold and each value of option
            tried obtained by one-layer of cross-validation on the ith repeat
            (1 value per class and value of option).}
            \item{Else}{
            Error signaling that the topic is not appropriate.}
        \item{if topic is \code{"genesSelected"}}{}
            \item{If genesType=\code{"freq"} or is missing}{
            list. Each elelement of the list corresponds to the genes selected for each model
            ordered by frequency.}
            \item{If genesType=\code{"fold"}}{
            list. Each elelement of the list corresponds to a model and contains a list
            of which one element correspond to the genes selected in a particular fold.}
            \item{Else}{
            Error signaling that the topic is not appropriate.}
        \item{if topic is \code{"bestOptionValue"}}{
        numeric. Size of subset (for RFE) or threshold (for NSC) corresponding to the minimum
        cross-validated error rate in the ith repeat of the one-layer cross-validation.}
        \item{if topic is \code{"executionTime"}}{
        Time in second to perform this repeat of one-layer cross-validation.}

    If \code{layer} is 2
        \item{General}{Get the results of the repeated two-layers cross-validation corresponding to
        the \code{object} of class assessment. If the two-layer cross-validation has not been performed and the user try to access it
        then the function return an error indicating that he must call \code{runTwoLayerExtCV} first.}
        \item{if topic is 'errorRate'}{}
            \item{If errorType=\code{"all"} or is \code{missing}}{All the following error rates}
            \item{If errorType=\code{"cv"}}{
            numeric. Cross-validated error-rate obtained by two-layers of
            cross-validation (1 value).}
            \item{If errorType=\code{"se"}}{
            numeric. Standard error on cross-validated error-rate obtained by
            two-layers of cross-validation (1 value).}
            \item{If errorType=\code{"class"}}{
            numeric. Class cross-validated error rate obtained by two-layers (1 value per class)}
            \item{Else}{
            Error signaling that the topic is not appropriate.}
        \item{if topic is \code{"bestOptionValue"}}{
        numeric. Average best number of genes for SVM-RFE of threshold for NSc obtained among the folds.}
        \item{if topic is \code{"executionTime"}}{
        Time in second to perform this two-layers cross-validation.}
    
    If \code{layer} is c(2,i)
        \item{General}{Get the results of the ith repeated of the two-layers cross-validation corresponding to
        the \code{object} of class assessment. If the two-layer cross-validation has not been performed and the user try to access it
        then the function return an error indicating that he must call \code{runTwoLayerExtCV} first.}
        \item{if topic is 'errorRate'}{   }
            \item{If errorType=\code{"all"} or is \code{missing}}{All the following error rates}
            \item{If errorType=\code{"cv"}}{
            numeric. Cross-validated error-rate obtained by two-layers of
            cross-validation in this repeat. (1 value).}
            \item{If errorType=\code{"se"}}{
            numeric. Standard error on cross-validated error-rate obtained by
            two-layers of cross-validation in this repeat (1 value).}
            \item{If errorType=\code{"class"}}{
            numeric. Class cross-validated error rate obtained by two-layers in this repeat}
            \item{If errorType=\code{"fold"}}{
            numeric. Error rate obtained on each of the folds in the second
            layer in this repeat(1 value per fold).
            of cross-validation (value per class).}
            \item{Else}{
            Error signaling that the topic is not appropriate.}
        \item{if topic is \code{"genesSelected"}}{  }
            \item{If genesType=\code{"fold"} or is missing}{
            list. Each elelement of the list corresponds to a fold and contains a list
            of the genes selected in this particular fold.}
            \item{Else}{
            Error signaling that the topic is not appropriate.}
        \item{if topic is \code{"bestOptionValue"}}{
        numeric. Average best number of genes obtained among the folds in this repeat.}
        \item{if topic is \code{"executionTime"}}{
        Time in second to perform this repeat of two-layers cross-validation.}
    \item{If \code{layer} is c(2,i,j)}{This layer corresponds to the jth inner layer of one-layer
    cross-validation performed inside the ith repeat of the two-layers cross-validation. The returned
    values are similar to the one returned by a repeated one-layer cross-validation.}
    \item{If \code{layer} is c(2,i,j,k)}{This layer corresponds to the kth repeat of the jth inner layer
    of one-layer cross-validation performed inside the ith repeat. The returned values are
    similar to the one returned by a repeat of one-layer cross-validation.}

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

mySubsets <- new("geneSubsets", optionValues=c(1,2,4,8,16,32,64,70))
myassessment <- new("assessment", dataset=vV70genes,
                                   noFolds1stLayer=5,
                                   noFolds2ndLayer=4,
                                   classifierName="svm",
                                   typeFoldCreation="original",
                                   svmKernel="linear",
                                   noOfRepeat=2,
                                   featureSelectionOptions=mySubsets)

myassessment <- runOneLayerExtCV(myassessment)
myassessment <- runTwoLayerExtCV(myassessment)

# --- Access to one-layer CV ---
# errorRate
# 1-layer CV: error Rates
getResults(myassessment, 1, 'errorRate')
# 1-layer CV: error Rates - all")
getResults(myassessment, 1, 'errorRate', errorType='all')
# 1-layer CV: error Rates - cv
getResults(myassessment, 1, 'errorRate', errorType='cv')
# 1-layer CV: error Rates - se
getResults(myassessment, 1, 'errorRate', errorType='se')
# 1-layer CV: error Rates - class
getResults(myassessment, 1, 'errorRate', errorType='class')

# genesSelected
# 1-layer CV: genes Selected
getResults(myassessment, 1, 'genesSelected')
# 1-layer CV: genes Selected - frequ
getResults(myassessment, 1, 'genesSelected', genesType='frequ')
# 1-layer CV: genes Selected - model 7
getResults(myassessment, 1, 'genesSelected', genesType='frequ')[[7]]
getResults(myassessment, 1, 'genesSelected')[[7]]

# bestOptionValue
# 1-layer CV: best number of genes
getResults(myassessment, 1, 'bestOptionValue')

# executionTime
# 1-layer CV: execution time
getResults(myassessment, 1, 'executionTime')

# --- Access to 2nd repeat of one-layer CV ---
# Error rates
# 1-layer CV repeat 2: error Rates
getResults(myassessment, c(1,2), 'errorRate')
# 1-layer CV repeat 2: error Rates - all
getResults(myassessment, c(1,2), 'errorRate', errorType='all')
# 1-layer CV repeat 2: error Rates - cv
getResults(myassessment, c(1,2), 'errorRate', errorType='cv')
# 1-layer CV repeat 2: error Rates - se
getResults(myassessment, c(1,2), 'errorRate', errorType='se')
# 1-layer CV repeat 2: error Rates - fold
getResults(myassessment, c(1,2), 'errorRate', errorType='fold')
# 1-layer CV repeat 2: error Rates - noSamplesPerFold
getResults(myassessment, c(1,2), 'errorRate', errorType='noSamplesPerFold')
# 1-layer CV repeat 2: error Rates - class
getResults(myassessment, c(1,2), 'errorRate', errorType='class')

# genesSelected
# 1-layer CV repeat 2: genes Selected
getResults(myassessment, c(1,2), 'genesSelected')
# 1-layer CV repeat 2: genes Selected - frequ
getResults(myassessment, c(1,2), 'genesSelected', genesType='frequ')
# 1-layer CV repeat 2: genes Selected - model 7 (twice)
getResults(myassessment, c(1,2), 'genesSelected', genesType='frequ')[[7]]
getResults(myassessment, c(1,2), 'genesSelected')[[7]]
# 1-layer CV repeat 2: genes Selected - fold
getResults(myassessment, c(1,2), 'genesSelected', genesType='fold')

# 1-layer CV repeat 2: best number of genes
getResults(myassessment, c(1,2), 'bestOptionValue')

# 1-layer CV repeat 2: execution time
getResults(myassessment, c(1,2), 'executionTime')

# --- Access to two-layers CV ---
# Error rates
# 2-layer CV: error Rates
getResults(myassessment, 2, 'errorRate')
# 2-layer CV: error Rates - all
getResults(myassessment, 2, 'errorRate', errorType='all')
# 2-layer CV: error Rates - cv
getResults(myassessment, 2, 'errorRate', errorType='cv')
# 2-layer CV: error Rates - se
getResults(myassessment, 2, 'errorRate', errorType='se')
# 2-layer CV: error Rates - class
getResults(myassessment, 2, 'errorRate', errorType='class')

# bestOptionValue
# 2-layer CV: best number of genes (avg)
getResults(myassessment, 2, 'bestOptionValue')

# executionTime
# 2-layer CV: execution time
getResults(myassessment, 2, 'executionTime')

# --- Access to two-layers CV access to repeats ---
# Error rates
# 2-layer CV repeat 1: error Rates
getResults(myassessment, c(2,1), 'errorRate')
# 2-layer CV repeat 1: error Rates - all
getResults(myassessment, c(2,1), 'errorRate', errorType='all')
# 2-layer CV repeat 1: error Rates - cv
getResults(myassessment, c(2,1), 'errorRate', errorType='cv')
# 2-layer CV repeat 1: error Rates - se
getResults(myassessment, c(2,1), 'errorRate', errorType='se')
# 2-layer CV repeat 1: error Rates - fold
getResults(myassessment, c(2,1), 'errorRate', errorType='fold')
# 2-layer CV repeat 1: error Rates - noSamplesPerFold
getResults(myassessment, c(2,1), 'errorRate', errorType='noSamplesPerFold')
# 2-layer CV repeat 1: error Rates - class
getResults(myassessment, c(2,1), 'errorRate', errorType='class')

# genesSelected
# 2-layer CV repeat 1: genes Selected
getResults(myassessment, c(2,1), 'genesSelected')
# 2-layer CV repeat 1: genes Selected - fold
getResults(myassessment, c(2,1), 'genesSelected', genesType='fold')

# 2-layer CV repeat 1: best number of genes
getResults(myassessment, c(2,1), 'bestOptionValue')

# 2-layer CV repeat 1: execution time
getResults(myassessment, c(2,1), 'executionTime')

# --- Access to one-layer CV inside two-layers CV ---
# errorRate
# 2-layer CV repeat 1 inner layer 3: error Rates
getResults(myassessment, c(2,1,3), 'errorRate')
# 2-layer CV repeat 1 inner layer 3: error Rates - all
getResults(myassessment, c(2,1,3), 'errorRate', errorType='all')
# 2-layer CV repeat 1 inner layer 3: error Rates - cv
getResults(myassessment, c(2,1,3), 'errorRate', errorType='cv')
# 2-layer CV repeat 1 inner layer 3: error Rates - se
getResults(myassessment, c(2,1,3), 'errorRate', errorType='se')
# 2-layer CV repeat 1 inner layer 3: error Rates - class
getResults(myassessment, c(2,1,3), 'errorRate', errorType='class')

# genesSelected
# 2-layer CV repeat 1 inner layer 3: genes Selected
getResults(myassessment, c(2,1,3), 'genesSelected')
# 2-layer CV repeat 1 inner layer 3: genes Selected - frequ
getResults(myassessment, c(2,1,3), 'genesSelected', genesType='frequ')
# 2-layer CV repeat 1 inner layer 3: genes Selected - model 7
getResults(myassessment, c(2,1,3), 'genesSelected', genesType='frequ')[[7]]
getResults(myassessment, c(2,1,3), 'genesSelected')[[7]]

# bestOptionValue
# 2-layer CV repeat 1 inner layer 3: best number of genes
getResults(myassessment, c(2,1,3), 'bestOptionValue')

# executionTime
# 2-layer CV repeat 1 inner layer 3: execution time
getResults(myassessment, c(2,1,3), 'executionTime')

 # --- two-layers CV access to repeat 1, inner layer 2 repeat 2 ---
# Error rates
# 2-layer CV inner layer 3 repeat 2: error Rates
getResults(myassessment, c(2,1,3,1), 'errorRate')
# 2-layer CV repeat 1 inner layer 3 repeat 1: error Rates - all
getResults(myassessment, c(2,1,3,1), 'errorRate', errorType='all')
# 2-layer CV repeat 1 inner layer 3 repeat 1: error Rates - cv
getResults(myassessment, c(2,1,3,1), 'errorRate', errorType='cv')
# 2-layer CV repeat 1 inner layer 3 repeat 1: error Rates - se
getResults(myassessment, c(2,1,3,1), 'errorRate', errorType='se')
# 2-layer CV repeat 1 inner layer 3 repeat 1: error Rates - class
getResults(myassessment, c(2,1,3,1), 'errorRate', errorType='class')
# 2-layer CV repeat 1 inner layer 3 repeat 1: error Rates - fold
getResults(myassessment, c(2,1,3,1), 'errorRate', errorType='fold')
# 2-layer CV repeat 1 inner layer 3 repeat 1: error Rates - noSamplesPerFold
getResults(myassessment, c(2,1,3,1), 'errorRate', errorType='noSamplesPerFold')

# genesSelected
# 2-layer CV repeat 1 inner layer 3 repeat 1: genes Selected
getResults(myassessment, c(2,1,3,1), 'genesSelected')
# 2-layer CV repeat 1 inner layer 3 repeat 1: genes Selected - fold
getResults(myassessment, c(2,1,3,1), 'genesSelected', genesType='fold')
# 2-layer CV repeat 1 inner layer 3 repeat 1: genes Selected - model 3 fold 1(twice)
getResults(myassessment, c(2,1,3,1), 'genesSelected', genesType='fold')[[3]][[1]]
# 2-layer CV repeat 1 inner layer 3 repeat 1: genes Selected frequ - model 3
getResults(myassessment, c(2,1,3,1), 'genesSelected')[[3]]

# 2-layer CV repeat 1 inner layer 3 repeat 1: best number of genes
getResults(myassessment,  c(2,1,3,1), 'bestOptionValue')

# 2-layer CV repeat 1 inner layer 3 repeat 1: execution time
getResults(myassessment,  c(2,1,3,1), 'executionTime')
}


\keyword{methods}
