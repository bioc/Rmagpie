# Copyright (c) 2008 Camille Maumet
# GPL >= 3 (LICENSE.txt) licenses.

#--------------------------- loadData-methods --------------------------------
# Function and methods to load the microarray data from the files to the
# ExpressionSet (slot eset) of a dataset object.
#
# Author: Camille Maumet
# Creation: March 2008
# Last Modified: 17 Jul. 2008
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Read data from file provided by <object>, create an ExpressionSet and store
# it in the <object> (dataset), optionnal miame and annotation can be added
#
# @param    object (dataset) Dataset containing the data to be loaded
#           miame (MIAME) Optional MIAME to describe the dataset
#           annotation (character) Optional annontation to describe the dataset
#
# @return   object (dataset) dataset containing the loaded expression set
#-------------------------------------------------------------------------------
loadData <- function(object, pathGeneExpr, pathClasses, miame, annotation){
    # ---- Test the validity of the parameters ----
    if (! missing(miame)) {
        if (! is(miame, "MIAME")){
            miame <- new("MIAME")
            warning("In loadData, 'miame' is ignored because it's not an object of class MIAME")
        }
    }  else {
        miame <- new("MIAME")
    }
    if (! missing(annotation)) {
        if (length(annotation) != 1 || mode(annotation) != "character"){
            annotation <- ""
            warning("In loadData, 'annotation' is ignored because it's not of class character")
        }
    } else {
        annotation <- ""
    }

    # ---- Core of the method ----
    #pathGeneExpr <- file.path(getDataPath(object), getGeneExprFile(object))
    #pathClasses <- file.path(getDataPath(object), getClassesFile(object))

    # Copy the gene expression levels for each sample
    geneExprs <- as.matrix(read.table(pathGeneExpr, sep = "", header = TRUE, row.names = 1))

    # Copy the class output for each sample
    pData <- data.frame(t(read.table(pathClasses, row.names = 1, header = TRUE, sep = "", colClasses="factor")))
    # Metadata to describe the data
    metadata <- data.frame(labelDescription = c("Case/control status"), row.names = c("type"))

    phenoData <- new("AnnotatedDataFrame", data = pData, varMetadata = metadata)

    # Creation of the ExpressionSet containing all the data of our assessment
    eset <- new("ExpressionSet", exprs = geneExprs, phenoData = phenoData, assessmentData = miame, annotation = annotation)

    return(eset)
}