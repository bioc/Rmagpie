data('vV70genesDataset')

expeOfInterest <- new("experiment", dataset=vV70genes,
                                   noFolds1stLayer=3,
                                   noFolds2ndLayer=2,
                                   classifierName="svm",
                                   typeFoldCreation="original",
                                   svmKernel="linear",
                                   noOfRepeat=10,
                                   featureSelectionOptions=new("geneSubsets", optionValues=c(1,2,3,4,5,6)))

expeOfInterest <- runOneLayerExtCV(expeOfInterest)
plotErrorsSummaryOneLayerCV(expeOfInterest)
plotErrorsRepeatedOneLayerCV(expeOfInterest)

expeOfInterest <- runTwoLayerExtCV(expeOfInterest)
plotErrorsFoldTwoLayerCV(expeOfInterest)