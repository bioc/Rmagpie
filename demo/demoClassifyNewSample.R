data('vV70genesDataset')

expeOfInterest <- new("experiment", dataset=vV70genes,
                                   noFolds1stLayer=10,
                                   noFolds2ndLayer=9,
                                   classifierName="svm",
                                   typeFoldCreation="original",
                                   svmKernel="linear",
                                   noOfRepeat=2,
                                   featureSelectionOptions=new("geneSubsets", optionValues=c(1,2,4,8,16,32,64,70)))

# Build the final classifier
expeOfInterest <- findFinalClassifier(expeOfInterest)

classifyNewSamples(expeOfInterest, "pathToFile/testSamples_geneExpr.txt", 4)

expeOfInterest <- runOneLayerExtCV(expeOfInterest)

classifyNewSamples(expeOfInterest, "pathToFile/testSamples_geneExpr.txt")
