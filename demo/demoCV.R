data('vV70genesDataset')

# Experiment with RFE and SVM
myExpe <- new("experiment", dataset=vV70genes,
                   noFolds1stLayer=9,
                   noFolds2ndLayer=10,
                   classifierName="svm",
                   typeFoldCreation="original",
                   svmKernel="linear",
                   noOfRepeat=2,
                   featureSelectionOptions=new("geneSubsets", optionValues=c(1,2,3,4,5,6)))

myExpe <- runOneLayerExtCV(myExpe)
myExpe <- runTwoLayerExtCV(myExpe)