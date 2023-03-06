#Alex Stocking and Sivakami Thinnappan
install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("edgeR")
library(limma)
library(edgeR)
#load in data
expresData <- read.delim("dgrp_expression_Female.txt", row.names = 1)
eyeSizeData <- read.delim("rpr.txt")
#remove unnecessary characters and split groups
colnames(expresData) <- gsub("line_", "", colnames(expresData))
group1names <- seq(1, length(expresData), 2)
splitData <- expresData[, group1names]
colnames(splitData) <- gsub(".1$", "", colnames(splitData))
#check what strains are in the data and create a finalized data frame for both
finalNames <- NULL
for (i in 1:nrow(eyeSizeData)) {
  target <- as.character(eyeSizeData[i,1])
  if(target %in% colnames(splitData)){
    finalNames <- append(finalNames, target)
  }
}
finalExpressData <- splitData[, finalNames]
droppedValues <- NULL
for (i in 1:nrow(eyeSizeData)) {
  target <- as.character(eyeSizeData[i, 1])
  if((target %in% finalNames)==FALSE){
    droppedValues <- append(droppedValues, i)
  }
}
finalSizeData <- eyeSizeData[-droppedValues,]
#visualize
hist(finalSizeData$Mean.Eye.Size, breaks = 20, xlab = "Mean Eye Size", main = "Histogram of Mean Eye Sizes")
expressMeans <- rowMeans(finalExpressData)
hist(expressMeans, breaks = 20, xlab = "Expression Value", main = "Histogram of Mean Expression Values")
#calculate quantiles
quantiles <- quantile(finalSizeData$Mean.Eye.Size, probs = c(0, 0.2, 0.8, 1))
groups <- c("small", "medium", "large")
eyeSizes <- cut(finalSizeData$Mean.Eye.Size, breaks = quantiles, include.lowest = TRUE, labels = groups)
#create a DGElist and recalculate normalization factors
expresList <- DGEList(as.matrix(finalExpressData))
expresList <- calcNormFactors(expresList)
#Filter out low expression genes
droppedGenes <- which(apply(cpm(expresList), 1, max) < 30)
#no dropped genes
#create a multidimensional scaling plot using the col names as labels
plotMDS(expresList, col = as.numeric(eyeSizes))
#create a model matrix
modelMat <- model.matrix(~0 + eyeSizes)
#use voom and apply it to a table to be read
vObj <- voom(expresList, modelMat, plot = T)
vFit <- lmFit(vObj, design = modelMat)
vFitBayes <- eBayes(vFit)
fitTable <- topTable(vFitBayes)
#use contrasts to get better results
contrasts <- makeContrasts(eyeSizeslarge - eyeSizesmedium, eyeSizesmedium - eyeSizessmall, eyeSizeslarge - eyeSizessmall, levels = modelMat)
contrFit <- contrasts.fit(vFitBayes, contrasts)
contrFitBayes <- eBayes(contrFit)
contrFitTable <- topTable(contrFitBayes, resort.by = "p")
#print out the table data
fitTable
contrFitTable
#transpose the data so ggscatter stops complaining
library(ggpubr)
library(data.table)
transData <- rbind(finalExpressData, as.vector(finalSizeData[, 2]))
transCol <- colnames(transData)
transRow <- rownames(transData)
transData <- transpose(transData)
rownames(transData) <- transCol
colnames(transData) <- transRow
colnames(transData)[colnames(transData) == "18141"] <- "eyeSize"
#do the plots for the top 5 genes found in the normal fitted table
ggscatter(data = transData, x = "eyeSize", y = rownames(fitTable[1,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(fitTable[1,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(fitTable[2,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(fitTable[2,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(fitTable[3,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(fitTable[3,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(fitTable[4,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(fitTable[4,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(fitTable[5,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(fitTable[5,]))
#do the plots for the top 5 genes in the contrast fitted table
ggscatter(data = transData, x = "eyeSize", y = rownames(contrFitTable[1,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(contrFitTable[1,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(contrFitTable[2,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(contrFitTable[2,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(contrFitTable[3,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(contrFitTable[3,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(contrFitTable[4,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(contrFitTable[4,]))
ggscatter(data = transData, x = "eyeSize", y = rownames(contrFitTable[5,]), add = "reg.line", xlab = "Mean Eye Size", ylab = "Expression Value", title = rownames(contrFitTable[5,]))

