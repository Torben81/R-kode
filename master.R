#packages
library(MASS)
library(gRim)
library(mlbench)
library(glasso)
library(pROC)
library(ggplot2)

#Data 
data(Satellite)
set.seed(5)
sampleSize <- sample(nrow(Satellite),4435)
trainSat <- Satellite[sampleSize,]
testSat <- Satellite[-sampleSize,]

data(breastcancer)
BC <- breastcancer
BC <- BC[sort(rownames(BC)),]
set.seed(5)
samplesize <- sample(nrow(BC), 125)
trainBC <- BC[samplesize,]
testBC <- BC[-samplesize,]


source("LDAQDA.R")
source("HLM.R")
source("RandomForest.R")