# Packages
library(MASS)
library(gRim)
library(mlbench)
library(glasso)
library(pROC)
library(ggplot2)

# Loading data 
data(Satellite)
set.seed(5)
sampleSizeSat <- sample(nrow(Satellite),4435)
trainSat <- Satellite[sampleSizeSat,]
testSat <- Satellite[-sampleSizeSat,]

data(breastcancer)
BC <- breastcancer
BC <- BC[sort(rownames(BC)),]
set.seed(5)
sampleSizeBC <- sample(nrow(BC), 125)
trainBC <- BC[sampleSizeBC,]
testBC <- BC[-sampleSizeBC,]

# Set working directory
setwd("C:/Users/Torben/Dropbox/speciale/R-kode/git/Speciale")

# Loading R-scripts. Takes about 30 minutes.
Sys.time() -> start
source("LDAQDASatellite.R")
source("HLM.R")
source("HLMSatellite.R")
source("HLMBCancer.R")
source("glassoBCancer.R")
source("ROCCurves.R")
Sys.time()-start
