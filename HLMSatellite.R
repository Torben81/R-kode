# In the follwing we make an example of model selection with the HLM compared with gRim.
d <- trainSat[trainSat[,37]=="red soil",-37]
m0 <- init(d)
set.seed(5)
system.time(
  m <- step(m0)
)
m$it
# 374 iterations in 0.28 seconds
fitm <- cmod(edgeL(m,d),d)
system.time(
  forwm <- forward(fitm,criterion="test",alpha=.01,type="unrestricted",search="headlong", steps=1)
)

# In the following we make Classification on the Satellite data using the HLM.
typesSat <- levels(Satellite$classes) # Types of classes.

# To obtain some sort of Random Forrest we apply the HLM R times on the training data.
R <- 100

fitmSat <- list()
fitModelsSat <- list()
muTrainSat <- list()
for(t in typesSat){  # Performing model selection R times for each class using the HLM.
  d <- trainSat[trainSat$classes==t, -37]
  m0 <- init(d)
  for(i in 1:R){
    mSat <- step(m0)
    mSatEdgeL <- edgeL(mSat, d)  # Extracting the edges of a model.
    fitmSat[[i]] <- cmod(mSatEdgeL,d)  # Fitting the model.
    cat(i,"\n")
  }
  fitModelsSat[[t]] <- fitmSat
  mu <- as.matrix(apply(d,2,mean))  # The empirical mean of the variabels for class t.
  muTrainSat[[t]] <- mu
}

# Estimating the conditional density of every observation in the test set given each class.
densModelSat <- v()
densTypesSat <- list()
densTestSat <- list()
for(ii in 1:nrow(testSat)){
  x <- t.default(as.matrix(testSat[ii,-37])) 
  for(t in typesSat){
    for (i in 1:R){
      J <- fitModelsSat[[t]][[i]]$fitinfo$K
      detJ <- fitModelsSat[[t]][[i]]$fitinfo$detK
      h <- J%*%muTrainSat[[t]]
      densModelSat[i] <- sqrt(detJ)*exp(t.default(h)%*%x - 0.5*t.default(x)%*%J%*%x -0.5*t.default(h)%*%muTrainSat[[t]])
    }
    densTypesSat[[t]] <- sum(densModelSat)/R
  }
  densTestSat[[ii]] <- densTypesSat
  cat(ii,"\n")
}

# Estimating the posterior densities where the probablity of belonging 
# to one of the six classes is equal.
postTypesSat <- list()
postTestSat <- list()
for(ii in 1:nrow(testSat)){
  for(t in typesSat){
    postTypesSat[[t]] <- densTestSat[[ii]][[t]]/do.call(sum,densTestSat[[ii]])
  }
  postTestSat[[ii]] <- postTypesSat
}

# Verifying that the posterior densities for each observation summarise to one. 
sumPostSat <- v()
for(ii in 1:nrow(testSat)){
  sumPostSat[ii] <- do.call(sum,postTestSat[[ii]])
}
sumPostSat

# The estimated classification of the test set when using the HLM. 
predTest <- v()
for(ii in 1:nrow(testSat)){
  predTest[ii] <- which(postTestSat[[ii]]==do.call(max,postTestSat[[ii]]))
}
predTest <- factor(predTest, label=typesSat)

# The following results may be slightly different than the results in the project 
# due to the randomness in the HLM.
confMatrixHLM <- table(testSat[,37],predTest) #HLM (Headlong method)
confMatrixHLM
diag(prop.table(confMatrixHLM, 1))
accSatHLm <- sum(diag(prop.table(confMatrixHLM)))
accSatHLm

#Standard error and confidence interval.
SESatHLM <- sqrt((accSatHLm*(1-accSatHLm)/nrow(testSat)))
c(accSatHLm - 2*SESatHLM, accSatHLm + 2*SESatHLM)
