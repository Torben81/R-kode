# In the follwing we make an example of model selection with the HLM compared with gRim
d <- trainSat[trainSat[,37]=="red soil",-37]
m0 <- init(d)
set.seed(5)
system.time(
  m <- step(m0)
)
m$it
#374 iterations in 0.28 seconds
fitm <- cmod(edgeL(m),d)
system.time(
  forwm <- forward(fitm,criterion="test",alpha=.01,type="unrestricted",search="headlong", steps=1)
)

# In the following we make Classification on the Satellite data using the HLM 
types <- levels(Satellite$classes)

models <- list()
typesModels <- list()
system.time(
for(t in types){
  d <- trainSat[trainSat$classes==t, -37]
  m0 <- init(d)
  for(i in 1:100){
    models[[i]] <- step(m0)
    cat(i,"\n")
  }
  typesModels[[t]] <- models
})

modelsEdgeL <- list()
typesEdgeL <- list()
for(t in types){
  for(i in 1:100){
    modelsEdgeL[[i]] <- edgeL(typesModels[[t]][[i]])
  }
  typesEdgeL[[t]] <- modelsEdgeL
}

fitModels <- list()
fitModelsTypes <- list()
system.time(
  for(t in types){
    d <- trainSat[trainSat$classes==t, -37]
    for(i in 1:100){
      fitModels[[i]] <- cmod(typesEdgeL[[t]][[i]],d)
      cat(i,"\n")
    }
    fitModelsTypes[[t]] <- fitModels
  }
)

modelsDensity <- rep(0,100)
typesDensity <- list()
testDensity <- list()
system.time(
  for(ii in 1:nrow(testSat)){
    x <- t(as.matrix(testSat[ii,-37]))
    for(t in types){
      d <- trainSat[trainSat$classes==t, -37]
      mu <- as.matrix(apply(d,2,mean))
      for (i in 1:100){
        #fitm <- cmod(typesEdgeL[[t]][[i]],d)
        J <- fitModelsTypes[[t]][[i]]$fitinfo$K
        detJ <- fitModelsTypes[[t]][[i]]$fitinfo$detK
        h <- J%*%mu
        modelsDensity[i] <- sqrt(detJ)*exp(t(h)%*%x - 0.5*t(x)%*%J%*%x -0.5*t(h)%*%mu)
      }
      typesDensity[[t]] <- modelsDensity
    }
    testDensity[[ii]] <- typesDensity
    cat(ii,"\n")
  }
)

avDenTypes <- rep(0,6)
names(avDenTypes) <- types
testTypes <- list()
for(ii in 1:nrow(testSat)){
  for(t in types){
    avDenTypes[t] <- sum(testDensity[[ii]][[t]])/100
  }
  testTypes[[ii]] <- avDenTypes
}

sumTest <- rep(0,nrow(testSat))
for(ii in 1:nrow(testSat)){
  sumDensity <- 0
  for(t in 1:6){
    sumDensity <- sumDensity + unname(testTypes[[ii]][t]/sum(testTypes[[ii]]))
    print(sumDensity)
  }
  sumTest[ii] <- sumDensity
}
sumDensity
sumTest

predTest <- rep(0,nrow(testSat))
for(ii in 1:nrow(testSat)){
  predTest[ii] <- which(testTypes[[ii]]/sum(testTypes[[ii]])==max(testTypes[[ii]]/sum(testTypes[[ii]])))
}
predTest <- factor(predTest, label=types)

confMatrixHLM <- table(testSat[,37],predTest) #HLM (Headlong method)
confMatrixHLM
diag(prop.table(confMatrixHLM, 1))
accSatHLm <- sum(diag(prop.table(confMatrixHLM)))

#Standard error
SESatHLM <- sqrt((accSatHLm*(1-accSatHLm)/nrow(testSat)))
c(accSatHLm - 2*SESatHLM, accSatHLm + 2*SESatHLM)
