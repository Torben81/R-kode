types <- levels(Satellite$classes)

models <- list()
typesModels <- list()
system.time(
for(t in types){
  d <- train[train$classes==t, -37]
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
    d <- train[train$classes==t, -37]
    for(i in 1:100){
      fitModels[[i]] <- cmod(typesEdgeL[[t]][[i]],d)
      cat(i,"\n")
    }
    fitModelsTypes[[t]] <- fitModels
  }
)

x <- t(as.matrix(test[1,-37]))

modelsDensity <- rep(0,100)
typesDensity <- list()
system.time(
  for(t in types){
    d <- train[train$classes==t, -37]
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
)

avDenTypes <- rep(0,6)
names(avDenTypes) <- types
for(t in types){
  avDenTypes[t] <- sum(typesDensity[[t]])/100
}

for(t in types){
  print(avDenTypes[t]/sum(avDenTypes))
}

sumTypes <- 0
for(t in types){
  sumTypes <- sumTypes + unname(avDenTypes[t]/sum(avDenTypes))
}

modelsDensity <- rep(0,100)
typesDensity <- list()
testDensity <- list()
system.time(
  for(ii in 1:nrow(test)){
    x <- t(as.matrix(test[ii,-37]))
    for(t in types){
      d <- train[train$classes==t, -37]
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
for(ii in 1:nrow(test)){
  for(t in types){
    avDenTypes[t] <- sum(testDensity[[ii]][[t]])/100
  }
  testTypes[[ii]] <- avDenTypes
}

sumTest <- rep(0,nrow(test))
for(ii in 1:nrow(test)){
  sumDensity <- 0
  for(t in 1:6){
    sumDensity <- sumDensity + unname(testTypes[[ii]][t]/sum(testTypes[[ii]]))
    print(sumDensity)
  }
  sumTest[ii] <- sumDensity
}
sumDensity
sumTest

predTest <- rep(0,nrow(test))
for(ii in 1:nrow(test)){
  predTest[ii] <- which(testTypes[[ii]]/sum(testTypes[[ii]])==max(testTypes[[ii]]/sum(testTypes[[ii]])))
}
predTest <- factor(predTest, label=types)

confMatrixNNM <- table(test[,37],predTest) #NNM (nearest neighbor method)
confMatrixNNM
diag(prop.table(confMatrixNNM, 1))
sum(diag(prop.table(confMatrixNNM)))

