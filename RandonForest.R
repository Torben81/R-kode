source("forBack.R")

set.seed(123)
sampleSize <- sample(nrow(Satellite),3217)
train <- Satellite[sampleSize,]
test <- Satellite[-sampleSize,]

types <- levels(Satellite$classes)

models <- list()
allModels <- list()
system.time(
for(t in types){
  d <- train[train$classes==t, -37]
  m0 <- init(d)
  for(i in 1:100){
    models[[i]] <- step(m0)
    cat(i,"\n")
  }
  allModels[[t]] <- models
})

modelsEdgeL <- list()
typesEdgeL <- list()
for(t in types){
  for(i in 1:100){
    modelsEdgeL[[i]] <- edgeL(allModels[[t]][[i]])
  }
  typesEdgeL[[t]] <- modelsEdgeL
}

x <- t(as.matrix(test[1,-37]))


modelsDensity <- rep(0,100)
typesDensity <- list()
system.time(
  for(t in types){
    d <- train[train$classes==t, -37]
    mu <- as.matrix(apply(d,2,mean))
    for (i in 1:100){
      fitm <- cmod(typesEdgeL[[t]][[i]],d)
      J <- fitm$fitinfo$K
      detJ <- fitm$fitinfo$detK
      h <- J%*%mu
      modelsDensity[i] <- sqrt(detJ)*exp(t(h)%*%x - 0.5*t(x)%*%J%*%x -0.5*t(h)%*%mu)
      cat(i,"\n")
    }
    typesDensity[[t]] <- modelsDensity
  })

avDenTypes <- rep(0,6)
names(avDenTypes) <- types
for(t in types){
  avDenTypes[t] <- sum(typesDensity[[t]])/100
}

sumDensity <- 0
for(i in 1:6){
  sumDensity <- sumDensity + unname(avDenTypes[i]/(sum(avDenTypes)))
  print(sumDensity)
}
sumDensity
