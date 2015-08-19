typesBC <- levels(BC$code)
muTrainBC <- list()
glassoTrainBC <- list()
for(t in typesBC){
  d <- trainBC[trainBC[,1001]==t,-1001]
  muTrainBC[[t]] <- apply(d,2,mean)
  covTrainBC <- cor(d)
  glassoTrainBC[[t]] <- glasso(covTrainBC, rho=0.4)
  cat(t,"\n")
}

densGlassoBC <- c()
densGlassoTestBC <- list()
for(ii in 1:nrow(testBC)){
  x <- t.default(as.matrix(testBC[ii,-1001]))
  for(t in typesBC){
    mu <- muTrainBC[[t]]
    J <- glassoTrainBC[[t]]$wi
    detJ <- det(glassoTrainBC[[t]]$wi)
    h <- J%*%mu
    densGlassoBC[t] <- sqrt(detJ)*exp(t.default(h)%*%x - 0.5*t.default(x)%*%J%*%x -0.5*t.default(h)%*%mu)
  }
  densGlassoTestBC[[ii]] <- densGlassoBC
  cat(ii,"\n")
}

postTypesGlasso <- c()
postGlasso <- list()
for(ii in 1:nrow(testBC)){
  for(t in typesBC){
  postTypesGlasso[t] <- densGlassoTestBC[[ii]][t]/sum(densGlassoTestBC[[ii]])
  }
  postGlasso[[ii]] <- postTypesGlasso
}

sumPostGlasso <- c()
for(ii in 1:nrow(testBC)){
    sumPostGlasso[ii] <- sum(postGlasso[[ii]])
}
sumPostGlasso

predglassoBC <- c()
for(ii in 1:nrow(testBC)){
  predglassoBC[ii] <- which(postGlasso[[ii]]==max(postGlasso[[ii]]))
}
predglassoBC <- factor(predglassoBC, labels=typesBC)

glassoconfM <- table(testBC[,1001], predglassoBC)
glassoconfM
diag(prop.table(glassoconfM, 1))
accGlasso <- sum(diag(prop.table(glassoconfM)))
accGlasso

#Standard error
SEglasso <- sqrt((accGlasso*(1-accGlasso)/nrow(testBC)))
c(accGlasso - 2*SEglasso, accGlasso + 2*SEglasso)
