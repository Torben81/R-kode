typesBC <- levels(BC$code)

# Applying the glasso method to the traning Breastcancer data for each class.
muTrainBC <- list()
glassoTrainBC <- list()
for(t in typesBC){
  d <- trainBC[trainBC[,1001]==t,-1001]
  muTrainBC[[t]] <- apply(d,2,mean)
  corTrainBC <- cor(d)
  glassoTrainBC[[t]] <- glasso(corTrainBC, rho=0.4)
  cat(t,"\n")
}

# Estimating the conditional density of every observation in the test set given the fitted models above.
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

# Estimating the posterior densities where the probablity of belonging 
# to one of the six classes is equal.
postTypesGlasso <- c()
postGlasso <- list()
for(ii in 1:nrow(testBC)){
  for(t in typesBC){
  postTypesGlasso[t] <- densGlassoTestBC[[ii]][t]/sum(densGlassoTestBC[[ii]])
  }
  postGlasso[[ii]] <- postTypesGlasso
}

# Verifying that the posterior densities for each observation summarise to one. 
sumPostGlasso <- c()
for(ii in 1:nrow(testBC)){
    sumPostGlasso[ii] <- sum(postGlasso[[ii]])
}
sumPostGlasso

# The estimated classification of the test set when using glasso on the Breastcancer data.
predglassoBC <- c()
for(ii in 1:nrow(testBC)){
  predglassoBC[ii] <- which(postGlasso[[ii]]==max(postGlasso[[ii]]))
}
predglassoBC <- factor(predglassoBC, labels=typesBC)

# Estimating the confusion matrix, sensitivity, specificity and accuracy of the glasso classifier.
glassoconfM <- table(testBC[,1001], predglassoBC)
glassoconfM
diag(prop.table(glassoconfM, 1))
accGlasso <- sum(diag(prop.table(glassoconfM)))
accGlasso

# Standard error and confidence interval.
SEglasso <- sqrt((accGlasso*(1-accGlasso)/nrow(testBC)))
c(accGlasso - 2*SEglasso, accGlasso + 2*SEglasso)
