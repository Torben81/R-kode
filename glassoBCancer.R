types <- levels(trainBC$code)
muBC <- list()
glassotrainBC <- list()
for(t in types){
  d <- trainBC[trainBC[,1001]==t,-1001]
  muBC[[t]] <- apply(d,2,mean)
  covtrainBC <- cor(d)
  system.time(glassotrainBC[[t]] <- glasso(covtrainBC, rho=0.4))
}
glassotrainBC[[1]]
covtrainBC[1:10,1:10] 
names(glassotrainBC[[1]])

densTypesTestBC <- c(0,0)
names(densTypesTestBC) <- types
densTestBC <- list()
for(ii in 1:nrow(testBC)){
  x <- t(as.matrix(testBC[ii,-1001]))
  for(t in types){
    mu <- muBC[[t]]
    J <- glassotrainBC[[t]]$wi
    detJ <- det(glassotrainBC[[t]]$wi)
    h <- J%*%mu
    densTypesTestBC[t] <- sqrt(detJ)*exp(t(h)%*%x - 0.5*t(x)%*%J%*%x -0.5*t(h)%*%mu)
  }
  densTestBC[[ii]] <- densTypesTestBC
}
densTestBC[[1]][1]/sum(densTestBC[[1]])
testBC[51,1001]
classes <- testBC[,1001]

caseglasso <- rep(0,125)
controlglasso <- rep(0,125)
classglasso <- list()
for(ii in 1:nrow(testBC)){
  caseglasso[ii] <- densTestBC[[ii]][1]/sum(densTestBC[[ii]])
  controlglasso[ii] <- densTestBC[[ii]][2]/sum(densTestBC[[ii]])
  classglasso[[ii]] <- c(caseglasso[ii], controlglasso[ii])
}

sumTestBC <- rep(0,nrow(testBC))
for(ii in 1:nrow(testBC)){
  sumDenBC <- 0
  for(t in types){
    sumDenBC <- sumDenBC + unname(densTestBC[[ii]][t]/sum(densTestBC[[ii]]))
  }
  sumTestBC[ii] <- sumDenBC
}
sumDenBC
sumTestBC

predglassoBC <- rep(0,nrow(testBC))
for(ii in 1:nrow(testBC)){
  predglassoBC[ii] <- which(classglasso[[ii]]==max(classglasso[[ii]]))
}
predglassoBC <- factor(predglassoBC, labels=types)

glassoconfM <- table(testBC[,1001], predglassoBC)
glassoconfM
(13+85)/125
diag(prop.table(glassoconfM, 1))
accGlasso <- sum(diag(prop.table(glassoconfM)))

#Standard error
SEglasso <- sqrt((accGlasso*(1-accGlasso)/nrow(testBC)))
c(accGlasso - 2*SEglasso, accGlasso + 2*SEglasso)
