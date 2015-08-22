# Model selection on the Breastcancer data using the HLM.
controlBC <- trainBC[trainBC[,1001]=="control",-1001]
initBC <- init(controlBC) 
system.time(
  modelBC <- step(initBC)
)
modelBC$it

# In following we perform classification on the Breastcancer data using the HLM.
# Since the number of variabels is large they are divided into several parts.
parts <- 50  # Number of parts the variabels are divided into.
colSize <- ncol(trainBC[,-1001])/parts  # Number of variabels in each part (must be an integer).
set.seed(5)
colSample <- sample(rep(1:parts, colSize))
trainParts <- split(colnames(trainBC[,-1001]), colSample)
testParts <- split(colnames(testBC[,-1001]), colSample)

classTrain <- trainBC[,1001]
classTest <- testBC[,1001]

trainColParts <- list()
testColParts <- list()
for(i in 1:parts){  # Training and test sets both contaning the same parts of variabels.
  trainColParts[[i]] <- cbind(trainBC[,trainParts[[i]]],classTrain)
  testColParts[[i]] <- cbind(testBC[,testParts[[i]]],classTest)
}

R <- 100   # To obtain so sort of Random forest the HLM is repited R times on the training data.
typesBC <- levels(trainBC$code)
muPartsBC <- list()
muPartsTrainBC <- list()
fitmBC <- list()
fitmTypesBC <- list()
fitmColPartsBC <- list()
for(i in 1:parts){  # Performing model selection R times on every part of variabels for each class.
  for(t in typesBC){
    d <- trainColParts[[i]][trainColParts[[i]]$classTrain==t,-(colSize+1)]
    muPartsBC[[t]] <- as.matrix(apply(d,2, mean)) # Emperical mean of a part of variabels.
    for(j in 1:R){
      m0 <- init(d)
      mBC <- step(m0)
      mBCEdgeL <- edgeL(mBC, d)  # Extrating the edges from the model.
      fitmBC[[j]] <- cmod(mBCEdgeL,d)  # Fitting the model.
    }
    fitmTypesBC[[t]] <- fitmBC
  }
  fitmColPartsBC[[i]] <- fitmTypesBC
  muPartsTrainBC[[i]] <- muPartsBC
  cat(i,"\n")
}

# Estimating the conditional density of every observation in the test set given the fitted model above.
densModelBC <- c()
densTypesBC <- list()
densPartsBC <- list()
densTestBC <- list()
for(ii in 1:nrow(testBC)){  
  for(i in 1:parts){
    x <- t.default(as.matrix(testColParts[[i]][ii,-(colSize+1)]))
    for(t in typesBC){
      for(j in 1:R){
        J <- fitmColPartsBC[[i]][[t]][[j]]$fitinfo$K
        detJ <- fitmColPartsBC[[i]][[t]][[j]]$fitinfo$detK
        mu <- muPartsTrainBC[[i]][[t]]
        h <- J%*%mu
        densModelBC[j] <- sqrt(detJ)*exp(t.default(h)%*%x - 0.5*t.default(x)%*%J%*%x -0.5*t.default(h)%*%mu)
      }
      densTypesBC[[t]] <- densModelBC
    }
    densPartsBC[[i]] <- densTypesBC
  }
  densTestBC[[ii]] <- densPartsBC
  cat(ii,"\n")
}

# Estimating the posterior densities where the probablity of belonging 
# to one of the six classes is equal.
densControl <- 
densCase <- 
densPartsCase <- 
densPartsControl <- 
votesTestBC <- c()
postTestBC <- list()
for(ii in 1:125){
  votes <- 0
  for(i in 1:parts){
    for(j in 1:R){
      densCase[j] <- densTestBC[[ii]][[i]][[1]][j]
      densControl[j] <- densTestBC[[ii]][[i]][[2]][j]
      votes <- votes + 1*(densCase[j] > densControl[j])
    }
    densPartsCase[i] <- sum(densCase)/R
    densPartsControl[i] <- sum(densControl)/R
  }
  votesTestBC[ii] <- votes
  densCaseBC <- sum(densPartsCase)
  densControlBC <- sum(densPartsControl)
  postCaseBC <- densCaseBC/sum(densCaseBC,densControlBC)
  postControlBC <- densControlBC/sum(densCaseBC,densControlBC)
  postTestBC[[ii]] <- c(postCaseBC,postControlBC)
}

voteClass <- c()
for(ii in 1:nrow(testBC)){
  if(votesTestBC[ii]>2500){
    voteClass[ii] <- 1
  } else voteClass[ii] <- 2
}
voteClass <- factor(voteClass,labels=typesBC)

sumPostTestBC <- c()
for(ii in 1:nrow(testBC)){
  sumPostTestBC[ii] <- sum(postTestBC[[ii]])
}
sumPostTestBC

predBC <- c()
for(ii in 1:nrow(testBC)){
  predBC[ii] <- which(postTestBC[[ii]]==max(postTestBC[[ii]]))
}
predBC <- factor(predBC, label=typesBC)

# Confusion matrix, sensitivity, specificity and accuracy of the HLM classifier.
confMatrixBC <- table(testBC[,1001], predBC)
confMatrixBC
diag(prop.table(confMatrixBC, 1))
accBCHLM <- sum(diag(prop.table(confMatrixBC)))
accBCHLM

# Standard error #HLM and confidence interval.
SEBCHLM <- sqrt((accBCHLM*(1-accBCHLM)/nrow(testBC)))
c(accBCHLM - 2*SEBCHLM, accBCHLM + 2*SEBCHLM)

# Confusion matrix, sensitivity, specificity and accuracy of the HLMvotes classifier.
confMatrixBCVotes <- table(testBC[,1001], voteClass)
confMatrixBCVotes
diag(prop.table(confMatrixBCVotes,1))
accBCHLMVotes <- sum(diag(prop.table(confMatrixBCVotes)))
accBCHLMVotes

#Standard error and confidence interval.
SEBCHLMVotes <- sqrt((accBCHLMVotes*(1-accBCHLMVotes)/nrow(testBC)))
c(accBCHLMVotes - 2*SEBCHLMVotes, accBCHLMVotes + 2*SEBCHLMVotes)
