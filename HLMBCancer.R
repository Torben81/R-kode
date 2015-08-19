# Model selection on the Breastcancer data using the HLM.
controlBC <- trainBC[trainBC[,1001]=="control",-1001]
initBC <- init(controlBC) 
system.time(
  mBC <- step(initBC)
)
mBC$it

# In following we perform classification on the Breastcancer data using the HLM.
# Since the number of variabels is large they are divided into several parts.
split <- 50  # Number of parts the variabels are divided into.
colSize <- ncol(trainBC[,-1001])/split  # Number of variabels in each part (must be an integer).
set.seed(5)
colSample <- sample(rep(1:split, colSize))
trainSplit <- split(colnames(trainBC[,-1001]), colSample)
testSplit <- split(colnames(testBC[,-1001]), colSample)

classTrain <- trainBC[,1001]
classTest <- testBC[,1001]

trainColSplit <- list()
testColSplit <- list()
for(i in 1:split){  # Training and test sets both contaning the same parts of variabels.
  trainColSplit[[i]] <- cbind(trainBC[,trainSplit[[i]]],classTrain)
  testColSplit[[i]] <- cbind(testBC[,testSplit[[i]]],classTest)
}

R <- 100   # To obtain so sort of Random forest the HLM is repited R times on the training data.
typesBC <- levels(trainBC$code)
muSplitBC <- list()
muSplitTrainBC <- list()
fitmBC <- list()
fitmTypesBC <- list()
fitmColSplitBC <- list()
for(i in 1:split){  # Performing model selection R times on every part of variabels for each class.
  for(t in typesBC){
    d <- trainColSplit[[i]][trainColSplit[[i]]$classTrain==t,-(colSize+1)]
    muSplitBC[[t]] <- as.matrix(apply(d,2, mean)) # Emperical mean of a part of variabels.
    for(j in 1:R){
      m0 <- init(d)
      mBC <- step(m0)
      mBCEdgeL <- edgeL(mBC, d)  # Extrating the edges from the model.
      fitmBC[[j]] <- cmod(mBCEdgeL,d)  # Fitting the model.
    }
    fitmTypesBC[[t]] <- fitmBC
  }
  fitmColSplitBC[[i]] <- fitmTypesBC
  muSplitTrainBC[[i]] <- muSplitBC
  cat(i,"\n")
}

# Estimating the conditional density of every observation in the test set given the fitted model above.
densModelBC <- c()
densTypesBC <- list()
densSplitBC <- list()
densTestBC <- list()
for(ii in 1:nrow(testBC)){  
  for(i in 1:split){
    x <- t.default(as.matrix(testColSplit[[i]][ii,-(colSize+1)]))
    for(t in typesBC){
      for(j in 1:R){
        J <- fitmColSplitBC[[i]][[t]][[j]]$fitinfo$K
        detJ <- fitmColSplitBC[[i]][[t]][[j]]$fitinfo$detK
        mu <- muSplitTrainBC[[i]][[t]]
        h <- J%*%mu
        densModelBC[j] <- sqrt(detJ)*exp(t.default(h)%*%x - 0.5*t.default(x)%*%J%*%x -0.5*t.default(h)%*%mu)
      }
      densTypesBC[[t]] <- densModelBC
    }
    densSplitBC[[i]] <- densTypesBC
  }
  densTestBC[[ii]] <- densSplitBC
  cat(ii,"\n")
}

# Estimating the posterior densities where the probablity of belonging 
# to one of the six classes is equal.
densControl <- 
densCase <- 
densSplitCase <- 
densSplitControl <- 
votesTestBC <- c()
postTestBC <- list()
for(ii in 1:125){
  votes <- 0
  for(i in 1:split){
    for(j in 1:R){
      densCase[j] <- densTestBC[[ii]][[i]][[1]][j]
      densControl[j] <- densTestBC[[ii]][[i]][[2]][j]
      votes <- votes + 1*(densCase[j] > densControl[j])
    }
    densSplitCase[i] <- sum(densCase)/R
    densSplitControl[i] <- sum(densControl)/R
  }
  votesTestBC[ii] <- votes
  densCaseBC <- sum(densSplitCase)
  densControlBC <- sum(densSplitControl)
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
