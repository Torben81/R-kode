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
for(i in 1:split){
  trainColSplit[[i]] <- cbind(trainBC[,trainSplit[[i]]],classTrain)
  testColSplit[[i]] <- cbind(testBC[,testSplit[[i]]],classTest)
}

R <- 100   # To obtain so sort of Random forest the HLM is repited R times on the training data.
types <- levels(trainBC$code)
mu <- list()
muSplit <- list()
dataInfo <- list()
modelsTypes <- list()
modelsColSplit <- list()
modelsEdgeL <- list()
edgeLColSplit <- list()
for(j in 1:R){
  for(i in 1:split){
    for(t in types){
      d <- trainColSplit[[i]][trainColSplit[[i]]$classTrain==t,-(colSize+1)]
      mu[[t]] <- as.matrix(apply(d,2, mean))
      m0 <- init(d)
      modelsTypes[[t]] <- step(m0)
      modelsEdgeL[[t]] <- edgeL(modelsTypes[[t]], d)
    }
    modelsColSplit[[i]] <- modelsTypes
    edgeLColSplit[[i]] <- modelsEdgeL
    muSplit[[i]] <- mu
  }
  dataInfo[[j]] <- list(modelsColSplit, edgeLColSplit, muSplit)
  cat(j,"\n")
}

fitmTypes <- list()
fitmColSplit <- list()
fitmodels <- list()
for(j in 1:R){
  for(i in 1:split){
    for(t in types){
      d <- trainColSplit[[i]][trainColSplit[[i]]$classTrain==t,-(colSize+1)]
      fitmTypes[[t]] <- cmod(dataInfo[[j]][[2]][[i]][[t]],d)
    }
    fitmColSplit[[i]] <- fitmTypes
  }
  fitmodels[[j]] <- fitmColSplit
  cat(j,"\n")
}

dens <- c(0,0)
names(dens) <- c("case", "control")
densCol <- list()
densModels <- list()
densObs <- list()
for(ii in 1:nrow(testBC)){
  for(j in 1:R){
    for(i in 1:split){
      x <- t(as.matrix(testColSplit[[i]][ii,-(colSize+1)]))
      for(t in types){
        J <- fitmodels[[j]][[i]][[t]]$fitinfo$K
        detJ <- fitmodels[[j]][[i]][[t]]$fitinfo$detK
        mu <- dataInfo[[j]][[3]][[i]][[t]]
        h <- J%*%mu
        dens[t] <- sqrt(detJ)*exp(t(h)%*%x - 0.5*t(x)%*%J%*%x -0.5*t(h)%*%mu)
      }
      densCol[[i]] <- dens
    }
    densModels[[j]] <- densCol
  }
  densObs[[ii]] <- densModels
  cat(ii,"\n")
}

densSplitControl <- rep(0,split)
densSplitCase <- rep(0,split)
avSplitControl <- rep(0,R)
avSplitCase <- rep(0,R)
avDensCaseBC <- rep(0,125)
avDensControlBC <- rep(0,125)
votesObs <- rep(0,125)
predCaseBC <- rep(0,125)
predControlBC <- rep(0,125)
predCaseControlBC <- list()

for(ii in 1:125){
  votes <- 0
  for(j in 1:R){
    for(i in 1:split){
      densSplitCase[i] <- densObs[[ii]][[j]][[i]][1]
      densSplitControl[i] <- densObs[[ii]][[j]][[i]][2]
      votes <- votes + 1*(densSplitCase[i] > densSplitControl[i])
    }
    avSplitCase[j] <- sum(densSplitCase)
    avSplitControl[j] <- sum(densSplitControl)
  }
  votesObs[ii] <- votes
  avDensCaseBC[ii] <- sum(avSplitCase)/R
  avDensControlBC[ii] <- sum(avSplitControl)/R
  predCaseBC[ii] <- avDensCaseBC[ii]/sum(avDensCaseBC[ii],avDensControlBC[ii])
  predControlBC[ii] <- avDensControlBC[ii]/sum(avDensCaseBC[ii],avDensControlBC[ii])
  predCaseControlBC[[ii]] <- c(predCaseBC[ii],predControlBC[ii])
}

voteClass <- rep(0,125)
for(ii in 1:nrow(testBC)){
  if(votesObs[ii]>2500){
    voteClass[ii] <- 1
  } else voteClass[ii] <- 2
}
voteClass <- factor(voteClass,labels=types)

sumCaseControl <- rep(0,125)
for(ii in 1:nrow(testBC)){
  sumCaseControl[ii] <- sum(predCaseBC[ii],predControlBC[ii])
}

predBC <- rep(0,nrow(testBC))
for(ii in 1:nrow(testBC)){
  predBC[ii] <- which(predCaseControlBC[[ii]]==max(predCaseControlBC[[ii]]))
}
predBC <- factor(predBC, label=types)

# Confusion matrix, sensitivity, specificity and accuracy of the HLM classifier.
confMatrixBC <- table(testBC[,1001], predBC)
confMatrixBC
diag(prop.table(confMatrixBC, 1))
accBCHLM <- sum(diag(prop.table(confMatrixBC)))

# Standard error #HLM and confidence interval.
SEBCHLM <- sqrt((accBCHLM*(1-accBCHLM)/nrow(testBC)))
c(accBCHLM - 2*SEBCHLM, accBCHLM + 2*SEBCHLM)

# Confusion matrix, sensitivity, specificity and accuracy of the HLMvotes classifier.
confMatrixBCVotes <- table(testBC[,1001], voteClass)
confMatrixBCVotes
diag(prop.table(confMatrixBCVotes,1))
accBCHLMVotes <- sum(diag(prop.table(confMatrixBCVotes)))

#Standard error and confidence interval.
SEBCHLMVotes <- sqrt((accBCHLMVotes*(1-accBCHLMVotes)/nrow(testBC)))
c(accBCHLMVotes - 2*SEBCHLMVotes, accBCHLMVotes + 2*SEBCHLMVotes)
