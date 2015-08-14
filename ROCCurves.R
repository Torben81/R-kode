setwd("C:/Users/Torben/Dropbox/speciale/Rapport/images")
# ROC curve glasso
rocglasso <- roc(testBC[,1001], controlglasso, levels=c("control", "case"))


ggplot()+geom_path(aes(x=1-rocglasso$specificities, y=rocglasso$sensitivities, colour="glasso"), size=0.8)+
  #geom_path(aes(x=ROC_cox_test_benchmark$times, y=ROC_cox_test_benchmark$AUC, colour="Benchmark"), size=0.8) + 
  scale_colour_discrete(name = "Model") + xlab("1 - Specificity") + ylab("Sensitivity") +
  theme(plot.background = element_rect(fill="transparent", color=NA), legend.background=element_rect(fill="transparent"))



#ROC curvesHLM 
rocBC <- roc(testBC[,1001], predControlBC, levels=c("control", "case"))
rocBC$auc

#HLMVotes
rocBCVotes <-roc(testBC[,1001], votesObs, levels=c("control", "case"))
rocBCVotes$auc

#All in one plot
pdf("ROCAll.pdf",7,5)
ggplot()+geom_path(aes(x=1-rocglasso$specificities, y=rocglasso$sensitivities, colour="glasso"), size=0.8)+
  geom_path(aes(x=1-rocBC$specificities, y=rocBC$sensitivities, colour="HLM"), size=0.8) + 
  geom_path(aes(x=1-rocBCVotes$specificities, y=rocBCVotes$sensitivities, colour="HLM votes"), size=0.8) + 
  scale_colour_discrete(name = "Model") + xlab("1 - Specificity") + ylab("Sensitivity") +
  theme(plot.background = element_rect(fill="transparent", color=NA), legend.background=element_rect(fill="transparent"))
dev.off()


#ROC curves for LDA applied on Satellite data
postLDA <- predictLDA$posterior
rocLDARS <- roc(testSat[,37]=="red soil", postLDA[,1])
rocLDARS
rocLDACC <- roc(testSat[,37]=="cotton crop", postLDA[,2])
rocLDACC$auc
rocLDAGS <- roc(testSat[,37]=="grey soil", postLDA[,3])
rocLDAGS$auc
rocLDADGS <- roc(testSat[,37]=="damp grey soil", postLDA[,4])
rocLDADGS$auc
rocLDAVS <- roc(testSat[,37]=="vegetation stubble", postLDA[,5])
rocLDAVS$auc
rocLDAVDGS <- roc(testSat[,37]=="very damp grey soil", postLDA[,6])
rocLDAVDGS$auc

#All ROC curves for LDA in one plot 
cols <- c(1,3:4,2,5:6)   # Colors of the six classes
pdf("LDArocSatellite.pdf",7,5)
ggplot()+geom_path(aes(x=1-rocLDACC$specificities, y=rocLDACC$sensitivities, colour="Cotten corp"), size=0.8) +
  geom_path(aes(x=1-rocLDARS$specificities, y=rocLDARS$sensitivities, colour="Red soil"), size=0.8)+ 
  geom_path(aes(x=1-rocLDAGS$specificities, y=rocLDAGS$sensitivities, colour="Grey soil"), size=0.8) + 
  geom_path(aes(x=1-rocLDADGS$specificities, y=rocLDADGS$sensitivities, colour="Damp grey soil"), size=0.8) +
  geom_path(aes(x=1-rocLDAVS$specificities, y=rocLDAVS$sensitivities, colour="Vegetation stubble"), size=0.8) +
  geom_path(aes(x=1-rocLDAVDGS$specificities, y=rocLDAVDGS$sensitivities, colour="Very damp grey soil"), size=0.8) +
  xlab("1 - Specificity") + ylab("Sensitivity") +
  scale_colour_manual(name = "Model", values = cols) +
  theme(plot.background = element_rect(fill="transparent", color=NA), legend.background=element_rect(fill="transparent"))
dev.off()

pMaxi <- list()
for(i in 1:6){
  pMaxi[[i]] <- apply(postLDA[,-i],1,max)
}

altPostRS <- postLDA[,1]/(postLDA[,1] +pMaxi[[1]])
rocAltLDARS <- roc(testSat[,37]=="red soil", altPostRS)
plot(rocAltLDARS, col=2)

altPostCC <- postLDA[,2]/(postLDA[,2] +pMaxi[[2]])
rocAltLDACC <- roc(testSat[,37]=="cotton crop", altPostCC)
plot(rocAltLDACC, add=TRUE, col=3)

altPostGS <- postLDA[,3]/(postLDA[,3] +pMaxi[[3]])
rocAltLDAGS <- roc(testSat[,37]=="grey soil", altPostGS)
plot(rocAltLDAGS, add=TRUE, col=1)

altPostDGS <- postLDA[,4]/(postLDA[,4] +pMaxi[[4]])
rocAltLDADGS <- roc(testSat[,37]=="damp grey soil", altPostDGS)
plot(rocAltLDADGS, add=TRUE, col=4)

altPostVS <- postLDA[,5]/(postLDA[,5] +pMaxi[[5]])
rocAltLDAVS <- roc(testSat[,37]=="vegetation stubble", altPostVS)
plot(rocAltLDAVS, add=TRUE, col=5)

altPostVDGS <- postLDA[,6]/(postLDA[,6] +pMaxi[[6]])
rocAltLDAVDGS <- roc(testSat[,37]=="very damp grey soil", altPostVDGS)
plot(rocAltLDAVDGS, add=TRUE, col=6)
  
#ROC curves for QDA applied on Satellite data
postQDA <- predictQDA$posterior
rocQDA <- roc(testSat[,37]=="red soil", postQDA[,1])
rocQDA
rocQDACC <- roc(testSat[,37]=="cotton crop", postQDA[,2])
rocQDACC$auc
rocQDAGS <- roc(testSat[,37]=="grey soil", postQDA[,3])
rocQDAGS$auc
rocQDADGS <- roc(testSat[,37]=="damp grey soil", postQDA[,4])
rocQDADGS$auc
rocQDAVS <- roc(testSat[,37]=="vegetation stubble", postQDA[,5])
rocQDAVS$auc
rocQDAVDGS <- roc(testSat[,37]=="very damp grey soil", postQDA[,6])
rocQDAVDGS$auc

#All ROC curves for QDA in one plot 
pdf("QDArocSatellite.pdf",7,5)
ggplot()+geom_path(aes(x=1-rocQDA$specificities, y=rocQDA$sensitivities, colour="Red soil"), size=0.8)+
  geom_path(aes(x=1-rocQDACC$specificities, y=rocQDACC$sensitivities, colour="Cotten corp"), size=0.8) + 
  geom_path(aes(x=1-rocQDAGS$specificities, y=rocQDAGS$sensitivities, colour="Grey soil"), size=0.8) + 
  geom_path(aes(x=1-rocQDADGS$specificities, y=rocQDADGS$sensitivities, colour="Damp grey soil"), size=0.8) +
  geom_path(aes(x=1-rocQDAVS$specificities, y=rocQDAVS$sensitivities, colour="Vegetation stubble"), size=0.8) +
  geom_path(aes(x=1-rocQDAVDGS$specificities, y=rocQDAVDGS$sensitivities, colour="Very damp grey soil"), size=0.8) +
  scale_colour_manual(name = "Model", values = cols) + xlab("1 - Specificity") + ylab("Sensitivity") +
  theme(plot.background = element_rect(fill="transparent", color=NA), legend.background=element_rect(fill="transparent"))
dev.off()


#satellite headlong
densRed <- rep(0,3218)
for(ii in 1:3218){
  densRed[ii] <- testTypes[[ii]][6]/sum(testTypes[[ii]])
}
rocSat <- roc(test[,37], densRed)
plot(rocSat, col="red")
