setwd("C:/Users/Torben/Dropbox/speciale/Rapport/images")
#glasso
rocglasso <- roc(testBC[,1001], controlglasso, levels=c("control", "case"))
pdf("ROCGlassoBC.pdf",5,5)
plot(rocglasso)
dev.off()


ggplot()+geom_path(aes(x=1-rocglasso$specificities, y=rocglasso$sensitivities, colour="glasso"), size=0.8)+
  #geom_path(aes(x=ROC_cox_test_benchmark$times, y=ROC_cox_test_benchmark$AUC, colour="Benchmark"), size=0.8) + 
  scale_colour_discrete(name = "Model") + xlab("1 - Specificity") + ylab("Sensitivity") +
  theme(plot.background = element_rect(fill="transparent", color=NA), legend.background=element_rect(fill="transparent"))



#Headlog
rocBC <- roc(testBC[,1001], predControlBC, levels=c("control", "case"))
rocBC$auc
pdf("ROCHeadlongBC.pdf",5,5)
plot(rocBC, col="red", main="ROC curve for headlong ")
dev.off()

#Votes
rocBCVotes <-roc(testBC[,1001], votesObs, levels=c("control", "case"))
rocBCVotes$auc
plot(rocBCVotes, col="red", main="ROC curve for headlong votes.")

#All in one plot
pdf("ROCAll.pdf",7,7)
plot(rocglasso)
plot(rocBC, add=TRUE, col="red")
plot(rocBCVotes, add=TRUE, col="blue")
legend("bottomright", legend=c("glasso", "headlong", "headlong votes"),col=c("black","blue","red"), lty=c(1,1,1))
dev.off()

pdf("ROCAll.pdf",7,5)
ggplot()+geom_path(aes(x=1-rocglasso$specificities, y=rocglasso$sensitivities, colour="glasso"), size=0.8)+
  geom_path(aes(x=1-rocBC$specificities, y=rocBC$sensitivities, colour="HLM"), size=0.8) + 
  geom_path(aes(x=1-rocBCVotes$specificities, y=rocBCVotes$sensitivities, colour="HLM votes"), size=0.8) + 
  scale_colour_discrete(name = "Model") + xlab("1 - Specificity") + ylab("Sensitivity") +
  theme(plot.background = element_rect(fill="transparent", color=NA), legend.background=element_rect(fill="transparent"))
dev.off()

names(rocVotes)
rocVotes$auc

#ROC curves for LDA applied on Satellite data
postLDA <- predictLDA$posterior
rocLDA <- roc(testSat[,37]=="red soil", postLDA[,1])
rocLDA
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
pdf("LDArocSatellite.pdf",7,7)
plot(rocLDA, col=2)
plot(rocLDACC, col=3, add=TRUE)
plot(rocLDAGS, add=TRUE)
plot(rocLDADGS, col=4, add=TRUE)
plot(rocLDAVS, col=5, add=TRUE)
plot(rocLDAVDGS, col=6, add=TRUE)
legend("bottomright", legend=levels(test[,37]), col=c(2,3,1,4:6), lty=rep(1,6))
dev.off()

# # Graph settings
# require(RColorBrewer)
# cols <- brewer.pal(7,"Set1")
# cols <- c(cols[2:3], cols[1], cols[4:5], cols[7]) # inverts red/blue order
cols <- c(1,3:4,2,5:6)

pdf("LDArocSatellite.pdf",7,5)
ggplot()+geom_path(aes(x=1-rocLDACC$specificities, y=rocLDACC$sensitivities, colour="Cotten corp"), size=0.8) +
  geom_path(aes(x=1-rocLDA$specificities, y=rocLDA$sensitivities, colour="Red soil"), size=0.8)+ 
  geom_path(aes(x=1-rocLDAGS$specificities, y=rocLDAGS$sensitivities, colour="Grey soil"), size=0.8) + 
  geom_path(aes(x=1-rocLDADGS$specificities, y=rocLDADGS$sensitivities, colour="Damp grey soil"), size=0.8) +
  geom_path(aes(x=1-rocLDAVS$specificities, y=rocLDAVS$sensitivities, colour="Vegetation stubble"), size=0.8) +
  geom_path(aes(x=1-rocLDAVDGS$specificities, y=rocLDAVDGS$sensitivities, colour="Very damp grey soil"), size=0.8) +
  xlab("1 - Specificity") + ylab("Sensitivity") +
  scale_colour_manual(name = "Model", values = cols) +
  theme(plot.background = element_rect(fill="transparent", color=NA), legend.background=element_rect(fill="transparent"))
dev.off()

head(postLDA)
pMaxi <- list()
for(i in 1:6){
  pMaxi[[i]] <- apply(postLDA[,-i],1,max)
}

altPostRS <- postLDA[,1]/(postLDA[,1] +pMaxi[[1]])
rocAltLDARS <- roc(test[,37]=="red soil", altPostRS)
plot(rocAltLDARS, col=2)

altPostCC <- postLDA[,2]/(postLDA[,2] +pMaxi[[2]])
rocAltLDACC <- roc(test[,37]=="cotton crop", altPostCC)
plot(rocAltLDACC, add=TRUE, col=3)

altPostGS <- postLDA[,3]/(postLDA[,3] +pMaxi[[3]])
rocAltLDAGS <- roc(test[,37]=="grey soil", altPostGS)
plot(rocAltLDAGS, add=TRUE, col=1)

altPostDGS <- postLDA[,4]/(postLDA[,4] +pMaxi[[4]])
rocAltLDADGS <- roc(test[,37]=="damp grey soil", altPostDGS)
plot(rocAltLDADGS, add=TRUE, col=4)

altPostVS <- postLDA[,5]/(postLDA[,5] +pMaxi[[5]])
rocAltLDAVS <- roc(test[,37]=="vegetation stubble", altPostVS)
plot(rocAltLDAVS, add=TRUE, col=5)

altPostVDGS <- postLDA[,6]/(postLDA[,6] +pMaxi[[6]])
rocAltLDAVDGS <- roc(test[,37]=="very damp grey soil", altPostVDGS)
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
pdf("QDArocSatellite.pdf",7,7)
plot(rocQDA, col=2)
plot(rocQDACC, col=3, add=TRUE)
plot(rocQDAGS, add=TRUE)
plot(rocQDADGS, col=4, add=TRUE)
plot(rocQDAVS, col=5, add=TRUE)
plot(rocQDAVDGS, col=6, add=TRUE)
legend("bottomright", legend=levels(test[,37]), col=c(2,3,1,4:6), lty=rep(1,6))
dev.off()

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
