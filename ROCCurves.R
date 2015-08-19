# ROC curves for LDA applied on Satellite data.
postLDA <- predictLDA$posterior  # Posterior probabilties for LDA.
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

# All ROC curves for LDA in one plot.
cols <- c(1,3:4,2,5:6)   # Colors of the six classes.
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
  
#ROC curves for QDA applied on Satellite data.
postQDA <- predictQDA$posterior  # Posterior probabilties for QDA.
rocQDARS <- roc(testSat[,37]=="red soil", postQDA[,1])
rocQDARS
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

# All ROC curves for QDA in one plot. 
pdf("QDArocSatellite.pdf",7,5)
ggplot()+geom_path(aes(x=1-rocQDARS$specificities, y=rocQDARS$sensitivities, colour="Red soil"), size=0.8)+
  geom_path(aes(x=1-rocQDACC$specificities, y=rocQDACC$sensitivities, colour="Cotten corp"), size=0.8) + 
  geom_path(aes(x=1-rocQDAGS$specificities, y=rocQDAGS$sensitivities, colour="Grey soil"), size=0.8) + 
  geom_path(aes(x=1-rocQDADGS$specificities, y=rocQDADGS$sensitivities, colour="Damp grey soil"), size=0.8) +
  geom_path(aes(x=1-rocQDAVS$specificities, y=rocQDAVS$sensitivities, colour="Vegetation stubble"), size=0.8) +
  geom_path(aes(x=1-rocQDAVDGS$specificities, y=rocQDAVDGS$sensitivities, colour="Very damp grey soil"), size=0.8) +
  scale_colour_manual(name = "Model", values = cols) + xlab("1 - Specificity") + ylab("Sensitivity") +
  theme(plot.background = element_rect(fill="transparent", color=NA), legend.background=element_rect(fill="transparent"))
dev.off()

# ROC curve for the HLM applyed on the Satellite data.
postControlBC <- as.numeric(lapply(postTestBC,"[[",2))
rocBC <- roc(testBC[,1001], postControlBC, levels=c("control", "case"))
rocBC$auc

# ROC curve for HLMvotes.
rocBCVotes <-roc(testBC[,1001], votesTestBC, levels=c("control", "case"))
rocBCVotes$auc

# ROC curve for glasso.
postControlGlasso <- as.numeric(lapply(postGlasso,"[[",2))
rocGlasso <- roc(testBC[,1001], postControlGlasso, levels=c("control", "case"))
rocGlasso$auc

# All in one plot.
pdf("ROCAll.pdf",7,5)
ggplot()+geom_path(aes(x=1-rocGlasso$specificities, y=rocGlasso$sensitivities, colour="glasso"), size=0.8)+
  geom_path(aes(x=1-rocBC$specificities, y=rocBC$sensitivities, colour="HLM"), size=0.8) + 
  geom_path(aes(x=1-rocBCVotes$specificities, y=rocBCVotes$sensitivities, colour="HLM votes"), size=0.8) + 
  scale_colour_discrete(name = "Model") + xlab("1 - Specificity") + ylab("Sensitivity") +
  theme(plot.background = element_rect(fill="transparent", color=NA), legend.background=element_rect(fill="transparent"))
dev.off()