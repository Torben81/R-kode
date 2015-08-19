# To illustrate the LDA decision boundaries, we apply LDA to the first two principal components 
# of the Satellite data.

# Principal components analysis using correlation matrix.
pca <- princomp(Satellite[,-37], cor=T, prior=rep(1,6)/6) 
pcComp <- pca$scores
pcComp1 <- pcComp[,1] # principal component 1 scores 
pcComp2 <- pcComp[,2] # principal component 2 scores

# The cumulative proportion of variance explained by the first two principal components.
sum(pca$sdev[1:2]^2/sum(pca$sdev^2))

# Applying LDA to the first two principal components.
modelLDA <- lda(Satellite[,37] ~ pcComp1+pcComp2, prior=rep(1,6)/6)
modelLDA

# Calculate the proportion of trace.
# (modelLDA$svd)^2/sum(modelLDA$svd^2)

# Predict the class using the estimated modelLDA, and calculating the confusion matrix.
predictPCALDA <- predict(modelLDA)$class
cofMatPCALDA <- table(Satellite[,37], predictPCALDA)
diag(prop.table(cofMatPCALDA, 1))     # The proportion of correct predicted observations within each class.
sum(diag(prop.table(cofMatPCALDA)))   # The proportion of correct predicted observations over all classes.

# This function estimate the LDA decision boundaries.
ldaBoundary <- function(mu.l,mu.k,sigma,pi.l=0.5,pi.k=0.5){
  isigma <- solve(sigma)
  b <- as.numeric(0.5*matrix(mu.l-mu.k,nrow=1)%*%isigma%*%matrix(mu.l+mu.k,ncol=1) - log(pi.k/pi.l))
  A <- as.numeric(matrix(mu.l-mu.k,ncol=2)%*%isigma)
  c(intercept=b/A[2],slope=-A[1]/A[2])
}

# This function estimate the sum of all the elements in a list. 
listSum <- function(x){
  stopifnot(is.list(x))
  if(length(x)==2){
    stopifnot(class(x[[1]])==class(x[[2]]))
    return(x[[1]]+x[[2]])
  }
  else listSum(c(list(x[[1]]+x[[2]]),x[-(1:2)]))
}

# A list contaning the empirical mean for each class.
musSat <- lapply(split(as.data.frame(pcComp[,1:2]),Satellite[,37]),colMeans)

# A list containing the empirical covariance matrix for each class.
sigmaSat <- lapply(split(as.data.frame(pcComp[,1:2]),Satellite[,37]),function(x)
  listSum(lapply(split(sweep(x,2,colMeans(x)),1:nrow(x)),function(y) as.vector(unlist(y))%*%t(as.vector(unlist(y))))))

tableSat <- table(Satellite[,37])  # Number of observations within each class.
SigmaSat <- listSum(sigmaSat)/(sum(tableSat)-6)  # The unbiased empirical covariance matrix.

# Linaer Discriminant Analysis on the first two principal components for class k and l.
bound12 <- ldaBoundary(mu.l=musSat[[1]],mu.k=musSat[[2]],sigma=SigmaSat)
bound13 <- ldaBoundary(mu.l=musSat[[1]],mu.k=musSat[[3]],sigma=SigmaSat)
bound14 <- ldaBoundary(mu.l=musSat[[1]],mu.k=musSat[[4]],sigma=SigmaSat)
bound15 <- ldaBoundary(mu.l=musSat[[1]],mu.k=musSat[[5]],sigma=SigmaSat)
bound25 <- ldaBoundary(mu.l=musSat[[2]],mu.k=musSat[[5]],sigma=SigmaSat)
bound34 <- ldaBoundary(mu.l=musSat[[3]],mu.k=musSat[[4]],sigma=SigmaSat)
bound46 <- ldaBoundary(mu.l=musSat[[4]],mu.k=musSat[[6]],sigma=SigmaSat)
bound56 <- ldaBoundary(mu.l=musSat[[5]],mu.k=musSat[[6]],sigma=SigmaSat)

# LDA -- boundary intersection points.
x0Sat134 <- (bound34[1]-bound13[1])/(bound13[2]-bound34[2])
y0Sat134 <- bound34[1] + bound34[2]*x0Sat134

x0Sat125 <- (bound12[1]-bound25[1])/(bound25[2]-bound12[2])
y0Sat125 <- bound12[1] + bound12[2]*x0Sat125

x0Sat145 <- (bound15[1]-bound14[1])/(bound14[2]-bound15[2])
y0Sat145 <- bound15[1] + bound15[2]*x0Sat145

x0Sat456 <- (bound46[1]-bound56[1])/(bound56[2]-bound46[2])
y0Sat456 <- bound46[1] + bound46[2]*x0Sat456

# Plotting the LDA decision boundaries.
cols <- c(2,3,1,4:6)   # Colors of the classes.
pdf("LDAPCA.pdf",height=7,width=9)
ggplot()+geom_point(aes(pcComp1,pcComp2,colour=Satellite[,37]),shape=1)+
  geom_line(aes(c(-10,x0Sat134),c(bound13[1]+bound13[2]*-10,y0Sat134)),colour=1,size=0.8)+
  geom_line(aes(c(-10,x0Sat134),c(bound13[1]+bound13[2]*-10,y0Sat134)),colour=2,size=0.8,linetype=2)+
  geom_line(aes(c(-2.5,x0Sat134),c(bound34[1]+bound34[2]*-2.5,y0Sat134)),colour=1,size=0.8)+
  geom_line(aes(c(-2.5,x0Sat134),c(bound34[1]+bound34[2]*-2.5,y0Sat134)),linetype=2,colour=4,size=0.8)+
  geom_line(aes(c(x0Sat134,x0Sat145),c(bound14[1]-bound14[2]*-x0Sat134,y0Sat145)),size=0.8, colour=2)+
  geom_line(aes(c(x0Sat134,x0Sat145),c(bound14[1]-bound14[2]*-x0Sat134,y0Sat145)),size=0.8,linetype=2,colour=4)+
  geom_line(aes(c(-10,x0Sat125),c(bound12[1]+bound12[2]*-10,y0Sat125)),size=0.8,colour=3)+
  geom_line(aes(c(-10,x0Sat125),c(bound12[1]+bound12[2]*-10,y0Sat125)),size=0.8,linetype=2,colour=2)+
  geom_line(aes(c(x0Sat145,x0Sat125),c(bound15[1]-bound15[2]*-x0Sat145,y0Sat125)),size=0.8, colour=2)+
  geom_line(aes(c(x0Sat145,x0Sat125),c(bound15[1]-bound15[2]*-x0Sat145,y0Sat125)),size=0.8,linetype=2,colour=5)+
  geom_line(aes(c(10,x0Sat125),c(bound25[1]-bound25[2]*-10,y0Sat125)),size=0.8, colour=3)+
  geom_line(aes(c(10,x0Sat125),c(bound25[1]-bound25[2]*-10,y0Sat125)),size=0.8,linetype=2,colour=5)+
  geom_line(aes(c(0.7,x0Sat456),c(bound46[1]-bound46[2]*-0.7,y0Sat456)),size=0.8, colour=4)+
  geom_line(aes(c(0.7,x0Sat456),c(bound46[1]-bound46[2]*-0.7,y0Sat456)),size=0.8,linetype=2,colour=6)+
  geom_line(aes(c(10,x0Sat456),c(bound56[1]-bound56[2]*-10,y0Sat456)),size=0.8, colour=5)+
  geom_line(aes(c(10,x0Sat456),c(bound56[1]-bound56[2]*-10,y0Sat456)),size=0.8,linetype=2,colour=6)+
  scale_colour_manual(name = "Model", values = cols) +
  geom_polygon(aes(x=c(-10,x0Sat134,x0Sat145,x0Sat125,-10),y=c(bound13[1]-bound13[2]*10,y0Sat134,y0Sat145,y0Sat125,bound12[1]-bound12[2]*10)), alpha=0.2, fill=2)+
  geom_polygon(aes(x=c(-10,x0Sat125,10,10,-10),y=c(bound12[1]-bound12[2]*10,y0Sat125,bound25[1]+bound25[2]*10,-15.5,-15.5)),fill=3, alpha=0.2)+
  geom_polygon(aes(x=c(-10,x0Sat134,-2.5,-10),y=c(bound13[1]-bound13[2]*10,y0Sat134,bound34[1]+bound34[2]*-2.5,5.4)),fill=1, alpha=0.2)+
  geom_polygon(aes(x=c(0.7,x0Sat456,x0Sat145,x0Sat134,-2.5),y=c(bound46[1]-bound46[2]*-0.7,y0Sat456,y0Sat145,y0Sat134,bound34[1]+bound34[2]*-2.5)),fill=4,alpha=0.2)+
  geom_polygon(aes(x=c(10,x0Sat125,x0Sat145,x0Sat456,10),y=c(bound25[1]-bound25[2]*-10,y0Sat125,y0Sat145,y0Sat456,bound56[1]+bound56[2]*10)),fill=5,alpha=0.2)+
  geom_polygon(aes(x=c(10,x0Sat456,0.7,10),y=c(bound56[1]-bound56[2]*-10,y0Sat456,bound46[1]+bound46[2]*0.7,5.5)),fill=6,alpha=0.2)+
  xlab("PC1") + ylab("PC2") +
  theme(plot.background = element_rect(fill="transparent", color=NA), legend.background=element_rect(fill="transparent"))
dev.off()

# Applying LDA to the Satellite data including all 36 variabels. 
ldaSat <- lda(classes~.,trainSat)
ldaSat$prior
predictLDA <- predict(ldaSat, testSat)
confMatrixLDA <- table(testSat[,37], predictLDA$class)
confMatrixLDA
diag(prop.table(confMatrixLDA, 1))
ldaSatCorrect <- sum(diag(prop.table(confMatrixLDA))) 
ldaSatCorrect

# Standard error and confidence interval.
SEldaSat <- sqrt((ldaSatCorrect*(1-ldaSatCorrect)/nrow(testSat)))
c(ldaSatCorrect - 2*SEldaSat, ldaSatCorrect + 2*SEldaSat)

# Histogram of the values of the discriminant functions for the samples from different classes.
#ldahist(predictLDA$x[,1], g=testSat[,37])
#ldahist(predictLDA$x[,2], g=testSat[,37])

# Applying QDA to the Satellite data including all 36 variabels.     
qdaSat <- qda(classes~., trainSat)
qdaSat$prior
predictQDA <- predict(qdaSat, testSat)
confMatrixQDA <- table(testSat[,37], predictQDA$class)
confMatrixQDA
diag(prop.table(confMatrixQDA, 1))
qdaSatCorrect <- sum(diag(prop.table(confMatrixQDA)))
qdaSatCorrect

#Standard error and confidence interval.
SEqdaSat <- sqrt((qdaSatCorrect*(1-qdaSatCorrect)/nrow(testSat)))
c(qdaSatCorrect - 2*SEqdaSat, qdaSatCorrect + 2*SEqdaSat)


