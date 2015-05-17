d <- train[train[,37]=="vegetation stubble",-37]
init <- function(d,al=.01){#d: dataframe al:level of significance
  #create a list vL where element i carries attributes to vertex i
  R <- cor(d)       #Correlation matrix
  L <- dim(R)[1]    #Number of variabels
  N=dim(d)[1]       #Sampel size
  oneV <- function(x,R,N){#one element in vL
    L <- dim(R)[1]
    h <- R[x,-x]          #The row x from the correlation matrix without element x
    t <- rep(0,L)
    t[-x] <- abs(h*sqrt((N-2)/(1-h^2)))
    out <- rep(TRUE,L)    
    out[x] <- FALSE
    list(out=out,t=t,nei=NULL,w=NULL)
    #attributes to vertex i: out is edges not neighbors to i; t[j] is
    #t-score for removing (i,j); nei are the neighbors to i; w is
    #the inverse of R[nei,nei] needed to calculate t efficiently
  }
  vL <- lapply(1:L,oneV, R=R,N=N)  #cerating vertex list whith info about each vertex
  eli <- matrix(FALSE,L,L)#info about each edge. Is it eligable for
  #adding or removal?
  crit <- qnorm(1-al/2)
  for(i in 1:L) eli[i,vL[[i]]$t>crit] <- TRUE
  list(vL=vL,R=R,L=L,eli=eli,slut=FALSE,crit=crit,N=dim(d)[1])
  #model-objekt.
}
updT <- function(x,a,m){#vertex a has gathered or lost a connection 
  #Here x is the information about vertex a
  R <- m$R
  out <<- (1:m$L)[x$out]
  nei <<- x$nei
  dout <<- dim(out)
  dnei <<- dim(nei)
  r <- R[out,x$nei]
  w <<- x$w
  h <- R[out,x$nei]%*%x$w
  be <- x$w%*%R[x$nei,a]
  ssa <- 1-t.default(be)%*%R[x$nei,x$nei]%*%be
  dw <- dim(x$w)[1]
  vs <- x$w[1 + 0:(dw - 1) * (dw + 1)]
  sdev <- sqrt(ssa/(m$N-1-length(x$nei))*vs)
  if(sum(is.na(sdev))>0){
    assign("zzzm", m, envir = .GlobalEnv)
    stop("hell")
  }
  x$t[x$nei] <- abs(be/sdev)
  alA <- 1/(1-rowSums(h*R[out,x$nei]))
  Routa <- (h%*%R[x$nei,a])
  beA <- alA*(R[out,a]-Routa)[,1]
  ss <- 1-R[a,x$nei]%*%be-
    alA*((Routa[,1])^2+R[out,a]^2-2*(Routa[,1])*R[out,a])
  sdev <- sqrt(ss/(m$N-2-length(x$nei))*alA)
  if(sum(is.na(sdev))>0){
    assign("zzzm", m, envir = .GlobalEnv)
    stop("hell")
  }
  x$t[out] <- abs(beA/sdev)
  m$vL[[a]] <- x
  outok <- unlist(lapply(m$vL[x$out],function(x,a) x$t[a]>m$crit,a=a))
  m$eli[a,x$out] <- m$eli[x$out,a] <- (x$t[x$out]>m$crit) & outok
  neiOK <- unlist(lapply(m$vL[x$nei],function(x,a) x$t[a]<m$crit,a=a))
  m$eli[a,x$nei] <- m$eli[x$nei,a] <- (x$t[x$nei]<m$crit) & neiOK
  m
}
rmE <- function(a,b,x,m){#updates if we disconnect a from b
  i <- match(b,x$nei)
  dw <<- dim(x$w)[1]
  if(is.null(dw)) print(x)
  if(dw<3) {
    w <- if(dw==2) matrix(1) else NULL
  } else {
    w <- x$w[-i,-i]-x$w[-i,i,drop=FALSE]%*%x$w[i,-i]/x$w[i,i]
  }
  x$w <- w
  x$nei <- x$nei[-i]
  x$out[b] <- TRUE
  updT(x,a,m)
}
addE <- function(a,b,x,m){#updates if a is connected to b
  R <- m$R
  if(is.null(x$w)) w <- matrix(1) else {
    h <- x$w%*%R[x$nei,b]
    alAb <- 1/(1-sum(h*R[x$nei,b]))
    w <- x$w+alAb*h%*%t.default(h)
    h <- -alAb*h[,1]
    w <- rbind(w,h)
    w <- cbind(w,c(h,alAb))
  }
  x$w <- w
  x$nei <- c(x$nei,b)
  x$out[b] <- FALSE
  updT(x,a,m)
}

oneStep <- function(m){
  #"random" choise of "signif" edge
  #either eligable for add=TRUE or remove.
  eli <- (1:m$L^2)[m$eli]
  no <- if(length(eli)==1) eli else sample(eli,1)
  b <<- floor((no-1)/m$L)+1
  a <<- no-(b-1)*m$L
  add <<- match(b,m$vL[[a]]$nei,0)==0
  if(add) {
    m <- addE(a,b,m$vL[[a]],m)
    m <- addE(b,a,m$vL[[b]],m)
  } else {
    m <- rmE(a,b,m$vL[[a]],m)
    m <- rmE(b,a,m$vL[[b]],m)
  }
  m
}
step <- function(m,maxit=1000,al=.01){
  m$crit <- qnorm(1-al/2)
  #forward until noSignif or maxit edges
  for(i in 1:maxit){
    if(sum(m$eli)==0) break
    m <- oneStep(m)
  }
  m$it <- i
  m
}
edgeL <- function(m){
  res <- list()
  for(i in 1:m$L){
    x <- m$vL[[i]]$nei
    x <- sort(x)
    x <- x[x>i]
    if(length(x)>0) for(j in x)
      res <-
      c(res,list(colnames(d)[c(i,j)]))
  }
  res
}

set.seed(5)
m0 <- init(d)
Rprof(interval=.01)
set.seed(5)
m <- step(m0)
Rprof(NULL)
summaryRprof()
#586 iterations in 0.4 seconds
fitm <- cmod(edgeL(m),d);fitm
Rprof(interval=.01)
forwm <- forward(fitm,criterion="test",alpha=.01,type="unrestricted",search="headlong",steps=1)
Rprof(NULL)
summaryRprof()
#4.25 seconds for one iteration
testadd(fitm,c("x.3","x.6"))
#This cant be LRT. Something is wrong.

rm(zzzm)
