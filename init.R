library(mlbench)
data(Satellite)
d <- Satellite[Satellite[,37]=="red soil",-37]

init <- function(d,al=.01){#d: dataframe al:level of significance
    #create a list vL where element i carries attributes to vertex i
    R <- cor(d)
    L <- dim(R)[1]
    oneV <- function(x,R){#one element in vL
        L <- dim(R)[1]
        h <- R[x,-x]
        t <- rep(0,L)
        t[-x] <- abs(h*sqrt((L-2)/(1-h^2)))
        out <- rep(TRUE,L)
        out[x] <- FALSE
        list(out=out,t=t,nei=NULL,w=NULL)
        #attributes to vertex i: out is edges not neighbors to i; t[j] is
        #t-score for removing (i,j); nei are the neighbors to i; w is
        #the inverse of R[nei,nei] needed to calculate t efficiently
    }
    vL <- lapply(1:L,oneV, R=R)
    list(vL=vL,R=R,L=L,slut=FALSE,crit=qnorm(1-al/2),N=dim(d)[1])
    #model-objekt.
}
m0 <- m <- init(d)

addE <- function(a,b,x,R){
    #When a is connected to b we need to update vL[[a]]
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
    out <- (1:m$L)[x$out]
    h <- R[out,x$nei]%*%x$w
    alA <- 1/(1-rowSums(h*R[out,x$nei]))
    Routa <- (h%*%R[x$nei,a])
    beA <- alA*(R[out,a]-Routa)[,1]
    ss <- 1-R[a,x$nei]%*%x$w%*%R[x$nei,a]-
        alA*((Routa[,1])^2+R[out,a]^2-2*(Routa[,1])*R[out,a])
    x$t[out] <- abs(beA/sqrt(ss/(m$N-1-length(x$nei))*alA))
    x
}

set.seed(5)
oneStep <- function(m){
    #"random" choise of signif edge
    ant <- m$L
    ins <- 1:m$L
    while(ant>0){
        i <- sample(ant,1)
        a <- ins[i]
        x <- m$vL[[a]]
        OK <- x$t>m$crit & x$out
        if(sum(OK)>0) break
        ant <- ant-1
        ins <- ins[-i]
    }
    if(ant==0){
        m$slut <- TRUE
        return(m)
    }
    h <- (1:m$L)[OK]
    b <- if(length(h)==1) h else sample(h,1)
    m$vL[[a]] <- addE(a,b,x,m$R)
    m$vL[[b]] <- addE(b,a,m$vL[[b]],m$R)
    m
}
step <- function(m,maxit=1000){
    #forward until noSignif or maxit edges
    for(i in 1:maxit){
        if(m$slut) break
        m <- oneStep(m)
    }
    m
}
Rprof(interval=.01)
m <- step(m0)
Rprof(NULL)
summaryRprof()
#351 edges added in 0.14 seconds

edgeL <- list()
for(i in 1:m$L){
    x <- m$vL[[i]]$nei
    x <- sort(x)
    x <- x[x>i]
    if(length(x)>0) for(j in x) edgeL <- c(edgeL,list(colnames(d)[c(i,j)]))
}


library(gRim)
fitm <- cmod(edgeL,d);fitm
Rprof(interval=.01)
forwm <- forward(fitm,criterion="test",alpha=.01,type="unrestricted",search="headlong",steps=1)
Rprof(NULL)
summaryRprof()
#Adding one edge: 1.88 seconds



#m1 <- cmod(~.^1,d) This is prohibitive:
#forwm <- forward(cm1,criterion="test",alpha=.01,type="unrestricted",search="headlong")
