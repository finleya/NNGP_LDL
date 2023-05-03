rm(list=ls())
dyn.load("nn.so")
source("nn.R")
source("util.R")
library(fields)
library(viridis)
set.seed(1)

n <- 100
coords <- cbind(runif(n,0,1), runif(n,0,1))
m <- 15

coords <- coords[order(coords[,1]),]

##assume and exp with params
phi <- 3/10
sigma.sq <- 1

##C and C inv
C <- sigma.sq*exp(-phi*rdist(coords))
C.inv.1 <- chol2inv(chol(C))

##C inv from (I-A)^T D^{-1}(I-A)
A <- matrix(0, n, n)
D <- rep(0,n)

for(i in 1:(n-1)){
    A[i+1,1:i] <- solve(sigma.sq*exp(-phi*rdist(coords[1:i,,drop=FALSE])), sigma.sq*exp(-phi*rdist(coords[1:i,,drop=FALSE], coords[i+1,,drop=FALSE])))
    D[i+1] <- sigma.sq - sigma.sq*exp(-phi*rdist(coords[i+1,,drop=FALSE], coords[1:i,,drop=FALSE])) %*% A[i+1,1:i]
}

D[1] <- sigma.sq

C.inv.2 <- t(diag(1,n) - A)%*%diag(1/D)%*%(diag(1,n) - A)

##C inv NN approx using (I-A)^T D^{-1}(I-A)
A <- matrix(0, n, n)
D <- rep(0,n)
nn <- mkNNIndx(coords, m)
nn.indx <- mk.n.indx.list(nn$nnIndx, n, m)

for(i in 1:(n-1)){
    A[i+1,nn.indx[[i+1]]] <- solve(sigma.sq*exp(-phi*rdist(coords[nn.indx[[i+1]],,drop=FALSE])), sigma.sq*exp(-phi*rdist(coords[nn.indx[[i+1]],,drop=FALSE], coords[i+1,,drop=FALSE])))
    D[i+1] <- sigma.sq - sigma.sq*exp(-phi*rdist(coords[i+1,,drop=FALSE], coords[nn.indx[[i+1]],,drop=FALSE])) %*% A[i+1,nn.indx[[i+1]]]
}

D[1] <- sigma.sq

C.inv.3 <- t(diag(1,n) - A)%*%diag(1/D)%*%(diag(1,n) - A)


plot.matrix <- function(m, tol=0.01, main=""){
    m[m <= tol] <- NA ##only for visualization
    image.plot(m[,nrow(m):1], col=viridis(10), main=main)
}

png(file="nn-test.png")
par(mfrow=c(1,3))
plot.matrix(C.inv.1, main="Direct chol")
plot.matrix(C.inv.2, main="Full (I-A)^T D^{-1}(I-A)")
plot.matrix(C.inv.3, main="NN Approx. (I-A)^T D^{-1}(I-A)")
dev.off()
