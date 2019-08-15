library(fda)
source("nonpafda.R")


t <- seq(0,1,length.out=100)
beta <- function(x) 2*x^2 + 0.25*x + 1  
m <- length(t)
mean.fun <- t 
#eigenvalue of covariance function
eigenv <- 0.5^(0:3)
#eigenfunction, FPC
eigenf1 <- function(x) sin(2*pi*x)/sqrt(2)
eigenf2 <- function(x) cos(2*pi*x)/sqrt(2)
eigenf3 <- function(x) sin(4*pi*x)/sqrt(2)
eigenf4 <- function(x) cos(4*pi*x)/sqrt(2)

#slope function
beta <- function(x) 2*x^2 + 0.25*x + 1  


#beta(0); beta(1); beta(0.25); beta(c(0,0.25))
a <- c()
a[1] <- integrate(function(x) x*beta(x), 0, 1)$value
a[2] <- integrate(function(x) eigenf1(x)*beta(x), 0, 1)$value
a[3] <- integrate(function(x) eigenf2(x)*beta(x), 0, 1)$value
a[4] <- integrate(function(x) eigenf3(x)*beta(x), 0, 1)$value
a[5] <- integrate(function(x) eigenf4(x)*beta(x), 0, 1)$value

B <- cbind(eigenf1(t), eigenf2(t), eigenf3(t), eigenf4(t))

#consider ntrain = 500
ntrain <- 300
ntest <- 100
N <- ntrain + ntest
score <- matrix(0, nrow=N, ncol=length(eigenv))
X <- matrix(0, nrow=N, ncol=m)

#repeat 300 times
MSE <- numeric(300)
#ISE <- numeric(300)

for(k in 1:300)
{ # k = 1
  
  set.seed(k)
  cat("k = ", k, "\n")
  #k = 1
  # index1 <- (N*(k-1)+1):(N*(k-1)+ntrain)
  # index2 <- (N*(k-1)+ntrain+1):(N*k)
  # Xtrain <- X[index1,]
  #Xtest <- X[index2]
  for(j in 1:4)
    score[,j] <- rnorm(N,sd=sqrt(eigenv[j]))
  X <- t(matrix(rep(mean.fun,N),nrow=m,ncol=N))  + score %*% t(B)
  Xtrain <- X[1:ntrain,]
  Xtest <- X[(ntrain+1):N,]
  
  u <- as.vector(a[1] + score %*% a[-1])
  y0 <- 0.66*exp(u^2) + 0.85*rt(N,df=3)
  
  Ytrain <- y0[1:ntrain]
  utest <- u[(ntrain+1):N]
  
  result.reg <- funopare.knn.lcv(Ytrain, Xtrain, Xtest, 2, nknot=20, c(0,1))
  
  pred <-  result.reg$Predicted.values
  
  MSE[k] <- mean((0.66*exp(utest^2) - pred)^2)
  #MSE[k] <- mean((y0[(ntrain+1):N] - pred)^2)
  
}

boxplot(MSE)
MSE<-MSE[MSE<4]
mean(MSE)
sd(MSE)

#=========================================
#try chi-square distribution
library(fda)
source("nonpafda.R")


t <- seq(0,1,length.out=100)
beta <- function(x) 2*x^2 + 0.25*x + 1  
m <- length(t)
mean.fun <- t 
#eigenvalue of covariance function
eigenv <- 0.5^(0:3)
#eigenfunction, FPC
eigenf1 <- function(x) sin(2*pi*x)/sqrt(2)
eigenf2 <- function(x) cos(2*pi*x)/sqrt(2)
eigenf3 <- function(x) sin(4*pi*x)/sqrt(2)
eigenf4 <- function(x) cos(4*pi*x)/sqrt(2)

#slope function
beta <- function(x) 2*x^2 + 0.25*x + 1  


#beta(0); beta(1); beta(0.25); beta(c(0,0.25))
a <- c()
a[1] <- integrate(function(x) x*beta(x), 0, 1)$value
a[2] <- integrate(function(x) eigenf1(x)*beta(x), 0, 1)$value
a[3] <- integrate(function(x) eigenf2(x)*beta(x), 0, 1)$value
a[4] <- integrate(function(x) eigenf3(x)*beta(x), 0, 1)$value
a[5] <- integrate(function(x) eigenf4(x)*beta(x), 0, 1)$value

B <- cbind(eigenf1(t), eigenf2(t), eigenf3(t), eigenf4(t))

#consider ntrain = 500
ntrain <- 300
ntest <- 100
N <- ntrain + ntest
score <- matrix(0, nrow=N, ncol=length(eigenv))
X <- matrix(0, nrow=N, ncol=m)

#repeat 300 times
MSE <- numeric(300)
#ISE <- numeric(300)

for(k in 1:300)
{ # k = 1
  
  set.seed(k)
  cat("k = ", k, "\n")
  #k = 1
  # index1 <- (N*(k-1)+1):(N*(k-1)+ntrain)
  # index2 <- (N*(k-1)+ntrain+1):(N*k)
  # Xtrain <- X[index1,]
  #Xtest <- X[index2]
  for(j in 1:4)
    score[,j] <- rnorm(N,sd=sqrt(eigenv[j]))
  X <- t(matrix(rep(mean.fun,N),nrow=m,ncol=N))  + score %*% t(B)
  Xtrain <- X[1:ntrain,]
  Xtest <- X[(ntrain+1):N,]
  
  u <- as.vector(a[1] + score %*% a[-1])
  y0 <- 0.66*exp(u^2) + sqrt(2)*rnorm(N)
  
  Ytrain <- y0[1:ntrain]
  utest <- u[(ntrain+1):N]
  
  result.reg <- funopare.knn.lcv(Ytrain, Xtrain, Xtest, 2, nknot=20, c(0,1))
  
  pred <-  result.reg$Predicted.values
  
  MSE[k] <- mean((0.66*exp(utest^2) - pred)^2)
  #MSE[k] <- mean((y0[(ntrain+1):N] - pred)^2)
  
}

boxplot(MSE)
MSE<-MSE[MSE<5]
mean(MSE)
sd(MSE)


#==================================
#try exponential case

library(fda)
source("nonpafda.R")


t <- seq(0,1,length.out=100)
beta <- function(x) 2*x^2 + 0.25*x + 1  
m <- length(t)
mean.fun <- t 
#eigenvalue of covariance function
eigenv <- 0.5^(0:3)
#eigenfunction, FPC
eigenf1 <- function(x) sin(2*pi*x)/sqrt(2)
eigenf2 <- function(x) cos(2*pi*x)/sqrt(2)
eigenf3 <- function(x) sin(4*pi*x)/sqrt(2)
eigenf4 <- function(x) cos(4*pi*x)/sqrt(2)

#slope function
beta <- function(x) 2*x^2 + 0.25*x + 1  


#beta(0); beta(1); beta(0.25); beta(c(0,0.25))
a <- c()
a[1] <- integrate(function(x) x*beta(x), 0, 1)$value
a[2] <- integrate(function(x) eigenf1(x)*beta(x), 0, 1)$value
a[3] <- integrate(function(x) eigenf2(x)*beta(x), 0, 1)$value
a[4] <- integrate(function(x) eigenf3(x)*beta(x), 0, 1)$value
a[5] <- integrate(function(x) eigenf4(x)*beta(x), 0, 1)$value

B <- cbind(eigenf1(t), eigenf2(t), eigenf3(t), eigenf4(t))

#consider ntrain = 500
ntrain <- 100
ntest <- 100
N <- ntrain + ntest
score <- matrix(0, nrow=N, ncol=length(eigenv))
X <- matrix(0, nrow=N, ncol=m)

#repeat 300 times
MSE <- numeric(300)
#ISE <- numeric(300)

for(k in 1:300)
{ # k = 1
  
  set.seed(k)
  cat("k = ", k, "\n")
  #k = 1
  # index1 <- (N*(k-1)+1):(N*(k-1)+ntrain)
  # index2 <- (N*(k-1)+ntrain+1):(N*k)
  # Xtrain <- X[index1,]
  #Xtest <- X[index2]
  for(j in 1:4)
    score[,j] <- rexp(N,rate=1/sqrt(eigenv[j])) - sqrt(eigenv[j])
  X <- t(matrix(rep(mean.fun,N),nrow=m,ncol=N))  + score %*% t(B)
  Xtrain <- X[1:ntrain,]
  Xtest <- X[(ntrain+1):N,]
  
  u <- as.vector(a[1] + score %*% a[-1])
  y0 <- 0.66*exp(u^2) + 0.58*(rchisq(N, df=3) - 3)
  
  Ytrain <- y0[1:ntrain]
  utest <- u[(ntrain+1):N]
  
  result.reg <- funopare.knn.lcv(Ytrain, Xtrain, Xtest, 2, nknot=20, c(0,1))
  
  pred <-  result.reg$Predicted.values
  
  MSE[k] <- mean((0.66*exp(utest^2) - pred)^2)
  #MSE[k] <- mean((y0[(ntrain+1):N] - pred)^2)
  
}

boxplot(MSE)
MSE<-MSE[MSE<2]
mean(MSE)
sd(MSE)