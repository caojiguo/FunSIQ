library(fda)
source("qsreg.R")


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
ntrain <- 500
ntest <- 100
N <- ntrain + ntest
score <- matrix(0, nrow=N, ncol=length(eigenv))
X <- matrix(0, nrow=N, ncol=m)

#tuning parameter
n_bs <- 4:7
lamb <- exp(seq(-25, -20, length.out=20))


#repeat 300 times
MSE <- numeric(300)
ISE <- numeric(300)

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
  
  GACV2 <- matrix(NA, nrow=length(lamb), ncol=length(n_bs))
  for(i in 1:nrow(GACV2))
  {
    for(j in 1:ncol(GACV2))
    {
      # cat("(i, j) = ", c(i, j), "\n")
      # res <- optim(par=seq(1,3,length.out=n_bs[j]-1), nllk.qsreg, num_bs=n_bs[j], X=X, timepts=t, Y=y0,
      #              tau=0.5, lambda=lamb[i], method="BFGS")
      # i = 1; j = 1
      
      res <- try(nlm(nllk.qsreg, p=seq(1,3,length.out=n_bs[j]-1),num_bs=n_bs[j],X=Xtrain,timepts=t, Y=Ytrain, tau=0.5, lambda=lamb[i]), silent=T)
      if(inherits(res, "try-error")) {next} else
      {
        GACV2[i,j] <- GACV.qsreg(X=Xtrain, timepts=t, Y=Ytrain, theta=res$estimate, num_bs=n_bs[j], tau=0.5, lambda=lamb[i])
      }
      
    }
  }
  
  colmin <- apply(GACV2, 2, min, na.rm=T)
  bs.opt <- n_bs[which.min(colmin)] # bs.opt <- 6
  lam.opt <- lamb[which.min(GACV2[,which.min(colmin)])]  #lam.opt <- 8.170972e-11
  
  #bs.opt <- 5; lam.opt <- lamb[20]
  res <- nlm(nllk.qsreg, p=seq(1,3,length.out=bs.opt-1),num_bs=bs.opt,X=Xtrain,timepts=t, Y=Ytrain, tau=0.5, lambda=lam.opt)
  
  bsbasis=create.bspline.basis(rangeval=c(0,1),nbasis=bs.opt, norder=4)
  basismat=eval.basis(t, bsbasis)
  v <- matrix(NA, nrow=N, ncol=bs.opt)
  for(i in 1:N)
    v[i, ] <- apply(X[i,]*basismat, 2, trapz, x=t)
  z <- v %*% c(1, res$estimate)
  
  fit <- qsreg(z[1:ntrain],Ytrain,alpha=0.5,lam=lam.opt)
  pred <-  predict(fit, z[(ntrain+1):N])
  
  MSE[k] <- mean((0.66*exp(utest^2) - pred)^2)
  
  betahat <- basismat %*% c(1, res$estimate)
  
  ISE[k] <- trapz(x=t, y=(betahat-beta(t))^2)
}

