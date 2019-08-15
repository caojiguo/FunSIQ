
library(fda)
library(caTools)
library(mgcv)

even <- function(x) x%%2 == 0 

#profile function for f_tau
#theta: coefficient vector for beta
#num_bs: number of B-spline basis functions
#Tmax: right endpoint of the functional observations
#X: observations for functional covariates, an n*m matrix
#timepts: observational time points for X


#lam: smoothing parameter, 
#Y: response variable


nllk.fsim <- function(X,timepts,Y,theta,num_bs,lam)
{
  #theta <- seq(1,3,length.out=n_bs[j]-1); num_bs=n_bs[j],X=X[1:500,],timepts=t, Y=resp$y0[1:500],
  # tau=0.5, lambda=lamb[i]
  theta_int <- c(1, theta[1:(num_bs-1)])
  #timepts=seq(0,10,length.out=172); num_bs=25; Tmax=max(timepts)
  Tmax <- max(timepts)
  bsbasis=create.bspline.basis(rangeval=c(0,Tmax),nbasis=num_bs, norder=4)
  basismat=eval.basis(timepts, bsbasis)
  #dim(basismat) = length(timepts) * num_bs
  
  # X <- rbind(1:172, 2:173)/100
  u <- matrix(NA, nrow=nrow(X), ncol=num_bs)
  for(i in 1:nrow(X))
    u[i, ] <- apply(X[i,]*basismat, 2, trapz, x=timepts)
  z <- u %*% theta_int
  
  #fit <- smooth.spline(z, Y, lambda=lam)
  fit <- gam(Y ~ s(z, bs="cr"), method="REML", sp=lam)
  # sc = sqrt(var(y)) * 1e-05
  # temp <- mean(qsreg.rho(fit$residuals, alpha = tau, C = sc))
  if(any(abs(theta_int) > 50)) return (1e08)
  temp <- mean((fit$residuals)^2)
  temp
}


#use Fourier basis function

nllk.fsim.F <- function(X,timepts,Y,theta,num_bs,lam)
{
  #theta <- seq(1,3,length.out=n_bs[j]-1); num_bs=n_bs[j],X=X[1:500,],timepts=t, Y=resp$y0[1:500],
  # tau=0.5, lambda=lamb[i]
  evensum <- sum(theta[even(1:(num_bs-1))])
  theta_int <- c(1-evensum, theta[1:(num_bs-1)])
  #timepts=seq(0,10,length.out=172); num_bs=25; Tmax=max(timepts)
  Tmax <- max(timepts)
  bsbasis=create.fourier.basis(rangeval=c(0,Tmax),nbasis=num_bs)
  basismat=eval.basis(timepts, bsbasis)
  #dim(basismat) = length(timepts) * num_bs
  
  # X <- rbind(1:172, 2:173)/100
  u <- matrix(NA, nrow=nrow(X), ncol=num_bs)
  for(i in 1:nrow(X))
    u[i, ] <- apply(X[i,]*basismat, 2, trapz, x=timepts)
  z <- u %*% theta_int
  
  #fit <- smooth.spline(z, Y, lambda=lam)
  fit <- gam(Y ~ s(z, bs="cr"), method="REML", sp=lam)
  # sc = sqrt(var(y)) * 1e-05
  # temp <- mean(qsreg.rho(fit$residuals, alpha = tau, C = sc))
  if(any(abs(theta_int) > 50)) return (1e08)
  temp <- mean((fit$residuals)^2)
  temp
}

#use GACV to choose tuning parameters: lam and num_bs
GACV.fsim <- function(X,timepts,Y,theta,num_bs,tau, lam)
{
  #theta <- seq(0, 5, length.out=50)
  theta_int <- c(1, theta[1:(num_bs-1)])
  #timepts=seq(0,10,length.out=172); num_bs=25; Tmax=max(timepts)
  Tmax <- max(timepts)
  bsbasis=create.bspline.basis(rangeval=c(0,Tmax),nbasis=num_bs, norder=4)
  basismat=eval.basis(timepts, bsbasis)
  #dim(basismat) = length(timepts) * num_bs
  n <- nrow(X)
  # X <- rbind(1:172, 2:173)/100
  u <- matrix(NA, nrow=n, ncol=num_bs)
  for(i in 1:nrow(X))
    u[i, ] <- apply(X[i,]*basismat, 2, trapz, x=timepts)
  z <- u %*% theta_int
  
  fit <- gam(Y ~ s(z, bs="cr"), method="REML", sp=lam)
  temp <- as.numeric(fit$gcv.ubre)
  temp
}



#use GACV to choose tuning parameters: lam and num_bs, Fourier basis
GACV.fsim.F <- function(X,timepts,Y,theta,num_bs,tau, lam)
{
  #theta <- seq(0, 5, length.out=50)
  evensum <- sum(theta[even(1:(num_bs-1))])
  theta_int <- c(1-evensum, theta[1:(num_bs-1)])
  #timepts=seq(0,10,length.out=172); num_bs=25; Tmax=max(timepts)
  Tmax <- max(timepts)
  bsbasis=create.fourier.basis(rangeval=c(0,Tmax),nbasis=num_bs)
  basismat=eval.basis(timepts, bsbasis)
  #dim(basismat) = length(timepts) * num_bs
  n <- nrow(X)
  # X <- rbind(1:172, 2:173)/100
  u <- matrix(NA, nrow=n, ncol=num_bs)
  for(i in 1:nrow(X))
    u[i, ] <- apply(X[i,]*basismat, 2, trapz, x=timepts)
  z <- u %*% theta_int
  
  fit <- gam(Y ~ s(z, bs="cr"), method="REML", sp=lam)
  temp <- as.numeric(fit$gcv.ubre)
  temp
}

