#define the check function for quantile regression
#alpha: quantile

check.rho <- function(r, alpha) 
{
  temp <- ifelse(r < 0, alpha*r - r, alpha*r)
  temp
}

# check.rho(0.5,0.5)
# check.rho(-0.5,0.5)
# check.rho(0,0.2)
# check.rho(0.5,0.1)


#composite simpson rule for integral
# int.simp = function(x,y)
# {
#   a <- min(x)
#   b <- max(x)
#   n <- length(x)-1
#   h = (b-a)/n
#   xi0 = y[1]+y[n+1]
#   xi1 = 0
#   xi2 = 0
#   for (i in 1:n-1)
#   {
#     #x = a+i*h
#     if (i%%2==0)
#     {
#       xi2 = xi2+y[i+1]
#     } else{
#       xi1 = xi1+y[i+1];
#      }
#   }
#   xi = h*(xi0+2*xi2+4*xi1)/3
#   xi
# } 
# 
# tt <- seq(0,1,length.out=200)
# 
# 
# int.simp(x=tt, y=sin(tt))
#   
# trapz(x=tt,y=sin(tt))
#     
  
even <- function(x) x%%2 == 0 




library(fda)
library(caTools)
library(fields)


#profile function for f_tau
#theta: coefficient vector for beta
#num_bs: number of B-spline basis functions
#Tmax: right endpoint of the functional observations
#X: observations for functional covariates, an n*m matrix
#timepts: observational time points for X

#tau: quantile level;
#lambda: smoothing parameter 
#Y: response variable


nllk.qsreg <- function(X,timepts,Y,theta,num_bs,tau,lambda)
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
  
  fit <- qsreg(z,Y,alpha=tau,lam=lambda)
 # sc = sqrt(var(y)) * 1e-05
 # temp <- mean(qsreg.rho(fit$residuals, alpha = tau, C = sc))
  if(any(abs(theta_int) > 50)) return (1e08)
  temp <- mean(check.rho(fit$residuals, alpha=tau))
  temp
}

#use Fourier basis functions
nllk.qsreg.F <- function(X,timepts,Y,theta,num_bs,tau,lambda)
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
  
  fit <- qsreg(z,Y,alpha=tau,lam=lambda)
 # sc = sqrt(var(y)) * 1e-05
 # temp <- mean(qsreg.rho(fit$residuals, alpha = tau, C = sc))
  if(any(abs(theta_int) > 50)) return (1e08)
  temp <- mean(check.rho(fit$residuals, alpha=tau))
  temp
}


#use GACV to choose tuning parameters: lam and num_bs
GACV.qsreg <- function(X,timepts,Y,theta,num_bs,tau,lambda)
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
  
  fit <- qsreg(z,Y,alpha=tau,lam= lambda)
  temp <- (sum(check.rho(fit$residuals, alpha=tau)))/(n-fit$trace)
  temp
}


#use Fourier basis functions
GACV.qsreg.F <- function(X,timepts,Y,theta,num_bs,tau,lambda)
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
  
  fit <- qsreg(z,Y,alpha=tau,lam= lambda)
  temp <- (sum(check.rho(fit$residuals, alpha=tau)))/(n-fit$trace)
  temp
}
