require(ftsa)
source("FSiM.R")  #functional single index model
source("qsreg.R") #functional quantile single index model


data("pm_10_GR_sqrt")

#preprocessing first
t <- seq(0,1,length.out=48)
daybasis <- create.fourier.basis(c(0, 1), nbasis=15)
pm10_smooth <- Data2fd(t, y=pm_10_GR_sqrt$y, basisobj=daybasis)
X=t(eval.fd(t, pm10_smooth))
Y <- apply(X, 1, max) #maximum of daily PM10 (square root)
n <- nrow(X)
#deinfe functional covariate and scalar response
fun.cov <- X[-n,]
scalar.res <- Y[-1] #the maximum PM10 of the next day as the scalar response

#hist(scalar.res); plot(density(scalar.res))
#boxplot(scalar.res)

#use all data to figure out what are good candidates for mean regression and quantile regression
#mean regression
n_Four_m <- c(11,13,15)
lamb_m <- 10^(seq(1, 3, length.out=20))
REML <- matrix(0, nrow=length(lamb_m), ncol=length(n_Four_m))
for(i in 1:length(lamb_m))
{ # i = 39
  for(j in 1:length(n_Four_m))
  { # j = 3
    cat("(i, j) = ", c(i, j), "\n")
   
    res <- nlm(nllk.fsim.F, p=seq(1,3,length.out=n_Four_m[j]-1),num_bs=n_Four_m[j],X=fun.cov,timepts=t,Y=scalar.res,lam =lamb_m[i])
    #GACV.qsreg <- function(X,timepts,Y,theta,num_bs,tau,lambda)
    #GACV.fsim(X=X, timepts=t, Y=y0, theta=res$estimate, num_bs=n_bs[j], lam=lamb[i])
    REML[i,j] <- GACV.fsim.F(X=fun.cov, timepts=t, Y=scalar.res, theta=res$estimate, num_bs=n_Four_m[j], lam=lamb_m[i])
    
  }
}
#windows()
#matplot(REML, type="l")

#quantile regression
n_Four_q <- c(11,13,15)
lamb_q <- 10^(seq(-2, 2, length.out=50))
GACV2 <- matrix(0, nrow=length(lamb_q), ncol=length(n_Four_q))
for(i in 1:length(lamb_q))
{
  for(j in 1:length(n_Four_q))
  {
    cat("(i, j) = ", c(i, j), "\n")
    res <- nlm(nllk.qsreg.F, p=seq(1,3,length.out=n_Four_q[j]-1),num_bs=n_Four_q[j],X=fun.cov,timepts=t,Y=scalar.res, tau=0.5,lambda=lamb_q[i])
    #GACV.qsreg <- function(X,timepts,Y,theta,num_bs,tau,lambda)
    GACV2[i,j] <- GACV.qsreg.F(X=fun.cov, timepts=t, Y=scalar.res,theta=res$estimate,num_bs=n_Four_q[j], tau=0.5,lambda=lamb_q[i])
    
  }
}
#windows()
#matplot(GACV2, type="l")

#=================================
#five fold cross-validation, 100 splits, compute prediction errors with two methods: mean regression vs quantile regression
K <- 100
id <- numeric(n-1)
CV_MSPE_q <- numeric(K)  #based on quantile regression
CV_MSPE_m <- numeric(K)  #based on mean regression
L <- floor((n-1)/5)
err_in_q <- err_in_m <- numeric(n-1)

n_Four_m <- 15
lamb_m <- 10^(seq(1, 3, length.out=10))
REML <- matrix(0, nrow=length(lamb_m), ncol=length(n_Four_m))

n_Four_q <- c(11,13,15)
lamb_q <- 10^(seq(-2, 2, length.out=30))
GACV2 <- matrix(100, nrow=length(lamb_q), ncol=length(n_Four_q))


#set.seed(123)
for(k in 11:20)
{ cat("k = ", k, "\n")
  set.seed(k)
  id[1:(n-2)] <- base::sample(rep(1:5, each=L),size=n-2)
  id[n-1] <- sample(1:5, size=1)
  for(s in 1:5)
  {
    train <- which(id!=s)
    test <- which(id==s)
    
   # mean regression
    for(i in 1:length(lamb_m))
    { 
      for(j in 1:length(n_Four_m))
      { 
        res <- nlm(nllk.fsim.F, p=seq(1,3,length.out=n_Four_m[j]-1),num_bs=n_Four_m[j],X=fun.cov[train,],timepts=t,Y=scalar.res[train],lam =lamb_m[i])
        REML[i,j] <- GACV.fsim.F(X=fun.cov[train,],timepts=t,Y=scalar.res[train],theta=res$estimate,num_bs=n_Four_m[j],lam=lamb_m[i])

      }
    }

    colmin <- apply(REML, 2, min)
    bs.opt <- n_Four_m[which.min(colmin)] 
    lam.opt <- lamb_m[which.min(REML[,which.min(colmin)])] 
    bs.opt <- n_Four_m; lam.opt <- lamb_m[which.min(REML)] 
    res <- nlm(nllk.fsim.F, p=seq(1,3,length.out=bs.opt-1),num_bs=bs.opt,X=fun.cov[train,],timepts=t,Y=scalar.res[train],lam=lam.opt)
    bsbasis=create.fourier.basis(rangeval=c(0,1),nbasis=bs.opt)
    basismat=eval.basis(t, bsbasis)
    v <- matrix(NA, nrow=n-1, ncol=bs.opt)
    for(i in 1:(n-1))
       v[i, ]<-apply(fun.cov[i,]*basismat,2,trapz, x=t)
    evensum <- sum(res$estimate[even(1:(bs.opt-1))])
    z <- v %*% c(1-evensum, res$estimate)
    u <- z[train]
    fit <- gam(scalar.res[train] ~ s(u, bs="cr"), method="REML", sp=lam.opt)
    pred <-  predict(fit, newdata= data.frame(u = z[test]))
    err_in_m[test] <- mean((scalar.res[test]-pred)^2)
     
    #===================================
    #quantile regression
    for(i in 1:length(lamb_q))
    {
      for(j in 1:length(n_Four_q))
      {
        cat("(i, j) = ", c(i, j), "\n")
        res <- nlm(nllk.qsreg.F, p=seq(1,3,length.out=n_Four_q[j]-1),num_bs=n_Four_q[j],X=fun.cov[train,],timepts=t,Y=scalar.res[train], tau=0.5,lambda=lamb_q[i])
        GACV2[i,j] <- GACV.qsreg.F(X=fun.cov[train,],timepts=t,Y=scalar.res[train],theta=res$estimate,num_bs=n_Four_q[j], tau=0.5,lambda=lamb_q[i])
        
      }
    }
    
    colmin <- apply(GACV2, 2, min)
    bs.opt <- n_Four_q[which.min(colmin)]
    lam.opt <- lamb_q[which.min(GACV2[,which.min(colmin)])]  
    res <- nlm(nllk.qsreg.F, p=seq(1,3,length.out=bs.opt-1),num_bs=bs.opt,X=fun.cov[train,],timepts=t,Y=scalar.res[train],tau=0.5,lambda=lam.opt)
    
    bsbasis=create.fourier.basis(rangeval=c(0,1),nbasis=bs.opt)
    basismat=eval.basis(t, bsbasis)
    v <- matrix(NA, nrow=(n-1), ncol=bs.opt)
    for(i in 1:(n-1))
      v[i, ] <- apply(fun.cov[i,]*basismat,2,trapz,x=t)
    evensum <- sum(res$estimate[even(1:(bs.opt-1))])
    z <- v %*% c(1-evensum, res$estimate)
    #z <- v %*% c(1, res$estimate)
    fit <- qsreg(z[train],scalar.res[train],alpha=0.5,lam=lam.opt)
    pred <-  predict(fit, z[test])
    err_in_q[test] <- mean((scalar.res[test]-pred)^2)
    
  }
  CV_MSPE_m[k]<-mean(err_in_m)
  CV_MSPE_q[k]<-mean(err_in_q)
}


#save(CV_MSPE_m, file="pm10_MSPE_mean.RData")
#save(CV_MSPE_q, file="pm10_MSPE_quantile.RData")

# MSPE <- cbind(CV_MSPE_m[1:20], CV_MSPE_q)
# colnames(MSPE) <- c("mean", "quantile")
# windows()
# op <- par(mar=c(5,7,4,2) + 0.1)
# par(op)
# #par()
# #pdf("MSPE.pdf")
# boxplot(MSPE, axes = FALSE, ann = FALSE) 
# axis(1, at = c(1,2), labels = c("mean", "quantile"), cex.axis = 1.2)
# axis(2, cex.axis = 2)
# title(ylab="Prediction Errors",cex.lab = 1.5,line=3)
# box()
# dev.off()

#var(Y)=3.59

#======================================================================
#use all data to estimate the index and link function for tau=0.5 and tau=0.75
# tau=0.5
for(i in 1:length(lamb_q))
{
  for(j in 1:length(n_Four_q))
  {
    #i = 13; j=3
    cat("(i, j) = ", c(i, j), "\n")
    res <- try(nlm(nllk.qsreg.F, p=seq(1,3,length.out=n_Four_q[j]-1),num_bs=n_Four_q[j],X=fun.cov,timepts=t,Y=scalar.res, tau=0.5,lambda=lamb_q[i]), silent=T)
    if(!(inherits(res, "try-error")))
    {
      GACV2[i,j] <- GACV.qsreg.F(X=fun.cov,timepts=t,Y=scalar.res,theta=res$estimate,num_bs=n_Four_q[j], tau=0.5,lambda=lamb_q[i])
    }
    
  }
}

#save(GACV2,file="PM10_50GACV.RData")
#load("PM10_50GACV.RData", verbose=T)

colmin <- apply(GACV2, 2, min)
bs.opt <- n_Four_q[which.min(colmin)]
lam.opt <- lamb_q[which.min(GACV2[,which.min(colmin)])]  
res <- nlm(nllk.qsreg.F, p=seq(1,3,length.out=bs.opt-1),num_bs=bs.opt,X=fun.cov,timepts=t,Y=scalar.res,tau=0.5,lambda=lam.opt)

bsbasis=create.fourier.basis(rangeval=c(0,1),nbasis=bs.opt)
basismat=eval.basis(t, bsbasis)
v <- matrix(NA, nrow=(n-1), ncol=bs.opt)
for(i in 1:(n-1))
  v[i, ] <- apply(fun.cov[i,]*basismat,2,trapz,x=t)
evensum <- sum(res$estimate[even(1:(bs.opt-1))])
z <- v %*% c(1-evensum, res$estimate)
fit <- qsreg(z,scalar.res,alpha=0.5,lam=lam.opt)

pdf("Residuals.pdf")
plot(fit$residuals, xlab="Day", ylab="Residuals")
lines(lowess(fit$residuals), col=2)
dev.off()


#pred <-  predict(fit, z[test])
fitted.fsimq <- predict(fit, sort(z))
u_medest <- sort(z); y_medest <- fitted.fsimq

#bootstrap procedure for link function f

B <- 500
boot.f50 <- matrix(NA, nrow=B, ncol=length(y_medest))
set.seed(123)
for(b in 1:500)
{
  cat("b = ", b, "\n")
  index.sel <- base::sample(1:(n-1), replace=T)
  newx <- fun.cov[index.sel,]; newy <- scalar.res[index.sel]
  # v <- matrix(NA, nrow=(n-1), ncol=bs.opt)
  # for(i in 1:(n-1))
  #   v[i, ] <- apply(newx[i,]*basismat,2,trapz,x=t)
  res <- try(nlm(nllk.qsreg.F, p=seq(1,3,length.out=bs.opt-1),num_bs=bs.opt,X=newx,timepts=t,Y=newy,tau=0.5,lambda=lam.opt), silent=T)
  if(!(inherits(res, "try-error")))
  {
    evensum <- sum(res$estimate[even(1:(bs.opt-1))])
    z <- v %*% c(1-evensum, res$estimate)
    fit <- qsreg(z,scalar.res,alpha=0.5,lam=lam.opt)
    boot.f50[b,] <- predict(fit, u_medest)
  }
}


#save(boot.beta50,file="PM10_bootbeta50.RData")
#load("PM10_bootbeta50.RData", verbose=T)
beta1 <- na.omit(boot.beta50[1:500,])
#beta50.band <- apply(beta1, 2, stats::quantile, probs=c(0.005, 0.995))
beta50.uppband <- beta50 + qnorm(0.975)*apply(beta1,2,sd)
beta50.lowband <- beta50 - qnorm(0.975)*apply(beta1,2,sd)



#windows()
pdf("PM10_link50.pdf")
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(x=z, y=scalar.res, xlab=expression(integral()*X(t)*beta(t)*d*t), ylab="f", main=expression(paste("link function, ", tau, "=0.5")))
lines(x=u_medest, y=y_medest)
dev.off()


##############################################
#bootstrap procedure for index function beta
beta50 <- basismat %*% c(1-evensum, res$estimate)
B <- 500
boot.beta50 <- matrix(NA, nrow=B, ncol=length(beta50))
for(b in 1:500)
{
  cat("b = ", b, "\n")
  index.sel <- base::sample(1:(n-1), replace=T)
  newx <- fun.cov[index.sel,]; newy <- scalar.res[index.sel]
  # v <- matrix(NA, nrow=(n-1), ncol=bs.opt)
  # for(i in 1:(n-1))
  #   v[i, ] <- apply(newx[i,]*basismat,2,trapz,x=t)
  res <- try(nlm(nllk.qsreg.F, p=seq(1,3,length.out=bs.opt-1),num_bs=bs.opt,X=newx,timepts=t,Y=newy,tau=0.5,lambda=lam.opt), silent=T)
  if(!(inherits(res, "try-error")))
  {
    evensum <- sum(res$estimate[even(1:(bs.opt-1))])
    boot.beta50[b,] <- basismat %*% c(1-evensum, res$estimate)
  }
}
#save(boot.beta50,file="PM10_bootbeta50.RData")
#load("PM10_bootbeta50.RData", verbose=T)
beta1 <- na.omit(boot.beta50[1:500,])
#beta50.band <- apply(beta1, 2, stats::quantile, probs=c(0.005, 0.995))
beta50.uppband <- beta50 + qnorm(0.975)*apply(beta1,2,sd)
beta50.lowband <- beta50 - qnorm(0.975)*apply(beta1,2,sd)
#windows()
pdf("PM10_index50.pdf")
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(x=t, y=beta50, type="l", xlab="t", ylab= expression(beta(t)), main = expression(paste("index function, ", tau, "=0.5")), ylim=c(min(beta50.lowband), max(beta50.uppband)))
lines(x=t, y=beta50.uppband,lty=2, col="red")
lines(x=t, y=beta50.lowband, lty=3, col="blue")
dev.off()



#####################################################################
## tau=0.75
for(i in 1:length(lamb_q))
{
  for(j in 1:length(n_Four_q))
  {
    cat("(i, j) = ", c(i, j), "\n")
    res <- nlm(nllk.qsreg.F, p=seq(1,3,length.out=n_Four_q[j]-1),num_bs=n_Four_q[j],X=fun.cov,timepts=t,Y=scalar.res,tau=0.75,lambda=lamb_q[i])
    GACV2[i,j] <- GACV.qsreg.F(X=fun.cov,timepts=t,Y=scalar.res,theta=res$estimate,num_bs=n_Four_q[j], tau=0.75,lambda=lamb_q[i])
    
  }
}

#save(GACV2,file="PM10_75GACV.RData")
#load("PM10_75GACV.RData", verbose=T)

colmin <- apply(GACV2, 2, min)
bs.opt <- n_Four_q[which.min(colmin)]
lam.opt <- lamb_q[which.min(GACV2[,which.min(colmin)])]  
res <- nlm(nllk.qsreg.F, p=seq(1,3,length.out=bs.opt-1),num_bs=bs.opt,X=fun.cov,timepts=t,Y=scalar.res,tau=0.75,lambda=lam.opt)

bsbasis=create.fourier.basis(rangeval=c(0,1),nbasis=bs.opt)
basismat=eval.basis(t, bsbasis)
v <- matrix(NA, nrow=(n-1), ncol=bs.opt)
for(i in 1:(n-1))
  v[i, ] <- apply(fun.cov[i,]*basismat,2,trapz,x=t)
evensum <- sum(res$estimate[even(1:(bs.opt-1))])
z <- v %*% c(1-evensum, res$estimate)
fit <- qsreg(z,scalar.res,alpha=0.75,lam=lam.opt)
#pred <-  predict(fit, z[test])
fitted.fsimq <- predict(fit, sort(z))
u_75est <- sort(z); y_75est <- fitted.fsimq

#bootstrap procedure for the link function f
B <- 500
boot.f75 <- matrix(NA, nrow=B, ncol=length(y_75est))
set.seed(123)
for(b in 1:500)
{
  cat("b = ", b, "\n")
  index.sel <- base::sample(1:(n-1), replace=T)
  newx <- fun.cov[index.sel,]; newy <- scalar.res[index.sel]
  # v <- matrix(NA, nrow=(n-1), ncol=bs.opt)
  # for(i in 1:(n-1))
  #   v[i, ] <- apply(newx[i,]*basismat,2,trapz,x=t)
  res <- try(nlm(nllk.qsreg.F, p=seq(1,3,length.out=bs.opt-1),num_bs=bs.opt,X=newx,timepts=t,Y=newy,tau=0.75,lambda=lam.opt), silent=T)
  if(!(inherits(res, "try-error")))
  {
    evensum <- sum(res$estimate[even(1:(bs.opt-1))])
    z <- v %*% c(1-evensum, res$estimate)
    fit <- qsreg(z,scalar.res,alpha=0.75,lam=lam.opt)
    boot.f75[b,] <- predict(fit, u_75est)
  }
}

boot.f <- na.omit(boot.f75[1:500,])
#beta50.band <- apply(beta1, 2, stats::quantile, probs=c(0.005, 0.995))
f75.uppband <- y_75est + qnorm(0.975)*apply(boot.f,2,sd)
f75.lowband <- y_75est - qnorm(0.975)*apply(boot.f,2,sd)

#windows()
pdf("PM10_link75.pdf")
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(x=z, y=scalar.res, xlab=expression(integral()*X(t)*beta(t)*d*t), ylab="f", main=expression(paste("link function, ", tau, "=0.75")))
lines(x=u_75est, y=y_75est)
lines(x=u_75est, y=f75.uppband,lty=2, col="red")
lines(x=u_75est, y=f75.lowband, lty=3, col="blue")
dev.off()




beta75 <- basismat %*% c(1-evensum, res$estimate)
B <- 500
boot.beta75 <- matrix(NA, nrow=B, ncol=length(beta75))
for(b in 1:500)
{
  cat("b = ", b, "\n")
  index.sel <- base::sample(1:(n-1), replace=T)
  newx <- fun.cov[index.sel,]; newy <- scalar.res[index.sel]
  # v <- matrix(NA, nrow=(n-1), ncol=bs.opt)
  # for(i in 1:(n-1))
  #   v[i, ] <- apply(newx[i,]*basismat,2,trapz,x=t)
  res <- try(nlm(nllk.qsreg.F, p=seq(1,3,length.out=bs.opt-1),num_bs=bs.opt,X=newx,timepts=t,Y=newy,tau=0.75,lambda=lam.opt), silent=T)
  if(!(inherits(res, "try-error")))
  {
    evensum <- sum(res$estimate[even(1:(bs.opt-1))])
    boot.beta75[b,] <- basismat %*% c(1-evensum, res$estimate)
  }
}


#save(boot.beta75,file="PM10_bootbeta75.RData")
#load("PM10_bootbeta75.RData", verbose=T)


beta1 <- na.omit(boot.beta75[1:500,])
#beta75.band <- apply(beta1, 2, stats::quantile, probs=c(0.005, 0.995))
beta75.uppband <- beta75 + qnorm(0.975)*apply(beta1,2,sd)
beta75.lowband <- beta75 - qnorm(0.975)*apply(beta1,2,sd)


pdf("PM10_index75.pdf")
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(x=t, y=beta75, type="l", xlab="t", ylab=expression(beta(t)), main = expression(paste("index function, ", tau, "=0.75")), ylim=c(min(beta75.lowband), max(beta75.uppband)))
lines(x=t, y=beta75.uppband,lty=2, col="red")
lines(x=t, y=beta75.lowband, lty=3, col="blue")
dev.off()


