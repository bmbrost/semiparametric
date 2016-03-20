rm(list=ls())

library(splines)
library(lmer)
library(nlme)


T <- 100  # number of observations

# Define covariates
time <- c(0,cumsum(rgamma(T-1,shape=1.1,scale=4.5)))  # time covariate
hr <- ((time)-24*floor(time/24))
day <- ceiling(time/24)
X <- cbind(1,day,hr)
X[,-1] <- scale(X[,-1])
qX <- ncol(X)

beta <- c(0.25,1.5,0.5)  # Coefficients on X

# Define non-linear trend to model non-parametrically 
trend <- 0.5*sin(0.02*time)  # non-linear pattern
plot(time,trend,type="l")

# Simulate data
sigma <- 1
y <- rnorm(T,X%*%beta+trend,sigma)
plot(time,y,type="l",ylim=range(c(trend,y,X%*%beta)))
lines(time,X%*%beta,col=2)
lines(time,trend,col=3)


#########################################################
### Maximum likelihood estimation with lme
#########################################################

###
### Truncated lines basis expansion
###

K <- 50  # number of knots
knots <- seq(min(time),max(time),length.out=K+2)[-c(1,K+2)]

# Basis expansion
Z <- outer(time,knots,"-")
Z <- Z*(Z>0)
matplot(Z,type="l",lty=1)

Id <- factor(rep(1,T))
fit <- lme(y~1,random=list(Id=pdIdent(~Z-1)))
y.hat <- fit$fitted[,2]
plot(time,y,type="l")
lines(time,y.hat,type="l",col=6)

###
### B-splines basis expansion
###

int <- 10  # interval between knots
knots <- seq(0,max(time),by=int)

# Basis exapnsion
Z <- bs(time,knots=knots,degree=3,intercept=FALSE)  # cubic spline
matplot(time,Z,type="l",lty=1,col=2,add=FALSE)
points(time,rep(-0.025,T),pch="|")

time.tmp <- seq(0,max(time),length.out=1000)
Z.tmp <- bs(time.tmp,knots=knots,degree=3,intercept=FALSE)
matplot(time.tmp,Z.tmp,type="l",lty=1,col="darkgrey",add=TRUE)
abline(v=knots,lty=2,col="lightgrey")

mod <- lm(y~Z)
plot(time,y,type="l")
lines(time,predict(mod),col=2)

Id <- factor(rep(1,length(y)))
fit <- lme(y~1,random=list(Id=pdIdent(~Z-1)))
y.hat <- fit$fitted[,2]
plot(time,y,type="l")
lines(time,y.hat,type="l",col=4)


#########################################################
### Estimation using MCMC
#########################################################

###
### Confirm algorithm recovers parameters used for simulation
###

# beta <- c(0.5,0.25,0.15)
beta <- c(0.25,1.5,1.0)

# Basis expansion
int <- 50  # interval between knots
knots <- seq(0,max(time),by=int)
Z <- bs(time,knots=knots,degree=3,intercept=FALSE)  # cubic spline
qZ <- ncol(Z)

sigma.alpha <- 1
alpha <- rnorm(qZ,0,sigma.alpha)

sigma <- 1
y <- rnorm(T,X%*%beta+Z%*%alpha,sigma)

plot(X%*%beta,type="l")
plot(Z%*%alpha,type="l")
plot(time,y,type="l")

Id <- factor(rep(1,T))
fit <- lme(y~day+hr,random=list(Id=pdIdent(~Z-1)))
y.hat.lme <- fit$fitted[,2]
lines(time,y.hat.lme,col=2)
beta.hat.lme <- as.matrix(fixef(fit))
alpha.hat.lme <- t(as.matrix(ranef(fit)))

source('~/Documents/git/SemiPar/normal.semipar.mcmc.R', chdir = TRUE)
start <- list(beta=beta,alpha=alpha,sigma=sigma,sigma.alpha=sigma.alpha)
# hist(sqrt(1/rgamma(1000,1,,2)))
priors <- list(sigma.beta=10,r.sigma=2,q.sigma=1)
out1 <- normal.semipar.mcmc(y,X,Z,priors=priors,start=start,sigma.alpha=NULL,n.mcmc=50000)
out1$DIC

matplot(out1$beta,type="l",lty=1);abline(h=beta,col=1:3)

beta.hat <- apply(out1$beta,2,mean)
beta.quant <- t(apply(out1$beta,2,quantile,c(0.025,0.975)))
plot(beta.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(beta.quant)))
abline(h=0,col=2,lty=2)
segments(1:qX,beta.quant[,1],1:qX,beta.quant[,2],col="lightgrey")
points(beta.hat,pch=19,col=rgb(0,0,0,0.25))
points(beta,pch=19)
points(beta.hat.lme,pch=19,col=3)

idx <- 3
matplot(out1$alpha[,idx],type="l",lty=1);abline(h=alpha[idx],col=1:3)

alpha.hat <- apply(out1$alpha,2,mean)
alpha.quant <- t(apply(out1$alpha,2,quantile,c(0.025,0.975)))
sum(alpha>alpha.quant[,1]&alpha<alpha.quant[,2])/qZ  # coverage probability
plot(alpha.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(alpha.quant)))
abline(h=0,col=2,lty=2)
segments(1:qZ,alpha.quant[,1],1:qZ,alpha.quant[,2],col="lightgrey")
points(alpha.hat,pch=19,col=rgb(0,0,0,0.25))
points(alpha,pch=19)
points(alpha.hat.lme,pch=19,col=3)

matplot(out1$sigma,type="l");abline(h=sigma,col=2,lty=2)
matplot(out1$sigma.alpha,type="l");abline(h=sigma.alpha)

y.hat <- apply(out1$y.hat,1,mean)
plot(time,y,type="l")
lines(time,y.hat,col=3)
lines(time,y.hat.lme,col=2)


###
### Fit model with 'unknown' non-linear trend
###

beta <- c(0.25,1.5,1.0)

sigma <- 1
y <- rnorm(T,X%*%beta+trend,sigma)
plot(time,y,type="l",ylim=range(c(trend,y,X%*%beta)))
lines(time,X%*%beta,col=2)
lines(time,trend,col=3)

# Basis expansion
int <- 20  # interval between knots
knots <- seq(0,max(time),by=int)
Z <- bs(time,knots=knots,degree=3,intercept=FALSE)  # cubic spline
qZ <- ncol(Z)

Id <- factor(rep(1,T))
fit <- lme(y~day+hr,random=list(Id=pdIdent(~Z-1)))
y.hat.lme <- fit$fitted[,2]
lines(time,y.hat.lme,col=2)
beta.hat.lme <- as.matrix(fixef(fit))
alpha.hat.lme <- t(as.matrix(ranef(fit)))

source('~/Documents/git/SemiPar/normal.semipar.mcmc.R', chdir = TRUE)
start <- list(beta=beta,sigma=sigma,sigma.alpha=1)
# hist(sqrt(1/rgamma(1000,1,,2)))
priors <- list(sigma.beta=10,r.sigma=2,q.sigma=1)
out1 <- normal.semipar.mcmc(y,X,Z,priors=priors,start=start,sigma.alpha=NULL,n.mcmc=10000)
out1$DIC

matplot(out1$beta,type="l",lty=1);abline(h=beta,col=1:3)

beta.hat <- apply(out1$beta,2,mean)
beta.quant <- t(apply(out1$beta,2,quantile,c(0.025,0.975)))
plot(beta.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(beta.quant)))
abline(h=0,col=2,lty=2)
segments(1:qX,beta.quant[,1],1:qX,beta.quant[,2],col="lightgrey")
points(beta.hat,pch=19,col=rgb(0,0,0,0.25))
points(beta,pch=19)
points(beta.hat.lme,pch=19,col=3)

idx <- 3
matplot(out1$alpha[,idx],type="l",lty=1)

alpha.hat <- apply(out1$alpha,2,mean)
alpha.quant <- t(apply(out1$alpha,2,quantile,c(0.025,0.975)))
plot(alpha.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(alpha.quant)))
abline(h=0,col=2,lty=2)
segments(1:qZ,alpha.quant[,1],1:qZ,alpha.quant[,2],col="lightgrey")
points(alpha.hat,pch=19,col=rgb(0,0,0,0.25))
points(alpha.hat.lme,pch=19,col=3)

matplot(out1$sigma,type="l");abline(h=sigma,col=2,lty=2)
matplot(out1$sigma.alpha,type="l");abline(h=sigma.alpha)

y.hat <- apply(out1$y.hat,1,mean)
plot(time,y,type="l")
lines(time,y.hat,col=3)
lines(time,y.hat.lme,col=2)

par(mfrow=c(2,1))
plot(time,trend,type="l")
plot(alpha.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(alpha.quant)))
lines(alpha.hat,pch=19,col=rgb(0,0,0,0.35))
abline(h=0,col=2,lty=2)

#########################################################
### Regularization (selection of sigma.alpha)
#########################################################

coarse.grid <- seq(0.01,2,0.1)  # coarse grid for sigma.alpha
l.coarse <- length(coarse.grid)
DIC.coarse <- numeric(l.coarse)
for(i in 1:l.coarse){
	DIC.coarse[i] <- normal.semipar.mcmc(y,X,Z,priors=priors,start=start,
		sigma.alpha=coarse.grid[i],n.mcmc=1000)$DIC
}

plot(coarse.grid,DIC.coarse,type="l")
idx <- which.min(DIC.coarse)

l.fine <- 10
fine.grid <- seq(coarse.grid[idx-1],coarse.grid[idx+1],length.out=10)
DIC.fine <- numeric(l.fine)
for(i in 1:l.fine){
	DIC.fine[i] <- normal.semipar.mcmc(y,X,Z,priors=priors,start=start,
		sigma.alpha=fine.grid[i],n.mcmc=1000)$DIC
}

plot(fine.grid,DIC.fine,type="l")
fine.grid[which.min(DIC.fine)]