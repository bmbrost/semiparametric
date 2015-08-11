rm(list=ls())

library(splines)
# library(lmer)
# library(nlme)
# library(lme4)

T <- 100  # number of observations

# Define covariates
time <- c(0,cumsum(rgamma(T-1,shape=1.1,scale=4.5)))  # time covariate
hr <- ((time)-24*floor(time/24))
day <- ceiling(time/24)
X <- cbind(1,day,hr)
X[,-1] <- scale(X[,-1])
qX <- ncol(X)

# beta <- c(-0.1,1.25,0.5)  # Coefficients on X
beta <- c(-0.5,1.25,0.5)  # Coefficients on X


###
### Confirm MCMC algorithm recovers alphas used for simulation
###

# Basis expansion
int <- 10  # interval between knots
knots <- seq(0,max(time),by=int)
Z <- bs(time,knots=knots,degree=3,intercept=FALSE)  # cubic spline
qZ <- ncol(Z)

sigma.alpha <- 5
alpha <- rnorm(qZ,0,sigma.alpha)
trend <- Z%*%alpha

p <- pnorm(X%*%beta+trend)  # probability of being hauled-out
hist(p);summary(p)
y <- rbinom(T,1,p)  # haulout indicator variable: 1=hauled-out, 0=at-sea
table(y)

plot(time,y,ylim=range(c(trend,y,X%*%beta)))
lines(time,X%*%beta,col=2)
lines(time,trend,col=3)
lines(time,X%*%beta+trend,col=4)

source('~/Documents/git/SemiparametricRegression/probit.semipar.mcmc.R', chdir = TRUE)
start <- list(beta=beta,alpha=alpha)
# hist(sqrt(1/rgamma(1000,1,,2)))
priors <- list(mu.beta=rep(0,qX),sigma.beta=10)
out1 <- probit.semipar.mcmc(y,X,Z,priors=priors,start=start,sigma.alpha=sigma.alpha,
	n.mcmc=10000)
out1$DIC

matplot(out1$beta,type="l",lty=1);abline(h=beta,col=1:3,lty=2)

beta.hat <- apply(out1$beta,2,mean)
beta.quant <- t(apply(out1$beta,2,quantile,c(0.025,0.975)))
plot(beta.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(beta.quant)))
abline(h=0,col=2,lty=2)
segments(1:qX,beta.quant[,1],1:qX,beta.quant[,2],col="lightgrey")
points(beta.hat,pch=19,col=rgb(0,0,0,0.25))
points(beta,pch=19)

alpha.hat <- apply(out1$alpha,2,mean)
alpha.quant <- t(apply(out1$alpha,2,quantile,c(0.025,0.975)))
plot(alpha.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(alpha.quant)))
abline(h=0,col=2,lty=2)
segments(1:qZ,alpha.quant[,1],1:qZ,alpha.quant[,2],col="lightgrey")
points(alpha.hat,pch=19,col=rgb(0,0,0,0.25))
points(alpha,pch=19,col=3)

boxplot(pnorm(out1$u),col=8,outline=FALSE)
points(y,col=3,pch=19,cex=0.5)

u.inv <- matrix(pnorm(out1$u),,T)
u.inv.mean <- apply(u.inv,2,mean)
u.inv.quant <- t(apply(u.inv,2,quantile,c(0.025,0.975)))
plot(u.inv.mean,pch=19,col=rgb(0,0,0,0.25),ylim=c(0,1))
segments(1:T,u.inv.quant[,1],1:T,u.inv.quant[,2],col=rgb(0,0,0,0.15))
points(y,col=3,pch=19,cex=0.5)

par(mfrow=c(2,1))
plot(time,trend,type="l")
plot(alpha.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(alpha.quant)))
lines(alpha.hat,col=rgb(0,0,0,0.25))
abline(h=0,col=2,lty=2)


###
### Fit model with 'unknown' non-linear trend
###

# Define non-linear trend to model non-parametrically 
trend <- 0.5*sin(0.1*time)  # non-linear pattern
trend <- 2*sin(0.1*time)  # non-linear pattern
plot(time,trend,type="l")

# Simulate data
p <- pnorm(X%*%beta+trend)  # probability of being hauled-out
hist(p);summary(p)
y <- rbinom(T,1,p)  # haulout indicator variable: 1=hauled-out, 0=at-sea
table(y)

plot(time,y,ylim=range(c(trend,y,X%*%beta)))
lines(time,X%*%beta,col=2)
lines(time,trend,col=3)
lines(time,X%*%beta+trend,col=4)

# Basis expansion
int <- 10  # interval between knots
knots <- seq(0,max(time),by=int)
Z <- bs(time,knots=knots,degree=3,intercept=FALSE)  # cubic spline
matplot(Z,type="l")
qZ <- ncol(Z)

source('~/Documents/git/SemiparametricRegression/probit.semipar.mcmc.R', chdir = TRUE)
start <- list(beta=beta,alpha=rep(0,qZ))
# hist(sqrt(1/rgamma(1000,1,,2)))
priors <- list(mu.beta=rep(0,qX),sigma.beta=10)
out1 <- probit.semipar.mcmc(y,X,Z,priors=priors,start=start,sigma.alpha=1,n.mcmc=1000)
out1$DIC

matplot(out1$beta,type="l",lty=1);abline(h=beta,col=1:3,lty=2)

beta.hat <- apply(out1$beta,2,mean)
beta.quant <- t(apply(out1$beta,2,quantile,c(0.025,0.975)))
plot(beta.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(beta.quant)))
abline(h=0,col=2,lty=2)
segments(1:qX,beta.quant[,1],1:qX,beta.quant[,2],col="lightgrey")
points(beta.hat,pch=19,col=rgb(0,0,0,0.25))
points(beta,pch=19)

alpha.hat <- apply(out1$alpha,2,mean)
alpha.quant <- t(apply(out1$alpha,2,quantile,c(0.025,0.975)))
plot(alpha.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(alpha.quant)))
abline(h=0,col=2,lty=2)
segments(1:qZ,alpha.quant[,1],1:qZ,alpha.quant[,2],col="lightgrey")
points(alpha.hat,pch=19,col=rgb(0,0,0,0.25))

boxplot(pnorm(out1$u),col=8,outline=FALSE)
points(y,col=3,pch=19,cex=0.5)

u.inv <- matrix(pnorm(out1$u),,T)
u.inv.mean <- apply(u.inv,2,mean)
u.inv.quant <- t(apply(u.inv,2,quantile,c(0.025,0.975)))
plot(u.inv.mean,pch=19,col=rgb(0,0,0,0.25),ylim=c(0,1))
segments(1:T,u.inv.quant[,1],1:T,u.inv.quant[,2],col=rgb(0,0,0,0.15))
points(y,col=3,pch=19,cex=0.5)

par(mfrow=c(2,1))
plot(time,trend,type="l")
plot(alpha.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(alpha.quant)))
lines(alpha.hat,col=rgb(0,0,0,0.25))
abline(h=0,col=2,lty=2)


###
### Regularization (selection of sigma.alpha)
###

coarse.grid <- seq(0.01,2,0.1)  # coarse grid for sigma.alpha
l.coarse <- length(coarse.grid)
DIC.coarse <- numeric(l.coarse)
for(i in 1:l.coarse){
	DIC.coarse[i] <- probit.semipar.mcmc(y,X,Z,priors=priors,start=start,
		sigma.alpha=coarse.grid[i],n.mcmc=1000)$DIC
}

plot(coarse.grid,DIC.coarse,type="l")
idx <- which.min(DIC.coarse)

l.fine <- 10
fine.grid <- seq(coarse.grid[idx-1],coarse.grid[idx+1],length.out=10)
DIC.fine <- numeric(l.fine)
for(i in 1:l.fine){
	DIC.fine[i] <- probit.semipar.mcmc(y,X,Z,priors=priors,start=start,
		sigma.alpha=fine.grid[i],n.mcmc=1000)$DIC
}

plot(fine.grid,DIC.fine,type="l")
fine.grid[which.min(DIC.fine)]