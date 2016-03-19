###
### Simulate binary data and fit model using probit.mixed.mcmc.R
###

rm(list=ls())
library(mvtnorm)
library(splines)

T <- 100  # number of observations per group
J <- 10  # number of groups
g <- rep(1:J,each=T)  # grouping variable

# Define covariates
time <- c(replicate(J,c(rgamma(1,shape=10.1,scale=4.5),
	cumsum(rgamma(T-1,shape=1.1,scale=4.5)))))  # time covariate
X <- cbind(1,time)
qX <- ncol(X)
X[,2] <- scale(X[,2])

mu.beta <- matrix(c(-1,1.25),,1)  # mean of betas
rho <- -0.25  # correlation between betas
Lambda <- diag(qX)*1  # variance-covariance of betas
Lambda[1,2] <- Lambda[2,1] <- Lambda[1,1]*Lambda[2,2]*rho

beta <- t(rmvnorm(J,mu.beta,Lambda))  # betas for each group
plot(t(beta))


# Basis expansion
int <- 100  # interval between knots
tapply(time,g,range)

# List of basis exansions, 1 per group
W <- sapply(1:J,function(x) bs(time[g==x],degree=3,intercept=FALSE,
	knots=seq(min(time[g==x]),max(time[g==x]),by=int)),simplify=FALSE)
qW <- unlist(lapply(W,ncol))

sigma.alpha <- 1
alpha <- lapply(qW,function(x) t(rmvnorm(1,rep(0,x),sigma.alpha^2*diag(x))))
trend <- sapply(1:J,function(x) W[[x]]%*%alpha[[x]])
matplot(trend,type="l",lty=1)
trend <- c(trend)

# Simulate data
beta.tmp <- t(beta[,g])
p <- pnorm(rowSums(X*beta.tmp)+trend)
hist(p);summary(p)

y <- rbinom(T*J,1,p)
table(y)
plot(p,y)

plot(time,y,ylim=range(c(trend,y,rowSums(X*beta.tmp))))
matplot(matrix(time,T,J),matrix(rowSums(X*beta.tmp),T,J),col=2,type="l",lty=2,add=TRUE)
matplot(matrix(time,T,J),matrix(trend,T,J),col=3,type="l",lty=2,add=TRUE)
matplot(matrix(time,T,J),matrix(rowSums(X*beta.tmp)+trend,T,J),
	col=4,type="l",lty=2,add=TRUE)

# Fit model
source('~/Documents/git/SemiReg/probit.semireg.mixed.mcmc.R')
start <- list(beta=beta,mu.beta=mu.beta,Lambda=Lambda,sigma.alpha=sigma.alpha,alpha=alpha)
priors <- list(sigma.beta=10,S0=diag(qX),nu=qX+1,r=2,q=1)
tune <- list(xi=0.01)
out1 <- probit.semireg.mixed.mcmc(y,X,g,W,priors,start,tune,adapt=TRUE,10000)

# Examine estimates for beta_j and alpha_j
g.idx <- 4  # group idx for plotting beta_j
matplot(out1$beta[,,g.idx],type="l",lty=1);abline(h=beta[,g.idx],col=1:qX,lty=2)
idx.tmp <- sample(1:qW[g.idx],3)
matplot(out1$alpha[[g.idx]][,idx.tmp],type="l",lty=1)
abline(h=alpha[[g.idx]][idx.tmp,],col=1:5,lty=2)

# Examine estimates for mu.beta
matplot(out1$mu.beta,type="l");abline(h=mu.beta,col=1:qX,lty=2)

# Examine estimates for Lambda
matplot(cbind(out1$Lambda[1,1,],out1$Lambda[1,2,]),type="l")
abline(h=c(Lambda[1,1],Lambda[1,2]),lty=2,col=1:qX)
mean(out1$Lambda[1,1,])
mean(out1$Lambda[2,2,])
mean(out1$Lambda[1,2,])

# Examine estimates for sigma.alpha
matplot(out1$sigma.alpha,type="l")
plot(apply(out1$sigma.alpha,1,mean),type="l")

# Examine estimates for v
boxplot(pnorm(out1$v),col=8,outline=FALSE)
points(y,col=3)

###
### Compare parameter estimates to those from other models
### 

# Compare to mixed effect model (without a nonparametric component)
source('~/Documents/git/GLMM/probit.glmm.mcmc.R', chdir = TRUE)
start <- list(beta=t(beta),mu.beta=mu.beta,Lambda=Lambda)
priors <- list(sigma.beta=5,S0=diag(qX),nu=qX+1)
out2 <- probit.glmm.mcmc(y,X,g,priors,start,10000)
matplot(out2$beta[,,g.idx],type="l",lty=1);abline(h=beta[,g.idx],col=1:qX,lty=2)
matplot(out2$mu.beta,type="l");abline(h=mu.beta,col=1:qX,lty=2)

# Compare to fixed effect model (individual-level model with nonparametric component)
source('~/Documents/git/SemiReg/probit.semireg.mcmc.R', chdir = TRUE)
start <- list(beta=c(beta[,g.idx]),alpha=c(alpha[[g.idx]]))
# hist(sqrt(1/rgamma(1000,1,,2)),breaks=100)
priors <- list(mu.beta=rep(0,qX),sigma.beta=10)
out3 <- probit.semireg.mcmc(y[g==g.idx],X[g==g.idx,],W[[g.idx]],
	priors=priors,start=start,sigma.alpha=sigma.alpha,n.mcmc=5000)
matplot(out3$beta,type="l");abline(h=beta[,g.idx],col=1:qX)

# Compare to GLM with probit link (individual-level model without nonparametric component)
source('~/Documents/git/GLM/probit.glm.mcmc.R', chdir = TRUE)
start <- list(beta=beta[,g.idx])
priors <- list(mu.beta=rep(0,qX),Sigma.beta=diag(qX)*100)
out4 <- probit.glm.mcmc(y[g==g.idx],X[g==g.idx,],priors,start,10000)
matplot(out4$beta,type="l");abline(h=beta[,g.idx],col=1:qX)

