probit.semipar.mcmc <- function(y,X,Z,priors,start,sigma.alpha,n.mcmc){
	
	###
	### Brian M. Brost (10 AUG 2015)
	### Semiparametric regression for binary data using probit link
	###

	###
	### Model statement:
	### y=0,u_t<=0
	### y=1,u_t>0
	### u_t~N(x*beta+z*alpha,1)
	### beta~N(0,sigma.beta^2*I)
	### alpha~N(0,sigma.alpha^2*I)
	###	
	
	###
	### Libraries and Subroutines
	###

	truncnormsamp <- function(mu,sig2,low,high,nsamp){
	  flow=pnorm(low,mu,sqrt(sig2)) 
	  fhigh=pnorm(high,mu,sqrt(sig2)) 
	  u=runif(nsamp) 
	  tmp=flow+u*(fhigh-flow)
	  x=qnorm(tmp,mu,sqrt(sig2))
	  x
	}

	###
	###  Setup Variables 
	###
  
	n <- length(y)  # number of observations
	qX <- ncol(X)  # number of 'fixed' effects
	qZ <- ncol(Z)  # number of 'random' effects
	y1 <- (y==1)
	y0 <- (y==0)
	y1.sum <- sum(y1)
	y0.sum <- sum(y0)
	u <- numeric(n)
	
	###
	### Starting values and priors
 	###

	# browser()
	beta <- matrix(start$beta,qX)
	alpha <- matrix(start$alpha,qZ)

	mu.beta <- priors$mu.beta
	Sigma.beta <- diag(qX)*priors$sigma.beta^2
	Sigma.beta.inv <- solve(Sigma.beta)
	
	mu.alpha <- rep(0,qZ)
	Sigma.alpha <- diag(qZ)*sigma.alpha^2
	Sigma.alpha.inv <- solve(Sigma.alpha)

	###
	### Create receptacles for output
	###
  
	beta.save <- matrix(0,n.mcmc,qX)  # coefficients for 'fixed' effects
	alpha.save <- matrix(0,n.mcmc,qZ)  # coefficients for 'random' effects
	sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of parameter model
	u.save <- matrix(0,n.mcmc,n)
	D.bar.save <- numeric(n.mcmc)  # D.bar for DIC calculation
		
	###
	### Begin MCMC loop
	###
  
	for (k in 1:n.mcmc) {
    	if(k%%1000==0) cat(k,"");flush.console()

# browser()
###
		### Sample u (auxilliary variable for probit regression)
	  	###
	 
		linpred <- X%*%beta+Z%*%alpha
	  	u[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	u[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)

		###
		###  Sample alpha ('random' effects) 
		###
# browser()
		A.inv <- solve(t(Z)%*%Z+Sigma.alpha.inv)
		b <- t(Z)%*%(u-X%*%beta)  # +mu.alpha%*%Sigma.alpha.inv
		alpha <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qZ),qZ,1)

		###
		###  Sample beta ('fixed' effects) 
		###

	 	A.inv <- solve(t(X)%*%X+Sigma.beta.inv)
	  	b <- t(X)%*%(u-Z%*%alpha)  # +mu.beta%*%Sigma.beta.inv
	  	beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)
		
		###
		###  Save samples 
	    ###

		beta.save[k,] <- beta
		alpha.save[k,] <- alpha
		u.save[k,] <- u
	  	D.bar.save[k] <- -2*(sum(dbinom(y,1,pnorm(u),log=TRUE)))
	}

	#  Calculate DIC
	if(qX==1)  postbetamn <- mean(beta.save)
	if(qX>1)  postbetamn <- apply(beta.save,2,mean)
	postumn <- apply(u.save,2,mean)
	D.hat=-2*(sum(dbinom(y,1,pnorm(postumn),log=TRUE)))
	D.bar <- mean(D.bar.save)
	pD <- D.bar-D.hat
	DIC <- D.hat+2*pD

	###
	### Write output
	###
  
	list(beta=beta.save,alpha=alpha.save,u=u.save,DIC=DIC,n.mcmc=n.mcmc)
}