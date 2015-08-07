probit.semipar.mcmc <- function(y,X,Z,priors,start,sigma.alpha=NULL,n.mcmc){
	
	###
	### Brian M. Brost (05 AUG 2015)
	### 'Mixed' effects model for normally distributed data
	###

	###
	### Model statement:
	### y~N(X%*%beta+Z%*%alpha,sigma^2*I)
	### beta~N(0,sigma.beta^2*I)
	### alpha~N(0,sigma.alpha^2*I)
	###	sigma^2~IG(r.sigma,q.sigma)
	###	
	
	###
	### Libraries and Subroutines
	###

	# library(pscl)  # rigamma for IG prior on sigma

	###
	###  Setup Variables 
	###
  
	n <- length(y)  # number of observations
	qX <- ncol(X)  # number of 'fixed' effects
	qZ <- ncol(Z)  # number of 'random' effects

	###
	### Starting values and priors
 	###

	# browser()
	beta <- matrix(start$beta,qX)
	# alpha <- matrix(start$alpha,qZ)
	sigma2 <- start$sigma^2
	sigma2.alpha <- ifelse(is.null(sigma.alpha),start$sigma.alpha^2,sigma.alpha^2)

	# Sigma.inv <- solve(sigma2*diag(n))

	mu.beta <- rep(0,qX)
	Sigma.beta <- diag(qX)*priors$sigma.beta^2
	Sigma.beta.inv <- solve(Sigma.beta)
	
	mu.alpha <- rep(0,qZ)
	Sigma.alpha <- diag(qZ)*sigma2.alpha
	Sigma.alpha.inv <- solve(Sigma.alpha)

	###
	### Create receptacles for output
	###
  
	beta.save <- matrix(0,n.mcmc,qX)  # coefficients for 'fixed' effects
	alpha.save <- matrix(0,n.mcmc,qZ)  # coefficients for 'random' effects
	sigma.save <- numeric(n.mcmc)  # standard deviation of data model
	sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of parameter model
	y.hat.save <- matrix(0,n,n.mcmc)
	D.bar.save <- numeric(n.mcmc)  # D.bar for DIC calculation
		
	###
	### Begin MCMC loop
	###
  
	for (k in 1:n.mcmc) {
    	if(k%%1000==0) cat(k,"");flush.console()

		###
		###  Sample alpha ('random' effects) 
		###

		A.inv <- solve(t(Z)%*%Z/(sigma2)+Sigma.alpha.inv)
		b <- t(t(y-X%*%beta)%*%Z/(sigma2))  # +mu.alpha%*%Sigma.alpha.inv
		alpha <- A.inv%*%b+chol(A.inv)%*%matrix(rnorm(qZ),qZ,1)

		###
		###  Sample beta ('fixed' effects) 
		###

		# browser()
		# A.inv <- solve(t(X)%*%Sigma.inv%*%X+Sigma.beta.inv)
		# b <- t(t(y-Z%*%alpha)%*%Sigma.inv%*%X+mu.beta%*%Sigma.beta.inv)
		A.inv <- solve(t(X)%*%X/(sigma2)+Sigma.beta.inv)
		b <- t(t(y-Z%*%alpha)%*%X/(sigma2))  # +mu.beta%*%Sigma.beta.inv
		beta <- A.inv%*%b+chol(A.inv)%*%matrix(rnorm(qX),qX,1)

		###
		###  Sample sigma (standard deviation of observation model) 
		###
		
		# Using rigamma in {pscl}
		# r.tmp <- sum((y-X%*%beta-Z%*%alpha)^2)/2+1/priors$r.sigma
		# q.tmp <- n/2+priors$q.sigma
		# sigma2 <- rigamma(1,q.tmp,r.tmp)
		
		# Using rgamm {base}
		r.tmp <- 1/(sum((y-X%*%beta-Z%*%alpha)^2)/2+1/priors$r.sigma)
		q.tmp <- n/2+priors$q.sigma
		sigma2 <- 1/rgamma(1,q.tmp,,r.tmp)

		###
		###  Save samples 
	    ###

		beta.save[k,] <- beta
		alpha.save[k,] <- alpha
    	sigma.save[k] <- sigma2    
    	sigma.alpha.save[k] <- sigma2.alpha
		y.hat.save[,k] <- X%*%beta+Z%*%alpha
		D.bar.save[k] <- -2*sum(dnorm(y,X%*%beta+Z%*%alpha,sqrt(sigma),log=TRUE))
	}

	###
	### Write output
	###

	# Convert variance to standard deviation for output
	sigma.save <- sqrt(sigma.save)    
   	sigma.alpha.save <- sqrt(sigma.alpha.save)

	#  Calculate DIC
	# browser()
	if(qX==1)  postbetamn <- mean(beta.save)
	if(qX>1)  postbetamn <- apply(beta.save,2,mean)
	postalphamn <- apply(alpha.save[,],2,mean)
	postsigmamn <- mean(sigma.save)
	D.hat <- -2*(sum(dnorm(y,X%*%postbetamn+Z%*%postalphamn,postsigmamn,log=TRUE)))
	D.bar <- mean(D.bar.save[])
	pD <- D.bar-D.hat
	DIC <- D.hat+2*pD
  
	# keep$sigma <- keep$sigma/n.mcmc
	# cat(paste("\nsigma acceptance rate:",round(keep$sigma,2))) 
	list(beta=beta.save,alpha=alpha.save,sigma=sigma.save,sigma.alpha=sigma.alpha.save,
  		y.hat=y.hat.save,DIC=DIC,n.mcmc=n.mcmc)
}