probit.semireg.mixed.mcmc <- function(y,X,g,W,priors,start,n.mcmc){

	###
	### Brian M. Brost (09 MAR 2016)
	### Mixed effects regression for binary data using probit link
	###

	###
	### Model statement:
	### y_j(t)=0,v_j(t)<=0
	### y_j(t)=1,v_j(t)>0
	### v_ij ~ N(x_j(t)*beta_j+w_j(t)*alpha_j,1)
	### beta_j ~ N(mu.beta,Lambda)
	###	mu.beta ~ N(0,sigma.beta^2*I)
	### alpha ~ N(0,sigma.alpha.j^2*I)
	###	Lambda ~ Wish(S0,nu)
	### sigma.alpha.j ~ IG(r,q)
	###	
	
	###
	###  Libraries and Subroutines
	###

	# library(MCMCpack)  # for function rwish 
	
	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# a <- min(0.01,1/sqrt(k))
		a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}
	
	truncnormsamp <- function(mu,sig2,low,high,nsamp){
	  flow=pnorm(low,mu,sqrt(sig2)) 
	  fhigh=pnorm(high,mu,sqrt(sig2)) 
	  u=runif(nsamp) 
	  tmp=flow+u*(fhigh-flow)
	  x=qnorm(tmp,mu,sqrt(sig2))
	  x
	}
	
	###
	###  Preliminary Variables
	###

# browser()
	J <- length(unique(g))  # number of groups
	g <- as.numeric(g)
	g.idx <- sapply(sort(unique(g)),function(x) which(g==x),simplify=FALSE)
		# indexes of observations in y by group
	n.j <- unlist(lapply(g.idx,length))  # number of observations per group
	n <- length(y)  # total number of observations
	qX <- ncol(X)
	qW <- unlist(lapply(W,ncol))
	y1 <- (y==1)
	y0 <- (y==0)
	y1.sum <- sum(y1)
	y0.sum <- sum(y0)
	v <- numeric(n)

	###
	### Starting values and priors
 	###

# browser()
	beta <- matrix(start$beta,qX,J)
	mu.beta <- matrix(start$mu.beta,qX,1)
	Lambda <- start$Lambda
	Lambda.inv <- solve(Lambda)
	alpha <- lapply(qW,function(x) rep(0,x))
# xi <- rep(1.0,qX)
	Sigma.alpha <- lapply(qW,function(x) diag(x)*start$sigma.alpha^2)
	Sigma.alpha.inv <- lapply(Sigma.alpha,solve)
	W.cross <- lapply(W,crossprod)  # t(W)%*%Wcross product of W

	Sigma.beta <- priors$sigma.beta^2*diag(qX)
	Sigma.beta.inv <- solve(Sigma.beta)
	S0 <- priors$S0
	nu <- priors$nu
	r <- priors$r
	q <- priors$q

	###
	### Create receptacles for output
	###
  
	beta.save <- array(0,dim=c(n.mcmc,qX,J))
	mu.beta.save <- matrix(0,n.mcmc,qX)
# xi.save <- matrix(0,n.mcmc,qX)
	alpha.save <- lapply(qW,function(x) matrix(0,n.mcmc,x))
	Lambda.save <- array(0,dim=c(qX,qX,n.mcmc))
	sigma.alpha.save <- matrix(0,n.mcmc,J)
	v.save <- matrix(0,n.mcmc,n)

	keep <- list(xi=0)
# keep.tmp <- keep
# T.b <- 50  # frequency of adaptive tuning

	###
	###  Gibbs Loop
	###
	
	for(k in 1:n.mcmc){
		if(k%%1000==0) cat(k," ");flush.console()

		# if(adapt==TRUE & k%%T.b==0) {  # Adaptive tuning
			# keep.tmp$xi <- keep.tmp$xi/(T.b*qX)
			# tune$xi <- get.tune(tune$xi,keep.tmp$xi,k)
			# keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	# } 	

		for(i in 1:J){

			###
			### Update alpha_j
			###
# browser()		
			idx <- g.idx[[i]]
			b <- crossprod(W[[i]],(v[idx]-X[idx,]%*%beta[,i]))
			A.inv <- solve(W.cross[[i]]+Sigma.alpha.inv[[i]])
			alpha[[i]] <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qW[i]),qW[i],1)
			# alpha[[i]] <- t(rmvnorm(1,A.inv%*%b,A.inv))
			alpha.save[[i]][k,] <- alpha[[i]]
		  	
		  	###
		  	### Sample beta_j
	  		###
	  		
			A.inv <- solve(t(X[idx,])%*%X[idx,]+Lambda.inv) 
			b <- crossprod(X[idx,],(v[idx]-W[[i]]%*%alpha[[i]]))+Lambda.inv%*%mu.beta
			beta[,i] <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)
			# beta[,i] <- t(rmvnorm(1,A.inv%*%b,A.inv))
		}

		###
		### Sample v (auxilliary variable for probit regression)
	  	###

		linpred <- rowSums(X*t(beta[,g]))+
			unlist(sapply(1:J,function(x) W[[x]]%*%alpha[[x]]))
		v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)

		###
		### Update sigma.alpha_j
		###

		r.tmp <- 1/(unlist(lapply(alpha,function(x) sum(x^2)))/2+1/r)
		q.tmp <- qW/2+q
		sigma2.alpha <- 1/sapply(r.tmp,function(x) rgamma(1,q.tmp,,x))
		for(i in 1:J) diag(Sigma.alpha[[i]]) <- sigma2.alpha[i]
		Sigma.alpha.inv <- lapply(Sigma.alpha,solve)
			
	  	###
	  	### Sample mu_beta
	  	###

		beta.mean <- apply(beta,1,sum)
		A.inv <- solve(J*Lambda.inv+Sigma.beta.inv)
		b <- Lambda.inv%*%beta.mean
	    mu.beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)
		# mu.beta <- t(rmvnorm(1,A.inv%*%b,A.inv))
	  	
	  	###
	  	### Sample Lambda
	  	###
# browser()		
	  	Sn <- S0+crossprod(t(beta)-matrix(mu.beta,J,qX,byrow=TRUE))
		Lambda <- solve(rWishart(1,nu+J,solve(Sn))[,,1])
		Lambda.inv <- solve(Lambda)

		###
		### Sample Lambda from scaled-inverse Wishart 
		### 

		# Sample unscaled covariance matrix
		# Q <- solve(rWishart(1,nu+J,solve(Sn))[,,1])

		# Sample xi (scaling factors)
		# xi.star <- rnorm(qX,xi,tune$xi)
		# idx <- which(xi>0 & xi<100)
		# for(i in idx){
			# xi.tmp <- xi
			# xi.tmp[i] <- xi.star[i]
			# Lambda.star <- (diag(qX)*xi.tmp)%*%Q%*%(diag(qX)*xi.tmp)

			# # Q[1,1]*xi.star[1]^2
			# # Q[1,2]*xi.star[1]*xi[2]
			# # Q[2,2]*xi[2]^2
			# # Q[2,2]*xi.star[2]^2
			# # Q[1,2]*xi[1]*xi.star[2]

			# mh.star.xi <- sum(dmvnorm(t(beta),mu.beta,Lambda.star,log=TRUE))
			# mh.0.xi <- sum(dmvnorm(t(beta),mu.beta,Lambda,log=TRUE))

			# if(exp(mh.star.xi-mh.0.xi)>runif(1)){
	        	# xi[i] <- xi.star[i]
				# Lambda <- Lambda.star
				# keep$xi <- keep$xi+1
				# keep.tmp$xi <- keep.tmp$xi+1
	    	# } 
		# }		

		# Lambda.inv <- solve(Lambda)
				

	  	###
	  	### Save Samples 
	  	###

	  	beta.save[k,,] <- beta
	  	v.save[k,] <- v
		mu.beta.save[k,] <- mu.beta
		sigma.alpha.save[k,] <- sqrt(sigma2.alpha)
		Lambda.save[,,k] <- Lambda
		# xi.save[k,] <- xi
	}
	
	# keep$xi <- keep$xi/(qX*n.mcmc)
	# cat(paste("\nxi acceptance rate:",round(keep$xi,2))) 
	# cat(paste("\nxi tuning parameter:",round(tune$xi,2))) 
	
	###
	###  Write output 
	###
	
	list(beta=beta.save,v=v.save,mu.beta=mu.beta.save,Lambda=Lambda.save,
		alpha=alpha.save,sigma.alpha=sigma.alpha.save,#xi=xi.save,keep=keep,
		y=y,X=X,g=g,W=W,start=start,priors=priors,n.mcmc=n.mcmc)
}
