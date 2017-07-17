#Function for deterministic model
g=function(alpha, gamma,kappa, x){
	#the observed independet variable is x
	mu = alpha*inv.logit(gamma+kappa*x)
	return(mu)
}




#Function to simulate data
get_data=function(alpha, gamma,kappa, sigma,n){
	par(mfrow=c(1,1))
	x=sort(runif(n,min=800, max = 3500))
	mu = g(alpha=alpha, gamma=gamma, kappa=kappa, x=x)
	plot(x,mu, typ="l")
	y=rgamma(length(mu), mu^2/sigma^2, mu/sigma^2)
	plot(x,y)
	lines(x,mu)

	model=nls(y ~ g(alpha=alpha,gamma=gamma,kappa=kappa, x=x), start=list(alpha=100,gamma=-10,kappa=.005))
	s=summary(model)
	p=coef(model)
	y.hat=g(alpha=p[1],gamma=p[2],kappa=p[3], x=x)

	lines(x,y.hat,col="red")
	legend(40,18, c("generating", "nls fit"), lty=c("solid", "solid"), col=c("black", "red"), bty="n")
	return(list(x=x, y=y, nls.alpha=p[1], nls.gamma=p[2], nls.kappa=p[3], nls.sigma = s$sigma, gen.alpha=alpha, gen.gamma=gamma, gen.kappa=kappa, gen.sigma=sigma))
}


#function to stet up storage for chains and name them--entries must be made in the function body
setup=function(n.iter,n.chain,parameter.names, dim.x){
	#set up storage for chains
	x=list()
	for(i in 1:length(parameter.names)) {
		x[[i]]=array(NA,dim=c(dim.x[i],n.iter,n.chain))
	}
	#assign parameter names to elements of the list
	names(x)=parameter.names
	#enter initial guesses at parameters here
	x$alpha[1,1,1]=100
	x$kappa[1,1,1] = .001
	x$gamma[1,1,1]=5
	x$sigma[1,1,1]=50
	#enter tuning parameters here
	x$tune=list(
		alpha=20,
		kappa = .0001,
		gamma=.75,
		sigma=5
	) #end of tune list
	#Enter support for parameters here. Choices are "non-negative", "real", "zero-one".
	x$support=list(
	alpha="non-negative",
	kappa = "non-negative",
	gamma="real",
	sigma="non-negative"
		) #end of support list
	return(x)
} #end of setup function


#Function for priors
prior=function(param,theta){  
	if(param == "alpha") return( dunif(theta, min=0, max=500, log=TRUE))
	if(param=="kappa") return(dunif(theta,min = 0, max = 1, log=TRUE))
	if(param=="gamma") return(dunif(theta,min=-100,max=100,log=TRUE))
	if(param=="sigma" ) return(dunif(theta,min=0,max=200, log=TRUE))
	}


#Likelihood function	

Like=function(y,x,alpha,gamma,kappa ,sigma){
	mu=g(alpha=alpha, gamma=gamma,kappa=kappa, x=x)
	#You may need to modify the likelihood function if negative values done't make sense for the data.g
	#LogL = dnorm(y, mean=mu, sd=sigma, log=TRUE)
	LogL = dgamma(y,mu^2/sigma^2, mu/sigma^2, log=TRUE)
	#browser()
	return(sum(LogL[is.finite(LogL)]))
	}

#Function to choose current or proposed value
choose=function(x, z, Like_z, Like_x, param, tune, support){
    numerator = Like_z + prior(param,theta=z)  # these are both logs so we add rather than multiply
    denominator = Like_x + prior(param,theta=x)
    q.ratio = q(theta=x,mu=z,tune=tune, type="density", support=support) / 		q(theta=z,mu=x,tune=tune, type="density", support=support)
    R =  exp(numerator - denominator) * q.ratio #because these are logs, we exponetiate the difference to ge the ratio.
    if(is.finite(R)){
    	if(R > runif(1,min=0,max=1)) new =z
   	    else new = x
    	return(new)
    	}
    else return(x)
}


#Proposal function
q=function(theta,mu,tune, type, support){
	sigma=tune
	if(support == "non-negative"){
		if (type == "density") return (dgamma(theta,mu^2/sigma^2, mu/sigma^2))
		if (type == "draw") return (rgamma(1,mu^2/sigma^2, mu/sigma^2))
	} #end of non-negative support block
	if(support == "real"){
		if (type == "density") return (dnorm(theta,mu,sigma))
		if (type == "draw") return (rnorm(1,mu,sigma))
	} #end of real support block
	if(support == "zero-one"){
		#do moment matching for beta distribution
		a <-(mu^2-mu^3-mu*sigma^2)/sigma^2
		 b <- (mu-2*mu^2+mu^3-sigma^2+mu*sigma^2)/sigma^2
		if (type == "density") return (dbeta(theta,a,b))
		if (type == "draw") return (rbeta(1,a,b))
	} #end of zero-one support block


}

Run_MCMC = function(x,tune,support,n.iter){
	
#set up list to count number of time proposal is accepted.  For most efficient results, should be about 50% for models with few parameters, about 25% for models with many.  These are just rules of thumb.
accept=list(
	alpha=0,
	gamma=0,
	kappa=0,
	sigma=0
)

#if you had more than one chain, you would set up an outer loop for each chain here.
for( i in 2:n.iter){
	
	 if(i%%1000==0) cat(i,"\n");flush.console() #iteration counter, prints every 1000th i
	#Update alpha
	z=q(mu=x$alpha[i-1],tune=tune$alpha, support=support$alpha, type="draw")
	Like_z = Like(y=data$y, x=data$x, alpha=z, kappa=x$kappa[i-1], gamma=x$gamma[i-1], sigma=x$sigma[i-1])
	Like_x = Like(y=data$y, x=data$x, alpha=x$alpha[i-1], kappa=x$kappa[i-1], gamma=x$gamma[i-1], sigma=x$sigma[i-1])
	x$alpha[i] =choose(x=x$alpha[i-1], z=z, Like_z=Like_z, Like_x=Like_x, param="alpha", tune=tune$alpha, support=support$alpha)
	#numerator of acceptance ratio
	if(x$alpha[i] != x$alpha[i-1]) accept$alpha=accept$alpha+1
	
	# # #Update gamma
	z = q(mu=x$gamma[i-1],tune=tune$gamma, type="draw", support=support$gamma)
	Like_z = Like(y=data$y, x=data$x, alpha=x$alpha[i], kappa=x$kappa[i-1], gamma=z, sigma=x$sigma[i-1])
	Like_x = Like(y=data$y, x=data$x, alpha=x$alpha[i], kappa=x$kappa[i-1], gamma=x$gamma[i-1], sigma=x$sigma[i-1])
	#browser()
	 x$gamma[i] =choose(x=x$gamma[i-1], z=z, Like_z=Like_z, Like_x=Like_x, param="gamma", tune=tune$gamma, support=support$gamma)
	if(x$gamma[i] != x$gamma[i-1]) accept$gamma=accept$gamma+1
	
	# # #Update kappa
	z = q(mu=x$kappa[i-1],tune=tune$kappa, type="draw", support=support$kappa)
	Like_z = Like(y=data$y, x=data$x, alpha=x$alpha[i], kappa=z, gamma=x$gamma[i], sigma=x$sigma[i-1])
	Like_x = Like(y=data$y, x=data$x, alpha=x$alpha[i], kappa=x$kappa[i-1], gamma=x$gamma[i], sigma=x$sigma[i-1])
	x$kappa[i] =choose(x=x$kappa[i-1], z=z, Like_z=Like_z, Like_x=Like_x, param="kappa", tune=tune$kappa, support=support$kappa)
	if(x$kappa[i] != x$kappa[i-1]) accept$kappa=accept$kappa+1
		
	# # #Update sigma
	#z = proposal(mu=x$sigma,tune=tune$sigma)
    z = q(mu=x$sigma[i-1],tune=tune$sigma, type="draw", support=support$sigma)
	Like_z = Like(y=data$y, x=data$x, alpha=x$alpha[i], kappa=x$kappa[i], gamma=x$gamma[i], sigma=z)
	Like_x = Like(y=data$y, x=data$x, alpha=x$alpha[i], kappa=x$kappa[i], gamma=x$gamma[i], sigma=x$sigma[i-1])    
	x$sigma[i] =choose(x=x$sigma[i-1], z=z, Like_z=Like_z, Like_x=Like_x, param="sigma", tune=tune$sigma, support=support$sigma)
	if(x$sigma[i] != x$sigma[i-1]) accept$sigma=accept$sigma + 1
	
	
	#estimate derived quantities
	x$y.hat[,i,1] = g(alpha=x$alpha[i],gamma=x$gamma[i],kappa=x$kappa[i],x=data$x)
	} #end of iteration loop
	
	
	a=accept$alpha/n.iter; 	print("alpha acceptance"); print(a) 
	a=accept$gamma/n.iter; 	print("gamma acceptance"); print(a) 
    a=accept$kappa/n.iter; 	print("kappa acceptance"); print(a) 
	a=accept$sigma/n.iter; 	print("sigma acceptance"); print(a) 

	
	
	return(x) #return the filled chain
	
}  #end of Run_MCMC function


do_plots=function(data, x, n.iter, burnin){
	
	#trace plots
	par(mfrow=c(2,2))
	plot(x$gamma,typ="l", xlab="Iteration"); abline(h=mean(x$gamma[burnin:n.iter]),col="red")
	plot(x$alpha,typ="l",, xlab="Iteration");abline(h = mean(x$alpha[burnin:n.iter]), col="red")
	plot(x$kappa,typ="l",, xlab="Iteration");abline(h = mean(x$kappa[burnin:n.iter]), col="red")
	plot(x$sigma,typ="l", xlab="Iteration"); abline(h=mean(x$sigma[burnin:n.iter]), col="red")

	par(mfrow=c(1,1))
	q.y.hat=apply(x$y.hat[,burnin:n.iter,1],1, function(x) quantile(x, c(.025,.5,.975)))
	plot(data$x,data$y,xlab="x", ylab="y", main="Inverse logit model")
	lines(data$x,q.y.hat[2,], col="orange", lwd="4")
	lines(data$x,q.y.hat[1,], lty = "dashed", col="orange" )
	lines(data$x,q.y.hat[3,], lty = "dashed", col="orange" )
	lines(data$x, g(alpha=data$nls.alpha, gamma=data$nls.gamma, kappa=data$nls.kappa, x=data$x), col="blue")
	legend(40,18, c("Median", "2.5% quantile", "97.5% quantile", "nlsfit"), lty=c("solid", "dashed", "dashed"),col=c("orange", "orange", "orange", "blue"),bty="n")


plot_density = function(p,v1,v2, param, burnin, n.iter){
	hist(p[burnin:n.iter],breaks=100, xlab="Value of parameter", freq=FALSE, main=param)
	abline(v=v1,col="red", lwd=2)
	abline(v=v2, col="blue", lwd=2)
	abline(v=median(p[burnin:n.iter]), col="orange", lwd=2)
	
}

par(mfrow=c(2,2))
plot_density(p=x$alpha,v1=data$gen.alpha, v2=data$nls.alpha, param=expression(alpha), burnin=burnin, n.iter=n.iter)
legend(40,.20,c("generating", "nls", "median"), lty=c("solid","solid", "solid"), col=c("red", "blue", "orange"), bty="n")
plot_density(p=x$gamma,v1=data$gen.gamma, v2=data$nls.gamma, param=expression(gamma), burnin=burnin, n.iter=n.iter)
plot_density(p=x$kappa,v1=data$gen.kappa, v2=data$nls.kappa, param=expression(kappa), burnin=burnin, n.iter=n.iter)
plot_density(p=x$sigma,v1=data$gen.sigma, v2=data$nls.sigma, param=expression(sigma), burnin=burnin, n.iter=n.iter)


} #end of plotting function



