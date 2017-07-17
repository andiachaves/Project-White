#Function for deterministic model
g=function(alpha, gamma,c, L){
	#the observed independet variable is w
	mu =alpha*(L - c) / (alpha/gamma + (L-c))
	return(mu)
}

#Function to simulate data
get_data=function(alpha, gamma,c, sigma){
	set.seed(4)
	par(mfrow=c(1,1))
	L=sort(runif(50,min=12, max = 100))
	mu = g(alpha=alpha, gamma=gamma, c=c, L=L)
	plot(L,mu, typ="l")
	y=rgamma(length(mu), mu^2/sigma^2, mu/sigma^2)
	plot(L,y)
	lines(L,mu)

	model=nls(y ~ g(alpha=alpha,gamma=gamma,c=c, L=L), start=list(alpha=50,gamma=4,c=2))
	s=summary(model)
	p=coef(model)
	y.hat=g(alpha=p[1],gamma=p[2],c=p[3], L=L)

	lines(L,y.hat,col="red")
	legend(40,18, c("generating", "nls fit"), lty=c("solid", "solid"), col=c("black", "red"), bty="n")
	return(list(x=L, y=y, nls.alpha=p[1], nls.gamma=p[2], nls.c=p[3], nls.sigma = s$sigma, gen.alpha=alpha, gen.gamma=gamma, gen.c=c, gen.sigma=sigma))
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
	x$alpha[1,1,1]=60 
	x$c[1,1,1] = 10
	x$gamma[1,1,1]=3
	x$sigma[1,1,1]=5
	#enter tuning parameters here
	tune=list(
		alpha=10,
		c = 1,
		gamma=.3,
		sigma=2
	) #end of tune list
	x$tune = tune
	return(x)
} #end of setup function


#Function for priors
prior=function(param,theta){  
	if(param == "alpha") return( dunif(theta, min=0, max=500, log=TRUE))
	if(param=="c") return(dunif(theta,min = 0, max = 200, log=TRUE))
	if(param=="gamma") return(dunif(theta,min=0,max=200,log=TRUE))
	if(param=="sigma" ) return(dunif(theta,min=0,max=200, log=TRUE))
	}


#Likelihood function	
Like=function(y,L,alpha,gamma,c ,sigma){
	mu=g(alpha=alpha, gamma=gamma,c=c, L=L)
	LogL = dnorm(y, mean=mu, sd=sigma, log=TRUE)
	return(sum(LogL[is.finite(LogL)]))
	}

#Function to choose current or proposed value
choose=function(x, z, Like_z, Like_x, param, tune){
    numerator = Like_z + prior(param,theta=z)  # these are both logs so we add rather than multiply
    denominator = Like_x + prior(param,theta=x)
    q.ratio = q(theta=x,mu=z,tune=tune, type="density") / q(theta=z,mu=x,tune=tune, type="density")
    R =  exp(numerator - denominator) * q.ratio #because these are logs, we exponetiate the difference to ge the ratio.
    if(R > runif(1,min=0,max=1)) new =z
    else new = x
    	return(new)
}


#Proposal function
q=function(theta,mu,tune, type){
	sigma=tune
	if (type == "density") return (dgamma(theta,mu^2/sigma^2, mu/sigma^2))
	if (type == "draw") return (rgamma(1,mu^2/sigma^2, mu/sigma^2))
}

Run_MCMC = function(x,tune,n.iter){
for( i in 2:n.iter){
	 if(i%%1000==0) cat(i,"\n");flush.console() #iteration counter, prints every 1000th i
	#Update alpha
	#z = proposal(mu=x$alpha[i-1],tune=tune$alpha)
	z=q(mu=x$alpha[i-1],tune=tune$alpha, type="draw")
	Like_z = Like(y=data$y, L=data$x, alpha=z, c=x$c[i-1], gamma=x$gamma[i-1], sigma=x$sigma[i-1])
	Like_x = Like(y=data$y, L=data$x, alpha=x$alpha[i-1], c=x$c[i-1], gamma=x$gamma[i-1], sigma=x$sigma[i-1])
	x$alpha[i] =choose(x=x$alpha[i-1], z=z, Like_z=Like_z, Like_x=Like_x, param="alpha", tune=tune$alpha)
	
	# # #Update gamma
	 #z = proposal(mu=x$gamma[i-1],tune=tune$gamma)
	  z = q(mu=x$gamma[i-1],tune=tune$gamma, type="draw")
	Like_z = Like(y=data$y, L=data$x, alpha=x$alpha[i], c=x$c[i-1], gamma=z, sigma=x$sigma[i-1])
	Like_x = Like(y=data$y, L=data$x, alpha=x$alpha[i], c=x$c[i-1], gamma=x$gamma[i-1], sigma=x$sigma[i-1])
	 x$gamma[i] =choose(x=x$gamma[i-1], z=z, Like_z=Like_z, Like_x=Like_x, param="gamma", tune=tune$gamma)
	
	# # #Update c
	#z = proposal(mu=x$c[i-1],tune=tune$c)
	z = q(mu=x$c[i-1],tune=tune$c, type="draw")
	Like_z = Like(y=data$y, L=data$x, alpha=x$alpha[i], c=z, gamma=x$gamma[i], sigma=x$sigma[i-1])
	Like_x = Like(y=data$y, L=data$x, alpha=x$alpha[i], c=x$c[i-1], gamma=x$gamma[i], sigma=x$sigma[i-1])
	x$c[i] =choose(x=x$c[i-1], z=z, Like_z=Like_z, Like_x=Like_x, param="c", tune=tune$c)
	
	# # #Update sigma
	#z = proposal(mu=x$sigma,tune=tune$sigma)
    z = q(mu=x$sigma[i-1],tune=tune$sigma, type="draw")
	Like_z = Like(y=data$y, L=data$x, alpha=x$alpha[i], c=x$c[i], gamma=x$gamma[i], sigma=z)
	Like_x = Like(y=data$y, L=data$x, alpha=x$alpha[i], c=x$c[i], gamma=x$gamma[i], sigma=x$sigma[i-1])    
	x$sigma[i] =choose(x=x$sigma[i-1], z=z, Like_z=Like_z, Like_x=Like_x, param="sigma", tune=tune$sigma)
	
	#estimate derived quantities
	x$y.hat[,i,1] = g(alpha=x$alpha[i],gamma=x$gamma[i],c=x$c[i],L=data$x)
	x$growth.ratio[i]= x$alpha[i] / x$gamma[i]
	} #end of iteration loop
	return(x) #return the filled chain
	
}  #end of Run_MCMC function


do_plots=function(data, x, n.iter, burnin){
	
	#trace plots
	par(mfrow=c(2,2))
	plot(x$gamma,typ="l", xlab="Iteration"); abline(h=mean(x$gamma[burnin:n.iter]),col="red")
	plot(x$alpha,typ="l",, xlab="Iteration");abline(h = mean(x$alpha[burnin:n.iter]), col="red")
	plot(x$c,typ="l",, xlab="Iteration");abline(h = mean(x$c[burnin:n.iter]), col="red")
	plot(x$sigma,typ="l", xlab="Iteration"); abline(h=mean(x$sigma[burnin:n.iter]), col="red")

	par(mfrow=c(1,1))
	q.y.hat=apply(x$y.hat[,burnin:n.iter,1],1, function(x) quantile(x, c(.025,.5,.975)))
	plot(data$x,data$y,xlab="Light level", ylab="Growth rate", main="Prediction of growth rate")
	lines(data$x,q.y.hat[2,], col="orange", lwd="4")
	lines(data$x,q.y.hat[1,], lty = "dashed", col="orange" )
	lines(data$x,q.y.hat[3,], lty = "dashed", col="orange" )
	lines(data$x, g(alpha=data$nls.alpha, gamma=data$nls.gamma, c=data$nls.c, L=data$x), col="blue")
	legend(40,18, c("Median", "2.5% quantile", "97.5% quantile", "nlsfit"), lty=c("solid", "dashed", "dashed"),col=c("orange", "orange", "orange", "blue"),bty="n")


plot_density = function(p,v1,v2, param, burnin, n.iter){
	hist(p[burnin:n.iter],breaks=100, xlab="Value of parameter", freq=FALSE, main=param)
	abline(v=v1,col="red", lwd=2)
	abline(v=v2, col="blue", lwd=2)
	abline(v=median(p[burnin:n.iter]), col="orange", lwd=2)
}

par(mfrow=c(2,2))
plot_density(p=x$alpha,v1=data$gen.alpha, v2=data$nls.alpha, param=expression(alpha), burnin=burnin, n.iter=n.iter)
plot_density(p=x$gamma,v1=data$gen.gamma, v2=data$nls.gamma, param=expression(gamma), burnin=burnin, n.iter=n.iter)
plot_density(p=x$c,v1=data$gen.c, v2=data$nls.c, param=expression(c), burnin=burnin, n.iter=n.iter)
plot_density(p=x$sigma,v1=data$gen.sigma, v2=data$nls.sigma, param=expression(sigma), burnin=burnin, n.iter=n.iter)

} #end of plotting function



