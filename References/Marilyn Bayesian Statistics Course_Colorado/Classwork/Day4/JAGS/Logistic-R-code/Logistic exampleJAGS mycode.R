##JAGS script
model { 
##Specify priors: 
K ~ dgamma(.001,.001)
r ~ dgamma(.001,.001)
tau ~ dgamma(.001,.001) # precision
sigma<-1/sqrt(tau) #calculate sd from precision

#likelihood

for(i in 1:n){
	mu[i] <- r - r/K * x[i]
	y[i] ~ dnorm(mu[i],tau)
	}

} #end of model

