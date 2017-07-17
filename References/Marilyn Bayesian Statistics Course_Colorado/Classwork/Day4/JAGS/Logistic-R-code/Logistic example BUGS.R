##Logistic example for Primer
model{
#priors
#K~dgamma(4,.004)
K~dgamma(.001,.001)
#K~dunif(200,4000)  #produced ca identical results

r~dgamma(.001,.001)
tau~dgamma(.001,.001)
sigma<-1/sqrt(tau)

#likelihood
for(i in 1:n){
	mu[i] <- r - r/K * x[i]
	y[i] ~ dnorm(mu[i],tau)
	}
	
} #end of model



