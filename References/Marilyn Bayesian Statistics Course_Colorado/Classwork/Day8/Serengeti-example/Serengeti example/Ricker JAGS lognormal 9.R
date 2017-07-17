

model{
#Priors
for(j in 1:n.coef){
	b[j] ~ dnorm(0,.0001)
}
tau.p ~ dgamma(.01,.01)
sigma.p <- 1/sqrt(tau.p)

# N0 ~ dunif(100,1000)
# N[1] <- N0
# logN[1] ~ dunif(0,100)

z[1] ~ dunif(0,500)

#give starting value to epsilon.proc to allow monitoring
epsilon.proc[1] <- 0

##Process model
# for(t in 2:T){
	# mu[t] <- log(N[t-1]*exp(b[1] + b[2]*N[t-1] + b[3]*Rain[t] +b[4]*Rain[t]*N[t-1]))
	# logN[t] ~ dnorm(mu[t], tau.p)
	# N[t] <- min(exp(logN[t]),50000)
	# }
	
	
	for(t in 2:T){
	mu[t] <- log(z[t-1]*exp(b[1] + b[2]*z[t-1] + b[3]*Rain[t] +b[4]*Rain[t]*z[t-1]))
	z[t] ~ dlnorm(mu[t], tau.p)
	}
#Data model
for(j in 1:n.obs){
		N.obs[j] ~ dnorm(z[index[j]],tau.obs[j])
}

#Derived quantities for model evaluation

for(t in 2:(T-4)){  #T-4 to remove effect of initial conditions and forecasting
	epsilon.proc[t] <- log(z[t]/z[t-1])  #to test for auto correlation in process
}
for(i in 1:n.obs){
		epsilon.obs[i] <- N.obs[i] -z[index[i]]  #to test for auto correlation in data
		N.new[i] ~ dnorm(z[index[i]],tau.obs[i])  # simulate new data for N
		sq[i] <- (N.obs[i] - z[index[i]] )^2
		sq.new[i] <-(N.new[i] - z[index[i]]) ^2 
}
fit  <- sum(sq[])
fit.new <- sum(sq.new[])
pvalue <-step(fit.new-fit)



}#end of model
