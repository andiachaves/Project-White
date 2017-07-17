model
{
##priors
#Need priors for N1, phi, lambda, H, sigma.proc
sigma.proc ~ dunif(0, 5)
tau.proc <- 1/sigma.proc^2

lambda ~ dunif(0.5, 10)
phi ~ dbeta(y.a, y.b) # phi needs to be between 0 and 1 because it is #familygroups/N
					# phi a and b are calculated from moment matching because we know the 
					# mean and SD of phi 

#Initial conditions on first observation of family groups
fg[1] ~ dpois(y1)  # y1 is the first observation in the census data.
N[1] ~ dpois(fg[1]/phi) #The first N is the fg divided by phi

##Process model
for (t in 2:(y.endyr + 1)) {
	mu[t] <- log(max(0.001, lambda*(N[t-1]-y.H[t-1])))  ## This has a max function on it bc this equation could go
										#negative and you can't take the log of a negative. So it sets the
										# minimum at 0.001
	N[t] ~ dlnorm(mu[t], tau.proc)
	fg[t] <- N[t] * phi
} #end process model

##Data model
for(t in 2:y.endyr) {
	y[t] ~ dpois(phi*N[t])
}

#simulate new data for posterior predicitve check
for(t in 1:y.endyr){
	    y.rep[t] ~ dpois(phi*N[t])
	    #accumlate test statistics for posterior predictive check
	    sq[t] <- (y[t]-phi*N[t])^2
	    sq.rep[t] <-(y.rep[t] - phi*N[t])^2
}

#calculate Bayesian P value
fit <- sum(sq[])
fit.new <- sum(sq.rep[])
pvalue <- step(fit-fit.new)

	
##forecast effects of different harvest regeimes on next years
#number of family grops
	for(i in 1:length(harv)){
		#mu.hat is the forecast 1 year beyond y.endyr +1, i.e., 2011
		mu.hat[i] <- log(max(.001,lambda*(N[y.endyr+1]-harv[i]))) 
		N.hat[i] ~ dlnorm(mu.hat[i], tau.proc)	#Nhat forecasts 2 years out
		fg.hat[i] <- N.hat[i] * phi
    }

} #end of model

