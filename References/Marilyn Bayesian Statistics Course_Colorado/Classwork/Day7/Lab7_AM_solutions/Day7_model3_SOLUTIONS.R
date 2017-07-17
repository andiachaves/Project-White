#### SOLUTIONS:
#### Model 3: Population-level parameters are treated as root nodes.
#### No species-level parameters.


model{
	# Loop through each observation i in the dataset:
	for(i in 1:N){
		# The likelihood (sampling distribution for loglength and logmass) is defined by a
		# multivariate normal distribution with mean mu and precision matrix Omega. 
		# Note, JAGS parameterizes the normal distribution in terms of a precision matrix.
		# The data are Y[,1] = loglength; Y[,2] = logmass.
		Y[i,1:2] ~ dmnorm(mu[i,1:2], Omega[1:2,1:2])
		# Define replicated data, Yrep:
		Yrep[i,1:2] ~ dmnorm(mu[i,1:2], Omega[1:2,1:2])
		for(k in 1:Nvars){
			# Compute squared difference (or squared error) for posterior
			# predictive loss calculation:
			sqdiff[i,k] <- pow(Yrep[i,k] - Y[i,k],2)
			}
						
		# Define the mean vector (i.e., scaling model that relates the true or latent variables).
		# alpha is the species-specific normalizing constant, and beta is the species-specific
		# scaling exponent.
		for(k in 1:Nvars){
			mu[i,k] <- alpha[k] + beta[k]*(LogD[i]-mean(LogD[]))
			}
		} # close observation (i) loop

		# Compute posterior predictive loss for each trait:
		for(k in 1:Nvars){
			# Sum of squared diff for each trait variable (sum across observations):
			Dsum[k] <- sum(sqdiff[,k])
				}
		
		# Model specification for population-level scaling parameters, based on the scaling model.
		# Scaling exponents are allowed to vary (independently) by model
			# Hierarchical priors for species-specific scaling exponents:			
		for(k in 1:Nvars){
			beta[k] ~ dnorm(0,0.00001)
			alpha[k] ~ dnorm(0,0.00001)
			}	

		# Fairly non-informative, conjugate Wishart priors for precision matrix,
    # and compute the covariance matrix (Sigma):
		Omega[1:2,1:2] ~ dwish(R[1:2,1:2], 2)
		Sigma[1:2,1:2] <- inverse(Omega[1:2,1:2])
		# Compute correlations between loglength and logmass residual errors:
		rho[1] <- Sigma[1,2]/sqrt(Sigma[1,1]*Sigma[2,2]) 
		# Compute standard deviations of residuals (errors) for each variable
		sig[1] <- sqrt(Sigma[1,1])  # for loglength
		sig[2] <- sqrt(Sigma[2,2])  # for logmass

}