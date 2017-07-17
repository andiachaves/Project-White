# Day 6 first model (Model 1). Model for variance term in normal likelihood 
# combines the observation error and process error variances. Thus, we do not
# specify an explict stochastic process model for the "latent" (unknown) 
# hydraulic conductivity values (this will be done in Model 2).

model{
	# Define the likelihood of the relative conductivity data (Ks)
	# Each relative conductivity measurement on segment i is distributed
	# as normal with mean mu and precision tau (std dev = sig = 1/sqrt(tau)),
	# which varies at the level of the observation (i, m):
	for(i in 1:N){ # N=total number of segments
		for(m in 1:M){ # M=number of repeated measurements per stem
			# Likelihood for observed data (precision [or variance] term varies 
			# by i and m to acount for stochastic embolism)
			Ks[i,m] ~ dnorm(mu[i,m],tau[i,m])

			# Define plc according to Weibull cdf, and compute mu given our model for plc:
			plc[i,m] <- 1 - exp(-p1[i,m])
			p1[i,m] <- pow(P[i,m]/alpha[i], beta[i])
			mu[i,m] <- ksat[i]*(1-plc[i,m])
		
			# Define variance associated with the likelihood for the Ks data; 
			# vO is the observation or measurement error variance, and vP 
			# describes the maximum process error variance, 
			# which occurs when plc = 0.5 (or 50% loss of conductivity):
			varK[i,m] <- vO + vP*(4*plc[i,m]*(1-plc[i,m]))
			# Compute overall precision for normal likelihood:
			tau[i,m] <- 1/varK[i,m]
		}
		
		# Assign semi-informative prior to the initial, unknown pressure
		P[i,1] ~ dunif(0,2)		

		# Hierarhical gamma priors for tree-level parameters, where there
		# gamma distributions depend on location-level parameters (e.g.,
		# a.beta, a.alpha,...,b.alpha, b.ksat):
		beta[i] ~ dgamma(a.beta[loc[i]], b.beta[loc[i]])
		alpha[i] ~ dgamma(a.alpha[loc[i]], b.alpha[loc[i]])
		ksat[i] ~ dgamma(a.ksat[loc[i]], b.ksat[loc[i]])
		
		# Compute tree-level quantities of interest (P50, S50):		
		P50[i] <- alpha[i]*pow(log(2),1/beta[i])
		S50[i] <- (beta[i]/P50[i])*exp(-pow(S1[i],beta[i]))*pow(S1[i],beta[i])
		S1[i] <- P50[i]/alpha[i]
	}


	# Relatively non-informative, independent priors for location-level parameters:
	for(i in 1:Nloc){	# Loop through Nloc locations
		# Hierarchical, relativley non-informative gamma priors for location-level 
		# parameters representing the first moment (i.e., means) of corresponding 
		# tree-level parameters:
		E.beta[i] ~ dgamma(0.01,0.01)
		E.alpha[i] ~ dgamma(0.01,0.01)
		E.ksat[i] ~ dgamma(0.01,0.01)

		# Via moment matching, compute gamma shape (a) and rate (b) parameters for
		# the location-level parameters:		
		a.beta[i] <- pow(E.beta[i],2)/V.beta
		b.beta[i] <- E.beta[i]/V.beta

		a.alpha[i] <- pow(E.alpha[i],2)/V.alpha
		b.alpha[i] <- E.alpha[i]/V.alpha

		a.ksat[i] <- pow(E.ksat[i],2)/V.ksat
		b.ksat[i] <- E.ksat[i]/V.ksat
		
		# Compute location-level quantities of interest (P50, S50):		
		E.P50[i] <- E.alpha[i]*pow(log(2),1/E.beta[i])
		E.S50[i] <- (E.beta[i]/E.P50[i])*exp(-pow(E.S1[i],E.beta[i]))*pow(E.S1[i],E.beta[i])
		E.S1[i] <- E.P50[i]/E.alpha[i]			
	}

	# Assign non-informative, independent priors to all standard deviation terms
	# and compute associated precision and variance:
	vO <- pow(sigO,2)
	vP <- pow(sigP,2)
	sigO ~ dunif(0,10)
	sigP ~ dunif(0,10)	
	# Uniform priors for standard deviations associated with hierarhical gamma priors
	sig.beta ~ dunif(0,100)
	sig.alpha ~ dunif(0,100)
	sig.ksat ~ dunif(0,100)
	# Compute associated variances for hierarchical gamma priors; these will be
	# used with the momemt matching technique to compute the gamma shape and scale
	# scale parameters of interest (i.e., the a's and b's).	
	V.beta <- pow(sig.beta,2)
	V.alpha <- pow(sig.alpha,2)
	V.ksat <- pow(sig.ksat,2)	
	# Store all standard deviations in one vector:	
	sigs[1] <- sigO		
	sigs[2] <- sigP	
	sigs[3] <- sig.beta
	sigs[4] <- sig.alpha
	sigs[5] <- sig.ksat					
}
