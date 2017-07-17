#############################
#Jamaica bird occupancy model 
#"First Attempt"
#############################

model
{
	b0 ~ dnorm(0, .0001)
	b1 ~ dnorm(0, .0001)
	b2 ~ dnorm(0, .0001)
	b3 ~ dnorm(0, .0001)
	b4 ~ dnorm(0, .0001)
	b5 ~ dnorm(0, .0001)
	b6 ~ dnorm(0, .0001)

	a0 ~ dnorm(0, .0001)
	a1 ~ dnorm(0, .0001)
	a2 ~ dnorm(0, .0001)
	a3 ~ dnorm(0, .0001)

	for (i in 1:99){  #For each patch
		for (j in 1:3) { #For each survey
			p[i,j] <- ilogit(a0 + a1(effort[i]) + a2(srvy[2,j]) + a3 (srvy[3,j])) #Calculate detection probability
		}
	}
	
	for (i in 1:99){ #For each patch
		psi[i] <- ilogit(b0 + b1*ag[i] + b2*urban[i] + b3*mining[i] + b4*veb[i] + b5*size[i] + b6*isol[i])) #calculate occupancy probability
		z[i] ~ dbern(psi[i]) #use occupancy probability to calculate z
	}


	for (i in 1:99) {
		for (j in 1:3){
			y[i,j] ~ dbern(p[i,j]*z[i]) # Calculate y for each i and j.
		}
	}

}	# end of model


