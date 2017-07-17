model
{
	b0 ~ dnorm(0, .0001)
	b1 ~ dnorm(0, .0001)
	b2 ~ dnorm(0, .0001)
	b3 ~ dnorm(0, .0001)
	p.detect ~ dunif(0,1)

	for (i in 1:237)
	{
		phi[i] <- ilogit(b0 + b1*elev[i] + b2*forest[i] + b3*elev[i]^2)
		z[i] ~ dbern(phi[i])
		y[i] ~ dbin(p.detect*z[i], n[i])
	}

	
	#elev.scaled.opt <- -(1/2)*b1/b2
	#elev.opt <- elev.scaled.opt * sd.elev + mu.elev
	
	#for (i in 1:237) 
	#{
	#logit(phi.elev[j]) <- b0 + b1*(elev.x[j]-mu.elev)/sd.elev + b2*((elev.x[j]-mu.elev)/sd.elev)^2
	#}


}	# end of model


