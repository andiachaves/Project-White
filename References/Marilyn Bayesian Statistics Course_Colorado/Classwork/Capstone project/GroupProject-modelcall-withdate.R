setwd("C:\\Users\\MBrandt\\Dropbox\\Transfer documents\\conferences\\Bayesian Statistics\\Classwork\\Capstone project\\")

#Load jags
library(rjags)

#Read in data
y = read.csv("Jamaicabirddata.csv")

#Create data list
bird_data = list(N = 99, effort = y[,X], date = y[,X], ag = y[,X], urban  = y[,X], 
mining = y[,X], veg  = y[,X], size = y[,X], isol = y[,X])

#Code the ag, urban, mining as zeros and ones. 

#Set initial conditionals to initialize MCMC chains. 

inits = list(z=zep(1,N),
	b0=0.1,
	b1=0.1,
	b2=0.1,
	b3=0.1,
	b4=0.1,
	b5=0.1,
	b6=0.1,
	a0=0.1,
	a1=0.1,
	a2=0.1,
)

n.adapt = 1000
n.iter = 25000
n.chains = 1

jm1 = jags.model("GroupProject-JAGS-withdate.R", 
 data=bird_data,
 inits=inits,
 n.chains=n.chains,
 n.adapt=n.adapt)

load.module("dic")
zc1 = coda.samples(jm1,variable.names=c("p","z", "psi","b0", "b1", "b2", "b3", "b4", "b5", "b6", "a0", "a1", "a2"),
  n.iter=n.iter)



