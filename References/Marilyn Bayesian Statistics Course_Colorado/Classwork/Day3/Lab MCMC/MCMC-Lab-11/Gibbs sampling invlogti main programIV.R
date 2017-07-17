#Light limitation of trees example of Gibbs sampler using Metropolis Hastings steps

rm(list=ls())
library(boot)
set.seed(4)
setwd("/Users/Tom/Documents/Ecological Modeling Course/2013 Spring Semester/Labs/Lab 5 MCMC/2013 modifications/")

#Execute the function definitions.
source("Gibbs MH functions invlogitIV.R")
 
 #set sample size
 n=50
#simulate data
data=get_data(alpha=150, gamma=-11, sigma=5, kappa=.005, n=n)

#Set up the number of iterations. While you are getting things working, this can be small, ca 5000
n.iter=50000

#Set up storage list (x), assign initial values to parameters, and tunning values.  Initialize the chain
#dim.x gives the dimensions for the array of each estimated quantity.  For scalars = 1, for vectors > 1
#Initial values and tuning parameters need to be set in the function body
x=setup(n.iter=n.iter,n.chain=1, parameter.names=c("alpha","kappa","gamma","sigma", "y.hat"), dim.x=c(1,1,1,1,n))
#Run the MCMC.  Function returns the filled chains
x=Run_MCMC(x=x,tune=x$tune, support=x$support, n.iter=n.iter)

#Set value for burin--the elements that will be discarded
burnin=15000
#Look at the output
do_plots(data=data, x=x, n.iter=n.iter, burnin=burnin)
	



