#Light limitation of trees example of Gibbs sampler using Metropolis Hastings steps

rm(list=ls())
set.seed(4)
setwd("C:\\Users\\MBrandt\\Dropbox\\Transfer documents\\conferences\\Bayesian Statistics\\Classwork\\Day3\\Lab MCMC\\MCMC-Lab-11\\")

#Execute the function definitions.
source("Gibbs MH functions-mycode.R")

#Simulate data to resemble data use in likelihood lab
data=get_data(alpha=38.5, gamma=1.7, sigma=2, c=8)
# call "data" and it will give you all of the information generated.


#Set up the number of iterations. While you are getting things working, this can be small, ca 5000
n.iter=50000

#Set up storage list (x), assign initial values to parameters, and tunning valuesInitialize the chain
#dim.x gives the dimensions for the array of each estimated quantity.  For scalars = 1, for vectors > 1
#Initial values and tuning parameters need to be set in the function body
x=setup(n.iter=n.iter,n.chain=1, parameter.names=c("alpha","c","gamma","sigma", "y.hat","growth_ratio"), dim.x=c(1,1,1,1,50,1))
#Run the MCMC.  Function returns the filled chains
x=Run_MCMC(x=x,tune=x$tune, n.iter=n.iter)

#Set value for burin--the elements that will be discarded
burnin=15000
#Look at the output
do_plots(data=data, x=x, n.iter=n.iter, burnin=burnin)
	



