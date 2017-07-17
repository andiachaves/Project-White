model
{
sigma.proc ~ unif(0,5)
tau.proc <- 1/sigma.proc^2

lambda ~ dunif(0.5, 10)
phi ~ beta(a, b) #a and b are the shape parameters that get from moment matching
