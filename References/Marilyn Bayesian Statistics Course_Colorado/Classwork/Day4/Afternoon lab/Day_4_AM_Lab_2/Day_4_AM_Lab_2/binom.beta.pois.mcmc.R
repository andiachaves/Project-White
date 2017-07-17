binom.beta.pois.mcmc <- function(y,N.total,n.mcmc){

####
####  Mevin Hooten 
####
#######################
####  This code implements the fish transfer/survival model that Royle and Dorazio (2008) 
####  discuss in their book (Example 2.4.4.2) and pg. 59. 
####
####  y: vector of data
####  N.total: total number of fish before transfer
####  N: vector of fish transferred to each container (auxiliary variable).
####
#######################
####
####
#######################
####  Last Updated:
####
####  20110824: Created template for later use. 
####  20110831: Updated for class example.
####
####
#######################

####
####  Subroutines and Libraries
####



####
####  Setup Variables 
####

n=length(y)

N.save=matrix(0,n,n.mcmc)
phi.save=rep(0,n.mcmc)

####
####  Priors and Starting Values 
####

lambda=N.total/n

N=y+1
N.save[,1]=N

####
####  Begin Gibbs Loop 
####
  
for(k in 2:n.mcmc){
  cat(k," ")

  ####
  ####  Sample phi 
  ####

  phi=rbeta(1,sum(y)+1,sum(N-y)+1)

  ####
  ####  Sample N 
  ####

  N=y+rpois(n,(1-phi)*lambda)

  ####
  ####  Save Samples 
  ####

  phi.save[k]=phi
  N.save[,k]=N

}
cat("\n")

####
####  Write Output 
####
 
list(n.mcmc=n.mcmc,phi.save=phi.save,N.save=N.save)

}
