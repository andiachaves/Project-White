#Example 1: Seven sites, 52-92 dead moths were glued to tree trunks. After 24 hrs
# the number of moths that had been removed were counted. 
#Example 2: Want to estimate life-span (years) or needles of lodgepole pine. We
# have data on needle life-spans recorded for multiple (18) trees from multiple
# sites. 


#Need rootSolve and pscl packages

#Problems:
########1. Pick appropriate sampling distributions to describe the likelihoods associated
#with the dataset

#**Moth: binomial - Number removed(i) is distributed (~) binomially(N = place(i), prob of removal=theta)
#**Needle Lifespan(LS): Lifespan(i) ~ Exponential(alpha) E(Lifespan) = 1/alpha

########2. Identify conjugate priors for each example.

#**Moth: P(theta) = beta(a, b)
#**Needle Lifespan(LS): P(alpha) = gamma(c,d)

########3. Solve for the posterior for each example, in terms of the actual data

#**Moths: P(theta | removed) = beta (_ , _)
#		so = beta (a + sum of yi, b + sum of (ni - yi))  with sum of yi being sum of successes and sum of ni - yi being sum of failures

#**Lifespan: P(alpha | LS) = gamma (_, _)
#		so = gamma (N + c, sum of yi + d), where c is prior sample size and d is prior lifespan


setwd("C:\\Users\\MBrandt\\Dropbox\\Transfer documents\\conferences\\Bayesian Statistics\\Classwork\\Day3\\Lab x\\")

#########4. specify 3 different prior that will be used in the moth predation example:
#a. Prior M1: A relatively non-informative prior that contributes information equivalent to n = 0.1
# placed moths and an expected prior probability of predation of 0.5


moth <- read.table("MothData1.txt", header=TRUE)
#Sum Placed and Removed

sumN <- sum(moth$Placed)  	#sum of number of moths placed
sumY <- sum(moth$Removed)	#sum of number of moths removed

library(rootSolve)

#Define the functions for the "total prior sample size" (F1)
# and the prior mean for p (probability of predation/removal). 

##Kiona's approach


#Mevin's approach:
fcn.4a <- function(ab){
  ans1=0.1
  ans2=0.5
  a=ab[1]
  b=ab[2]
  tmp1=a+b
  tmp2=a/(a+b)
  tmp=(tmp1-ans1)^2+(tmp2-ans2)^2  #  Note: this is a squared error loss function
  tmp
}
 
optim(c(1,1),fcn.4a)$par # These are a and b


#what are the values to use for beta (a,b)?
#a is the prior number of moths removed







