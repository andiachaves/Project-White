##5/21/13 - Lab 3 - Likelihood and maximum likelihood
# Part 1: will simulate a data set and calculate the parameters of the distribution
#using both MOM and ML. This is to demonstrate the MOM is typically biased and MLE are
# preferred.

#Load MLE libraries:
library(stats4)
library(bbmle)


mu.true <- 1
k.true <- 0.4
x <- rnbinom(50, mu=mu.true,size = k.true)

#Review: rnbinom(n, size, prob, mu)
#n - number of observations
# size - target for number of successful trials or dispersion parameter (must be positive)
# prob - prob of success in each trial
# mu - alternative parameterization via mean.

#For negative binomial can parameterize in two ways: 
# one with clustering K 


plot(table(factor(x, levels = 0:max(x))), ylab = "Frequency", xlab = "x")

#Build a function that calculates the neg log-likelihood for the distribution
#given a set of parameters. The arguments for the function are p, the vector of 
#parameters (mu and k) and dat the vector of data.

NLLfun1 = function(dat = x, mu, k )
{
-sum(dnbinom(x, mu=mu, size = k, log = TRUE))
}

##Calculate the neg log-likelihood with the true values of the distribution. 
#Have to combine values into a vector to be able to pass them to the NLL function:

nll.true <- NLLfun1(mu = mu.true, k = k.true)
nll.true

#The output value from the function is 74.63333

#Find the method of moments estimate of the parameters. Know that: 

m = mean(x)
v = var(x)
mu.mom = m #the mu by the method of moments is the mean
k.mom = m/(v/m - 1) #the k by the method of moments incorporates the variance

# So the negative log-likelihood estimate for the MOM parameters is: 

nll.mom <- NLLfun1(mu = mu.mom, k = k.mom)
nll.mom

#The answer given for this is 75.43553

## What is the difference in likelihood of the two estimates?  The LRT test would 
# say that it has to be greater than a chi-square with two degrees of freedom (0.95)/2: 

Ldiff <- nll.true - nll.mom
Ldiff
qchisq(0.95, df = 2)/2 

#The Ldiff is not greater than the chi-square, so there is no significant difference between
# what's produced based on the MOM estimates and that based on the "true" parameters

###Now testing the MLE estimates. 
#Use mle2 with the default Nelder-Mead algorithm and use the MOM estimates as the 
#starting conditions.

sol1 = mle2(minuslogl = NLLfun1, start = list(mu = mu.mom, k = k.mom), hessian = TRUE)
summary(sol1)

#mle2 arguments: 
#minuslogl: Function to calculate negative log-likelihood, or a formula
#start: Named list. Initial values for optimizer (so this has mu and k)
#hessian: options for Hessian calculation, passed through to the hessian function)

# so estimates given for mu and k are: 
# mu = 1.24013
# k = 0.38809
# -2 log L = 148.5354 - so divide this by 2 = 74.2677 (compare with true: 74.63333, and
# with MOM: 75.43553) it is smaller so does a better job of estimating.

#Find likelihood surfaces, profiles and confidence intervals

confint(sol1) 	#confint computes confidence intervals for one or more parameters in a fitted model.
			#There is a default and a method for objects inheriting from class "lm"
plot(profile(sol1))




##Exercise 1:
#Generate data of your choice using one of the probability functions you learned during
# the last lab. Calculate the MOM estimates of the parameters. Calculate the MLE and
# estimate the difference in likelihood and in prediction. Calculate the support intervals
# for the parameters. 


#Will try for Poisson
lambda.true <- 1.7
x <- rpois(100,lambda = lambda.true)

plot(table(factor(x, levels = 0:max(x))), ylab = "Frequency", xlab = "x")

NLLfun2 = function(dat = x, lambda) 
{
-sum(dpois(x, lambda = lambda.true, log = TRUE))
}

nll.true <- NLLfun2(lambda = lambda.true)
nll.true

# Returned: 151.9238

#Find MOM estimate of parameter
# For Poisson lambda = mean and lambda = var

m = mean(x)
lambda.mom = m

#And neg log-likelihood for the method of moments parameter: 

nll.mom <- NLLfun2(lambda = lambda.mom)
nll.mom

# Returned: 151.9238
# Hmmmm.... no difference at all....

solpois = mle2(minuslogl = NLLfun2, start = list(lambda = lambda.mom), hessian = TRUE)
summary(solpois)

# Returned 303.8475 - divide by 2 and get 151.92375 - smaller by a very small amount

# ****THEY are all very close because we're dealing with simple distributions****


####Exercise 2

#Read in hemlock data from Excel file


#Converted excel file to a csv file and then read in data:

setwd("C:\\Users\\MBrandt\\Dropbox\\Transfer documents\\conferences\\Bayesian Statistics\\Classwork\\Day2\\")
data <- read.csv("Hemlock-light-data-forR.csv")

attach(data)


#Check results of Solver by using nls() function in R 
# nls is nonlinear least squares 

###  The process model is: a*(Light-c)/((a/b)+(Light-c))
# Input the process model directly into the nls function:

##starting points for the parameters - mine were 50, 1 and 5
start <- list (a=50, b=1, c=5)
m1 <- nls(ObsGrowth~a*(Light-c)/((a/b)+(Light-c)), start = start)

#plot the data
plot(Light,ObsGrowth)

#plot the line using estimated coefficients from the nls
a<-coef(m1)["a"]
b<-coef(m1)["b"]
c<-coef(m1)["c"]
x<-Light
curve(a*(x-c)/((a/b)+(x-c)), add=TRUE, col="blue")

pred <- fitted(m1) 	#fitted is a generic function that extracts fitted values from 
				#objects returned by modeling functions

plot(ObsGrowth, pred, xlim=c(0,50), ylim=c(0,50))
abline(0,1)

# Just for fun....you might also try

pred.growth<-function(a,b,c){a*(Light-c)/((a/b)+(x-c))}


NLLfun3<-function(a,b,c,sd)
	{
	pred<-pred.growth(a,b,c)
	-sum(dnorm(ObsGrowth,mean=pred,sd=sd,log=TRUE),na.rm=TRUE)
	}


m6<-mle2(NLLfun3,start=list(a=30,b=2, c=6,sd=sd(ObsGrowth)),method="L-BFGS-B",lower=0.01)

summary(m6)
a<-coef(m6)["a"]
b<-coef(m6)["b"]
c<-coef(m6)["c"]
plot(Light,ObsGrowth)
pred.grow<-pred.growth(a,b,c)
x<-Light
curve(a*(x-c)/((a/b)+(x-c)), add=TRUE, col="blue")
plot(ObsGrowth, pred.grow, xlim=c(0,50), ylim=c(0,50))
abline(0,1)

###In Maria's code: use the gamma function 

#Reef fish data


setwd("C:\\Users\\MBrandt\\Dropbox\\Transfer documents\\conferences\\Bayesian Statistics\\Classwork\\Day2\\")
data <- read.table("Reef_fish3.txt")
attach(data)

plot(settlers, recruits)






































