#Lab 2

##generate 800 values using the binomial distribution
t <- rbinom(800, 40, 0.20)

t
#plot the values in a histogram
hist(rbinom(800, 40, 0.20))

x <-rbinom(10000,prob=0.2,size=12)

##Tabulate the values and divide by the number of samples to get a probability distribution
tx <- table(factor(x,levels=0:12))/10000

## Creates a barplot of the values, that will extend the y limits to make room for the 
#theoretical values and saving the x locations at which the bars are drawn
b1 <- barplot(tx,ylim=c(0,0.3),ylab="Probability")

#Add the theoretical values, plotting them at the same x-locations as the centers 
#of the bars:
points(b1, dbinom(0:12, prob = 0.2, size = 12), pch = 16)

#######Exercise 1.2 - Pick 10,000 negative binomial deviates with 
#mean = 2, k = 0.5. Draw histogram. Check that the mean and variance
#agree reasonably well with the theoretical values. Add points that 
#represent the theoretical distribution to the plot.

#Neg binomial: rnbinom(n, size, prob, mu) Parameterized with prob&size OR 
#mu & k(size)


j <- rnbinom(10000,2,0.5)
tj <- table(factor(j,levels=0:12))/10000
b2 <- barplot(tj,ylim=c(0,0.3),ylab = "probability")
points(b2, dnbinom(0:12, mu = 2, size = 0.5), pch = 16)




##Exercise 2 - Creating new distributions
#define a probability distribution for a zero-inflated negative binomial: 

#density function for zero-inflated negative binomial:
dzinbinom = function(x, mu, size, zprob)
{
ifelse(x == 0, 
zprob + (1 - zprob) * dnbinom(0, mu = mu, size = size), 
(1 - zprob) * dnbinom(x, mu = mu, size = size))
}

#random generator for the zero-inflated negative binomial:
rzinbinom = function(n, mu, size, zprob) 
{
ifelse(runif(n) < zprob,
0,
rnbinom(n, mu = mu, size = size))
}


m <- rzinbinom(10000,2,0.5,0.5)
tm <- table(factor(m,levels=0:12))/10000
b3 <- barplot(tm,ylim=c(0,1),ylab = "probability")
points(b3, dzinbinom(0:12, mu = 2, size = 0.5, zprob = 0.5), pch = 16)


######Now write a density function and random deviate generator for:
#A zero-inflated Poisson distribution

#density function for zero-inflated Poisoon distribution:
#Need to know Poisson: rpois(n, lambda) n = number of rando draws, lambda is 
#expected value of the distribution


dzinpoisson = function(n, lambda, zprob)
{
ifelse(n == 0, 
zprob + (1 - zprob) * dpois(0, lambda = lambda),
(1 - zprob) * dpois(n, lambda = lambda))
}

rzinpoisson = function(n, lambda, zprob)
{
ifelse(runif(n) < zprob,
0,
rpois(n, lambda = lambda))
}

s <- rzinpoisson(10000,1.5,0.5)
ts <- table(factor(s,levels=0:12))/10000
b4 <- barplot(ts,ylim=c(0,1),ylab = "probability")
points(b4, dzinpoisson(0:12, lambda = 1.5, zprob = 0.5), pch = 16)



####Compound distributions
#In R, can generate random deviates using vector of different parameters instead
#of just designating one parameter

#Typical generation of random deviates, where sample 8 clutches, each of 10 hatchlings, with 
#prob of survival of 0.8 (80%): 

rbinom(8, size = 10, prob = 0.8)

#But if clutches were all of different sizes, can pick random values using a vector 
#of sizes: 

clutch_size = c(10,9,9,12,10,10,8,11)
rbinom(8, size = clutch_size, prob = 0.8)


##3. Determining the error structure in your data:

#Read data in 

setwd("C:\\Users\\MBrandt\\Dropbox\\Transfer documents\\conferences\\Bayesian Statistics\\Classwork\\Day1-Intro\\Lab 2\\")

tows<-read.table("HMtab43.txt", header=TRUE)
