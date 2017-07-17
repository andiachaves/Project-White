##Bayes lab

#1. Simulate 50 data points from a Poisson distribution with mean = 6.4


set.seed(3)

lambda.data <- 6.4

x <- rpois(50, lambda.data)
plot(table(factor(x, levels = 0:max(x))), ylab = "Frequency", xlab = "x")


#2 Set up a vector containing a sequence of values for theta, the mean number of invasive
#plants, using code like theta = seq(0,15,0.01)

theta = seq(0,15,0.01)

#3 Write a function for the prior

mu.prior <- 10.2
sigma.prior <- 1.0


prior <- function(theta, mu=mu.prior, sigma=sigma.prior)
{
alpha = mu^2/sigma^2
beta = mu/sigma^2
dgamma(theta, shape=alpha, rate=beta)
}

prior.data <- prior(theta = theta)
plot(prior.data, typ = "l")

#4 - Plot histogram
hist(prior.data, freq=FALSE, breaks = 10)

#5 - Create a function for the likelihood. The function must use all 50 observations
#to output the total likelihood (not the log likelihood) of the data for each value of theta.

likelihood <- function (theta, data = x) 
{
L=rep(0,length(theta))
for(i in 1:length(theta))
L[i] = prod(dpois(data, theta[i], log=FALSE))
return(L)
}

plot(theta,likelihood(theta), typ="l", xlab=expression(theta), ylab=expression(paste("P(y | ",theta,")")), main="Likelihood",xlim=c(0,15))


#6 Create a function for the joint distribution of the parameters and the data as the product of the prior and the likelihood function.

joint <- function(theta) 
{
likelihood(theta)*prior(theta)
}

#make the plot a line plot because it's continuous
plot(theta,joint(theta), typ="l")


#7 - Need to integrate the joint distribution using R's numerical integration function
#to obtain the normalization constant, P(y). 

Py = integrate(joint, 0, 30, abs.tol=1e-100)

Py

#8 Now estimate the posterior distribution by dividing each element of the vector of 
#the output produced by the joint function by the integratl of the joint function.

P.theta <- joint(theta)/Py$value

plot(theta,P.theta,typ="l")

par(mfrow=c(2,3))
plot(theta,prior(theta),typ="l")
hist(x, freq=FALSE, breaks=10)
plot(table(factor(x, levels = 0:max(x))), ylab = "Frequency", xlab = "x")
plot(theta,likelihood(theta),typ="l")
plot(theta,joint(theta),typ="l")
plot(theta,P.theta,typ="l")

#The posterior distribution goes above one because the area has to integrate to 1.

#=======scaled, overlay==============================
#Scale the liklelihood so it can be compared to prior and posterior on overlay plot
like.v=likelihood(theta)
like.scaled = like.v/max(like.v)*max(P.theta)
#pdf("Scaled_overlay.pdf", height=4, width=6)
par(mfrow=c(1,1))

plot(theta,like.scaled,typ="l",col="red",xlim=c(5,15),xlab=expression(theta),ylab= "Probability density", main="Scaled Overlay")
lines(theta,prior(theta),col="blue")

#Calculate the posterior using gamma-Poisson conjugate relationship and overlay on the integrated posterior
P.conj=dgamma(theta,sum(x)+mu.prior^2/sigma.prior^2, length(x)+mu.prior/sigma.prior^2)
lines(theta,P.conj,typ ="l",lwd=1,col="orange")
lines(theta,P.theta,col="black", lwd=4, typ="l", lty="dashed")

legend(11,1, legend= c("Scaled Likelihood", "Prior","Integrated posterior", "Conjugate posterior"), cex=.8,lwd=2, bty="n", col=c("red","blue","black","orange"), lty=c("solid","solid","dashed","solid"))



