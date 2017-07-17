# Components of Bayes theorem
rm(list=ls())
setwd("/Users/Tom/Documents/NSF Statistics Workshop/2013 Course/Hobbs course materials/Bayes lecture and lab/")
par(cex=.8, lwd=2,cex.lab=1.5,cex.main=1.5,mai=c(1.5,1,1,1), mfrow=c(2,3))
#set random number generator
set.seed(3)
#vector for theta:
theta=seq(0,15,.01)
#simulate the data
data=rpois(50, lambda=6.4)
#prior mean and standard deviation
mu.prior=10.2;sigma.prior=0.5


#=======function returning the probabilty density of theta based on prior knowledge==================
Prior=function(theta, mu=mu.prior, sigma=sigma.prior) dgamma(theta,mu^2/sigma^2,mu/sigma^2)
#check moment matching in prior:
sd(rgamma(100000,mu.prior^2/sigma.prior^2,mu.prior/sigma.prior^2))
mean(rgamma(100000,mu.prior^2/sigma.prior^2,mu.prior/sigma.prior^2))
#pdf(file="Bayes collage.pdf",width=7, height = 5)
par(mfrow=c(2,3))
plot(theta,Prior(theta), typ="l", ylab=expression(paste("P(",theta,")")), xlab=expression(theta),main="Prior", xlim=c(5,15))
#==========================================================================


#======Histogram of the data============================
hist(data,freq=FALSE,breaks=10)
discrete_hist=function(data) {
	library(plyr)
	y=count(data)
	dens=y$freq/length(data)
	y=cbind(y,dens)
	plot(y$x,y$dens,typ="h", ylab="Density",xlab="y",main="Improved histogram of data" ,frame=FALSE)
#lines(y$x,dpois(y$x,mu), col="red")
#plot(y$x,y$freq,typ="h",,ylab="Frequency",xlab="y", main=paste("n=",n), frame=FALSE)
}
discrete_hist(data)
#==================================================


#==========function for the likelihood: y is vector of data, theta a vector of parameter values================
Like=function(theta,y=data){
	L=rep(0,length(theta))
	for(i in 1:length(theta))	L[i] = prod(dpois(y,theta[i], log=FALSE))
	return(L)
	} #end of Likelihood function function
plot(theta,Like(theta), typ="l", xlab=expression(theta), ylab=expression(paste("P(y | ",theta,")")), main="Likelihood",xlim=c(5,15))
#============================================================================

#======function for the joint distribution of the parameter and the data====================
Joint=function(theta) Like(theta)*Prior(theta)
plot(theta,Joint(theta),typ="l",ylab=expression(paste("P(y | ",theta,") x P(",theta,")"   )), main="Joint", xlab=expression(theta), xlim=c(5,15) )
#===================================================================

#=======marginal distribution: integrate the joint distribtuion over theta======================
Py=integrate(Joint,0,30,abs.tol=1e-100)
Py  #Note that Py is a very small scalar
#========================================================================

#The posterior distribution==================
P.theta = Joint(theta)/Py$value 

plot(theta,P.theta,typ="l", xlab=expression(theta), ylab=expression(paste("P( ",theta," | y)")),xlim=c(5,15), main="Posterior")
#dev.off()
#===========================================

#=======scaled, overlay==============================
#Scale the liklelihood so it can be compared to prior and posterior on overlay plot
like.v=Like(theta)
like.scaled = like.v/max(like.v)*max(P.theta)
#pdf("Scaled_overlay.pdf", height=4, width=6)
par(mfrow=c(1,1))

plot(theta,like.scaled,typ="l",col="red",xlim=c(5,15),xlab=expression(theta),ylab= "Probability density", main="Scaled Overlay")
lines(theta,Prior(theta),col="blue")

#Calculate the posterior using gamma-Poisson conjugate relationship and overlay on the integrated posterior
P.conj=dgamma(theta,sum(data)+mu.prior^2/sigma.prior^2, length(data)+mu.prior/sigma.prior^2)
lines(theta,P.conj,typ ="l",lwd=1,col="orange")
lines(theta,P.theta,col="black", lwd=4, typ="l", lty="dashed")

legend(11,1, legend= c("Scaled Likelihood", "Prior","Integrated posterior", "Conjugate posterior"), cex=.8,lwd=2, bty="n", col=c("red","blue","black","orange"), lty=c("solid","solid","dashed","solid"))
#dev.off()
#======================================

