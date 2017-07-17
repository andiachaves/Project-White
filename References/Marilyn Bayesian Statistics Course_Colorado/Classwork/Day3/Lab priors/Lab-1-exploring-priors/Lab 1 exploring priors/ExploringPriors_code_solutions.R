####
#### Read in Moth and LifeSpan data
####

Moth <- read.table("MothData.txt",header=TRUE)
Needles <- read.table("NeedleLifeSpan.txt",header=TRUE)

# Look at data:
#Moth$Placed
#Needles

# Total Placed:
tot.placed <- sum(Moth$Placed)
# Total Removed
tot.removed <- sum(Moth$Removed)
# Total not removed
tot.remain <- sum(Moth$Placed)- sum(Moth$Removed)

# Sum of life-span data
tot.life <- sum(Needles$LifeSpan)
# Sample size:
n.life <- length(Needles$LifeSpan)


##########  Mevin's approach to Problem 4    #############

####
####  MH: 4 a
####

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

####
####  MH: 4 b
####

fcn.4b <- function(ab){
  ans1=0.1
  ans2=0.025
  a=ab[1]
  b=ab[2]
  tmp1=a/(a+b)
  tmp2=sqrt((a*b)/(((a+b)^2)*(a+b+1)))
  tmp=(tmp1-ans1)^2+(tmp2-ans2)^2  #  Note: this is a squared error loss function
  tmp
}
 
optim(c(10,100),fcn.4b)$par # These are a and b

####
####  MH: 4 c 
####

fcn.4c <- function(ab){
  ans1=0.1
  ans2=0.5
  a=ab[1]
  b=ab[2]
  tmp1=qbeta(.025,a,b)
  tmp2=qbeta(.975,a,b)
  tmp=(tmp1-ans1)^2+(tmp2-ans2)^2  #  Note: this is a squared error loss function
  tmp
}

beta.4c.out=optim(c(1,3),fcn.4c)$par # These are a and b
curve(dbeta(x,beta.4c.out[1],beta.4c.out[2]),lwd=2,xlim=c(0,1))
abline(v=qbeta(c(.025,.975),beta.4c.out[1],beta.4c.out[2]),col=8)
abline(v=.25,col=8)
beta.4c.out


########################################################################
##########
##########  Kiona's approach to Problem 4 (and the rest)   
##########
########################################################################

####
####  KO: 4 a
####	need package: rootSolve
#### 	to load: >library(rootSolve)
####

fcn.4a <- function(ab,val){
	a=ab[1]
	b=ab[2]
	F1 <- a+b - val["n"]
	F2 <- a/(a+b) - val["Ep"]
	c(F1=F1, F2=F2)
}

sol.4a <- multiroot(f=fcn.4a, c(0.045, 0.065),maxiter=100,rtol=1e-6,val=c(n=.1,Ep=.5))

####
####  KO: 4 b
####

fcn.4b <- function(ab,val){
	a=ab[1]
	b=ab[2]
	F1 <- a/(a+b) - val["Ep"]
	F2 <- sqrt(a*b/(((a+b)^2)*(a+b+1)))-val["SDp"]
	c(F1=F1, F2=F2)
}

sol.4b <- multiroot(f=fcn.4b, c(6, 12),maxiter=100,rtol=1e-6,positive=TRUE,val=c(Ep=0.1,SDp=0.025))


####
####  KO: 4 c
####

fcn.4c <- function(ab,val){
	a=ab[1]
	b=ab[2]
	F1 <- (a-1)/(a+b-2) - val["mode"]
	F2 <- pbeta(val["q"],a,b)-val["prob"]
	c(F1=F1, F2=F2)
}

sol.4c <- multiroot(f=fcn.4c, c(6, 12),maxiter=100,rtol=1e-6,positive=TRUE,val=c(mode=.25,q=.1,prob=0.01))


####
####  KO: 5
####
####	This is a little tricky because need to recognize that the mean life-span is given
#### 	by E(Y) = 1/beta given the rate parameterization of the exponential likelihood.
####	So, beta is the expected turn-over (1/yr), and if beta ~ gamma(a,b), then this is
#### 	equivalent to 1/beta ~ inv-gamma(a,b). Use moment matching with the inv-gamma
####	to solve for hyperparameters a and b (or, c and d in hand-out?).
####

fcn.5 <- function(ab,val){
	a=ab[1]
	b=ab[2]
	F1 <- b/(a-1) - val["Ep"]
	F2 <- sqrt((b^2)/(((a-1)^2)*(a-2)))-val["SDp"]
	c(F1=F1, F2=F2)
}

sol.5 <- multiroot(f=fcn.5, c(6, 12),maxiter=100,rtol=1e-6,positive=TRUE,val=c(Ep=12,SDp=2))


####
####  KO: 6
####

# Likehood (unnormalized) for moth data:
likeM <- function(x){
	(x^132)*((1-x)^352)
}

# Likehood (unnormalized) for needle life-span data:
likeL <- function(x){
	(x^18)*exp(-x*147.6)
}

par(mfrow=c(3,2))
# For prior M1:
curve(dbeta(x,0.05,0.05),from=0,to=1,lwd=3,ylim=c(0,20),n=5000,ylab="Density")
curve(dbeta(x,132+0.05,352+0.05),from=0,to=1,add=TRUE,col="cyan",lwd=3)
curve(likeM(x)*20/(likeM(132/(132+352))),from=0,to=1,add=TRUE,col="red",lwd=3)
#legend(x=0.5,y=18,legend=c("prior","posterior","likelihood"),lty=1,lwd=2,col=c("black","cyan","red"))

# For prior M2:
curve(dbeta(x,14,129),from=0,to=1,lwd=3,ylim=c(0,25),n=5000,ylab="Density")
curve(dbeta(x,132+14,352+129),from=0,to=1,add=TRUE,col="cyan",lwd=3)
curve(likeM(x)*23.5/(likeM(132/(132+352))),from=0,to=1,add=TRUE,col="red",lwd=3)
#legend(x=0.5,y=18,legend=c("prior","posterior","likelihood"),lty=1,lwd=2,col=c("black","cyan","red"))

# For prior M3:
curve(dbeta(x,6.9,18.6),from=0,to=1,lwd=3,ylim=c(0,22),n=5000,ylab="Density")
curve(dbeta(x,132+6.9,352+18.6),from=0,to=1,add=TRUE,col="cyan",lwd=3)
curve(likeM(x)*20/(likeM(132/(132+352))),from=0,to=1,add=TRUE,col="red",lwd=3)
#legend(x=0.5,y=18,legend=c("prior","posterior","likelihood"),lty=1,lwd=2,col=c("black","cyan","red"))

# For prior L1; for beta (turn-over rate, yr-1):
curve(dgamma(x,38,444),from=0,to=0.5,lwd=3,ylim=c(0,35),n=5000,ylab="Density")
curve(dgamma(x,18+38,147.6+444),from=0,to=0.5,add=TRUE,col="cyan",lwd=3)
curve(likeL(x)*32/(likeL(18/147.6)),from=0,to=0.5,add=TRUE,col="red",lwd=3)
#legend(x=0.3,y=34,legend=c("prior","posterior","likelihood"),lty=1,lwd=2,col=c("black","cyan","red"))

# For prior L1; for 1/beta (life-span, yr):
curve(densigamma(x,38,444),from=2.5,to=21,lwd=3,ylim=c(0,0.3),xlim=c(0,25), n=5000,ylab="Density")
curve(densigamma(x,18+38,147.6+444),from=2.5,to=21,add=TRUE,col="cyan",lwd=3)
curve(likeL(x)*0.28/(likeL(18/147.6)),from=2.5,to=21,add=TRUE,col="red",lwd=3)

plot(1,2,type="n",axes=FALSE,xlab=NA,ylab=NA)
legend(x=1,y=2,xjust=0.5,yjust=0.5,legend=c("prior","posterior","likelihood"),lty=1,lwd=2,col=c("black","cyan","red"),cex=2)



####
####  KO: 7
####	need package: pscl
####  to load: > library(pscl)

hpd.fcn <- function(LS,ab,prob){
	a=ab[1]
	b=ab[2]
	F1 <- densigamma(LS[1],a,b)-densigamma(LS[2],a,b)
	F2 <- pigamma(LS[2],a,b)-pigamma(LS[1],a,b)-prob
	c(F1=F1, F2=F2)
}

hpd.prior <- multiroot(f=hpd.fcn, c(8, 16),maxiter=100,rtol=1e-6,positive=TRUE,ab=c(38,444),prob=c(0.95))
hpd.post <- multiroot(f=hpd.fcn, c(8, 14),maxiter=100,rtol=1e-6,positive=TRUE,ab=c(38+18,444+147.6),prob=c(0.95))

# Central 95% CI prior
qigamma(0.025,38,444)
qigamma(0.975,38,444)


# Central 95% CI posterior
qigamma(0.025,38+18,444+147.6)
qigamma(0.975,38+18,444+147.6)

####
####  KO: 8
####	

# Define function to generate unobserved or new data (y) from a binomial
# distribution given that the probability parameter (theta) arised from
# a beta distribution with parameters a and b (defaults set to the prior
# values, but can change these to the posterior values since the beta
# is conjugate to the likelihood (and thus, posterior for theta is also a beta).

pred.y <- function(Nsamples,a=prior.a, b=prior.b, n){
	theta <- rbeta(Nsamples,shape1=a, shape2=b)
	#y <- rbinom(Nsamples,size=rep(n,length(theta)),prob=theta)
	y <- rbinom(Nsamples,n,prob=theta)
}

# Generate y from its prior predictive and posterior predictive distributions
# and plot these for the 3 different priors (M1, M2, M3)

par(mfrow=c(3,2))
Nsamples = 5000
n=100
# Prior predictive: Prior M1
y.prior.pred <- pred.y(Nsamples ,a=sol.4a$root[1],b=sol.4a$root[2],n)
# Posterior predictive: Prior M1
y.post.pred <- pred.y(Nsamples ,a=sol.4a$root[1]+tot.removed,b=sol.4a$root[2]+tot.remain,n=100)
hist(y.prior.pred,xlim=c(0,100),breaks=(-1:100)+0.5,freq=FALSE,ylab="Density")
hist(y.post.pred,,xlim=c(0,100),breaks=(-1:100)+0.5,freq=FALSE,ylab="Density")

# Prior predictive: Prior M2
y.prior.pred <- pred.y(Nsamples ,a=sol.4b$root[1],b=sol.4b$root[2],n)
# Posterior predictive: Prior M2
y.post.pred <- pred.y(Nsamples ,a=sol.4b$root[1]+tot.removed,b=sol.4b$root[2]+tot.remain,n=100)
hist(y.prior.pred,xlim=c(0,100),breaks=(-1:100)+0.5,freq=FALSE,ylab="Density")
hist(y.post.pred,xlim=c(0,100),breaks=(-1:100)+0.5,freq=FALSE,ylab="Density")

# Prior predictive: Prior M3
y.prior.pred <- pred.y(Nsamples ,a=sol.4c$root[1],b=sol.4c$root[2],n)
# Posterior predictive: Prior M2
y.post.pred <- pred.y(Nsamples ,a=sol.4c$root[1]+tot.removed,b=sol.4c$root[2]+tot.remain,n=100)
hist(y.prior.pred,xlim=c(0,100),breaks=(-1:100)+0.5,freq=FALSE,ylab="Density")
hist(y.post.pred,xlim=c(0,100),breaks=(-1:100)+0.5,freq=FALSE,ylab="Density")

