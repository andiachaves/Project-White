####
#### Read in Moth and LifeSpan data
####

Moth <- read.table("MothData.txt",header=TRUE)
Needles <- read.table("NeedleLifeSpan.txt",header=TRUE)

# Look at data:
Moth$Placed
Moth$Removed
Needles

# Total of Placed:
sum(Moth$Placed)
# Total of Removed
sum(Moth$Removed)
# Total of not removed
sum(Moth$Placed)- sum(Moth$Removed)

# Sum of life-span data
sum(Needles$LifeSpan)
# Sample of life-span size:
length(Needles$LifeSpan)

####
####  Problem 4a
####	get package: rootSolve
#### 	to load: >library(rootSolve)
####


# Define the functions for the "total prior sample size" (F1)
# and the prior mean for p (probability of predation/removal).
fcn.4a <- function(ab,val){
	a=ab[1]
	b=ab[2]
	F1 <- a+b - val["n"]
	F2 <- a/(a+b) - val["Ep"]
	c(F1=F1, F2=F2)
}

# Use multiroot function to solve for the values of a and b that
# satisfy the 2 equations for "prior sample size" and prior mean for p:
sol.4a <- multiroot(f=fcn.4a, c(0.045, 0.065),maxiter=100,rtol=1e-6,val=c(n=.1,Ep=.5))


####
####  Problem 5
####
####	This is a little tricky because need to recognize that the mean life-span is given
#### 	by E(Y) = 1/beta given the rate parameterization of the exponential likelihood.
####	So, beta is the expected turn-over (1/yr), and if beta ~ gamma(a,b), then this is
#### 	equivalent to 1/beta ~ inv-gamma(a,b). Use moment matching with the inv-gamma
####	to solve for hyperparameters a and b (or, c and d in hand-out?).
####

# Define the functions for the prior expected life-span (F1)
# and the prior standard deviation of life-span. Fill-in the ????

fcn.5 <- function(ab,val){
	a=ab[1]
	b=ab[2]
	F1 <- b/(a-1) - val["Ep"]
	F2 <- ???? - val["SDp"]
	c(F1=F1, F2=F2)
}

# Use multiroot to solve for a and b:
sol.5 <- multiroot(???)


####
####  Problem 6
####	
####	need package: pscl (for inverse-gamma pdf)
####  to load: > library(pscl)
####

# Define likehood (unnormalized) for moth data:
likeM <- function(x){
	(x^132)*((1-x)^352)
}

# Define likehood (unnormalized) for needle life-span data:
likeL <- function(x){
	(x^18)*exp(-x*147.6)
}

# Use par function to create a panel of 3 x 2 figures:
par(mfrow=c(3,2))

# Plots for prior M1:
curve(dbeta(x,0.05,0.05),from=0,to=1,lwd=3,ylim=c(0,20),n=5000,ylab="Density")
curve(dbeta(x,132+0.05,352+0.05),from=0,to=1,add=TRUE,col="cyan",lwd=3)
curve(likeM(x)*20/(likeM(132/(132+352))),from=0,to=1,add=TRUE,col="red",lwd=3)

# Plots for prior M2:
curve(dbeta(x,14,129),from=0,to=1,lwd=3,ylim=c(0,25),n=5000,ylab="Density")
curve(dbeta(x,132+14,352+129),from=0,to=1,add=TRUE,col="cyan",lwd=3)
curve(likeM(x)*23.5/(likeM(132/(132+352))),from=0,to=1,add=TRUE,col="red",lwd=3)
#legend(x=0.5,y=18,legend=c("prior","posterior","likelihood"),lty=1,lwd=2,col=c("black","cyan","red"))

# Plots for prior M3: FILL-IN
curve(????)
curve(????)
curve(????)

# Plots for prior L1; for beta (turn-over rate, yr-1):
curve(dgamma(x,38,444),from=0,to=0.5,lwd=3,ylim=c(0,35),n=5000,ylab="Density")
curve(????)
curve(????)

# For prior L1; for 1/beta (life-span, yr):
curve(densigamma(x,38,444),from=2.5,to=21,lwd=3,ylim=c(0,0.3),xlim=c(0,25), n=5000,ylab="Density")
curve(????)
curve(????)

# Use the 6th panel to insert a legend:
plot(1,2,type="n",axes=FALSE,xlab=NA,ylab=NA)
legend(x=1,y=2,xjust=0.5,yjust=0.5,legend=c("prior","posterior","likelihood"),lty=1,lwd=2,col=c("black","cyan","red"),cex=2)


####
####  Problem 7
####	

# Define function to solve for HPD interval; need 2 equations to solve
# for 2 unknowns (i.e., lower and upper limits of the HPD inverval, LS[1] 
# and LS[2], respectively)

hpd.fcn <- function(LS,ab,prob){
	a=ab[1]
	b=ab[2]
	F1 <- densigamma(LS[1],a,b)-densigamma(LS[2],a,b)
	F2 <- ????
	c(F1=F1, F2=F2)
}

# Compute 95% "HPD" interval based on prior:
hpd.prior <- multiroot(f=hpd.fcn, c(8, 16),maxiter=100,rtol=1e-6,positive=TRUE,ab=c(38,444),prob=c(0.95))

# Compute 95% "HPD" interval based on posterior:
hpd.post <- multiroot(????)

# Central 95% CI prior
qigamma(0.025,38,444)
qigamma(0.975,38,444)

# Central 95% CI posterior
qigamma(????)
qigamma(????)