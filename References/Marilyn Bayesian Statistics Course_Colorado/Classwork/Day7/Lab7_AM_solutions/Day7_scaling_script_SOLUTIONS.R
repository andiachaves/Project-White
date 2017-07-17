#### Load rjags and xtable libraries
library(rjags)
library(xtable)

#### Read in cavitation data:
data=(read.csv("MetabolicScalingData.csv"))
Nobs = length(data[,1])
Y = matrix(0,Nobs,2)
Y[,1] = data[,4]
Y[,2] = data[,3]
LogD = data[,5]
SP = data[,1]

#### Define data to be used in model, including "sample/group" sizes
scaling_data = list(N=256, Nsp=14, Nvars = 2, R=structure(.Data=c(1,0,0,1),.Dim=c(2,2)),
	LogD = LogD, Y = Y, SP = SP)

# The initial values for Model 1 are contained in the R file names "Model1.inits"
# and the list of initials are named inits1 for model 1. The source function will
# put inits1 in the workspace, to be used later.
source("Model1_inits.R")
n.chains = length(inits1)

#### Set-up the MCMC procedure for MODEL 1 by specifying the file that contains 
#### the model code, indicate the data, indicate the initials, indicate the number 
#### of chains, specify the length of the adapting phase.
#####################################################################
#### NOTE: need to specify the file name that contains the JAGS code
#####################################################################

#### Set n.adapt (i.e., number of initial iterations for adapting
#### phase, if JAGS uses, e.g., Metropolis-Hastings).
n.adapt=1000

jm1=jags.model("Day7_model1_SOLUTIONS.R",
               data=scaling_data,
               inits=inits1,
               n.chains=n.chains,
               n.adapt=n.adapt)

#### Run the MCMC chains and store the sample obtained for those variables.names
#### specified. Need to monitor quantities necessary for evaluating model fit and 
#### for comparing between models (i.e., need to monitor Yrep and Dsum in addition 
#### to deviance). Also need to load the dic module so we can monitor deviance and
#### compute DIC. Set initial number of n.iter to evaluate length of burn-in period:
load.module("dic")
n.iter=2000
zc1 = coda.samples(jm1,variable.names=c("deviance","Sigma","alpha","beta","mu.alpha",
                                        "mu.beta","sig.alpha","sig.beta","rho","sig",
                                        "Yrep","Dsum"),
                   n.iter=n.iter)

#### Now view trace plots to evaluate burn-in and convergence.
#### With below command, R will ask you to hit enter to show each page of 
#### output/graphics:
devAskNewPage(ask = TRUE)

#### For ease of coding, create variables that indicate the start and end
#### iterations for the MCMC chains associated with the zc1 coda object:
st=start(zc1)
en=end(zc1) 
th=thin(zc1)

#### View trace plots of monitored quantities, but you probably don't want to
#### create plots of Yrep (too many!!!):
xyplot(window(zc1, start=st, end=en, thin=1)[,c("deviance")])
xyplot(window(zc1, start=st, end=en, thin=1)[,c("mu.alpha[1]","mu.alpha[2]",
                                                "mu.beta[1]","mu.beta[2]",
                                                "sig.alpha[1]","sig.alpha[2]",
                                                "sig.beta[1]","sig.beta[2]")], layout=c(2,4))
xyplot(window(zc1, start=st, end=en, thin=1)[,paste("alpha[",1:14,",",1:2,"]",sep="")],layout=c(2,4))
xyplot(window(zc1, start=st, end=en, thin=1)[,paste("beta[",1:14,",",1:2,"]",sep="")],layout=c(2,4))


#### Update the JAGS model (jm1) for enough iterations to accomodate the burn-in
#### period. E.g., set n.update such that it is a bit longer than the burn-in 
#### period, then update the JAGS model object:
n.update=1000
update(jm1, n.iter=n.update)


#### Rerun the JAGS model via coda.samples to update the coda object. Run for sufficiently
#### long, with the goal of using these coda samples to compute posterior stats:
n.iter = 5000
zc1 = coda.samples(jm1,variable.names=c("deviance","Sigma","alpha","beta","mu.alpha",
                                        "mu.beta","sig.alpha","sig.beta","rho","sig",
                                        "Yrep","Dsum"),
                   n.iter=n.iter)

#### Recreate trace plots to verify that we've passed the burn-in:
st=start(zc1)
en=end(zc1) 
th=thin(zc1)
xyplot(window(zc1, start=st, end=en, thin=1)[,c("deviance")])
xyplot(window(zc1, start=st, end=en, thin=1)[,c("mu.alpha[1]","mu.alpha[2]",
                                                "mu.beta[1]","mu.beta[2]",
                                                "sig.alpha[1]","sig.alpha[2]",
                                                "sig.beta[1]","sig.beta[2]")], layout=c(2,4))
xyplot(window(zc1, start=st, end=en, thin=1)[,paste("alpha[",1:14,",",1:2,"]",sep="")],layout=c(2,4))
xyplot(window(zc1, start=st, end=en, thin=1)[,paste("beta[",1:14,",",1:2,"]",sep="")],layout=c(2,4))

#### Specify burn-in period, if any (could set bi=0 if burn-in has passed):
bi = 100

#### Look at within chain autocorrelation for deviance to determine if we need
#### to thin, especially since we may want to compute the variance of the 
#### deviance as part of our model selection criteria
autocorr.plot(window(zc1, start=st+bi, end=en, thin=1)[,c("deviance")],lag.max=30)
#### Re-plot autocorrelation function after thinning:
autocorr.plot(window(zc1, start=st+bi, end=en, thin=3)[,c("deviance")],lag.max=30)


#### Produce posterior summary statistics; discard burn-in and thin samples
#### Should thin some quantities by more than 10 due to excessive within
#### chain autocorrelation (thinning necessary for getting unbiased posterior
#### standard deviations or variances, but not necessary for getting posterior
#### means or quantiles). Thin by 10 here for demo purposes.
n.thin = 3
#### Replicated stats for the first response variable (logLength):
yrep1.stats=summary(window(zc1, start=st+bi, end=en, thin=n.thin)
               [,paste("Yrep[",1:scaling_data$N,",1]",sep="")])$stat
yrep1.quantiles=summary(window(zc1, start=st+bi, end=en, thin=n.thin)
               [,paste("Yrep[",1:scaling_data$N,",1]",sep="")])$quantiles
#### Replicated stats for the 2nd response variable (logMass):
yrep2.stats=summary(window(zc1, start=st+bi, end=en, thin=n.thin)
                    [,paste("Yrep[",1:scaling_data$N,",2]",sep="")])$stat
yrep2.quantiles=summary(window(zc1, start=st+bi, end=en, thin=n.thin)
                        [,paste("Yrep[",1:scaling_data$N,",2]",sep="")])$quantiles

#### Compute Posterior Predictive Loss (Dinf) and its components (G and P)
#### Compute for the 1st and 2nd response variables (logLength, LogMass)
G = matrix(0,1,2)
P = G
G[1] <- sum((yrep1.stats[,1]-Y[,1])^2)
G[2] <- sum((yrep2.stats[,1]-Y[,1])^2)
P[1] <- sum((yrep1.stats[,2])^2)
P[2] <- sum((yrep2.stats[,1])^2)
Dinf <- G + P

#### Call the dic.samples function, which operates on the updated JAGS model object
#### (it does not operate on the coda object, zc1). Run the dic function to obtain 
#### the DIC statistics. Run for sufficiently long.
dicj = dic.samples(jm1,n.iter=1000,thin=n.thin,type="pD")
dicj

#### As an alternative to the above DIC values, compute the mean deviance and
#### and pd* = (1/2)*Var(deviance)
dev.mean = summary(window(zc1, start=st+bi, end=en, thin=n.thin)[,"deviance"])$stat[1]
pd.star = 0.5*((summary(window(zc1, start=st+bi, end=en, thin=n.thin)[,"deviance"])$stat[2])^2)
dic.star = dev.mean + pd.star


#### Creat a scatter plot of predicted (y = posterior mean for replicated
#### data) versus the observed data (x = Y[,1]), and overlay a 1:1 line.
#### Create vertical "arrows" to represent the 95% CI for each replicated 
#### data point. The 2.5th and 97.5th percentiles are in the yrep1.quantiles
#### and yrep2.quantiles arrays
plot(Y[,1],yrep1.stats[,1],ylim=c(3,5),xlim=c(3,5),pch=20,cex=.75,
     ylab="Replicated Log(L) (mean and 95% CI)",xlab="Observed Log(L)")
arrows(x0=Y[,1], y0=yrep1.quantiles[,1], x1=Y[,1], y1=yrep1.quantiles[,5], code=3, angle=90, length=0,col="gray")
points(Y[,1],yrep1.stats[,1],pch=20,cex=.75)
abline(c(0,1))

plot(Y[,2],yrep2.stats[,1],ylim=c(1,7),xlim=c(1,7),pch=20,cex=.75,
     ylab="Replicated Log(M) (mean and 95% CI)",xlab="Observed Log(M)")
arrows(x0=Y[,2], y0=yrep2.quantiles[,1], x1=Y[,2], y1=yrep2.quantiles[,5], code=3, angle=90, length=0,col="gray")
points(Y[,2],yrep2.stats[,1],pch=20,cex=.75)
abline(c(0,1))