#####
##### Source MCMC function for N-mixture
#####

source("binom.beta.pois.mcmc.R")

#####
##### Fit model with R&D data 
#####

mcmc.out=binom.beta.pois.mcmc(c(15,11,12,5,12),250,50000)

#####
##### View Trace Plots 
#####

matplot(t(mcmc.out$N.save),type="l",lty=1)

plot(mcmc.out$phi.save,type="l")

#####
##### View Posterior Results  
#####

mean(mcmc.out$phi.save[-(1:1000)])
quantile(mcmc.out$phi.save[-(1:1000)],c(0.025,0.975))

apply(mcmc.out$N.save[,-(1:1000)],1,mean)
apply(mcmc.out$N.save[,-(1:1000)],1,quantile,c(0.025,0.975))

hist(mcmc.out$phi.save[-(1:1000)],breaks=30,col=8)

#####
##### Source MCMC function for N-mixture with Multinomial 
#####

source("binom.beta.MN.mcmc.R")

#####
##### Fit model with R&D data 
#####

mcmc.MN.out=binom.beta.MN.mcmc(c(15,11,12,5,12),250,50000)

#####
##### View Trace Plots 
#####

matplot(t(mcmc.MN.out$N.save),type="l",lty=1)

plot(mcmc.MN.out$phi.save,type="l")

mcmc.MN.out$mh.N/mcmc.MN.out$n.mcmc

#####
##### View Posterior Results  
#####

mean(mcmc.MN.out$phi.save[-(1:1000)])
quantile(mcmc.MN.out$phi.save[-(1:1000)],c(0.025,0.975))

apply(mcmc.MN.out$N.save[,-(1:1000)],1,mean)
apply(mcmc.MN.out$N.save[,-(1:1000)],1,quantile,c(0.025,0.975))

hist(mcmc.MN.out$phi.save[-(1:1000)],breaks=30,col=8,prob=TRUE,xlim=c(0,1))
lines(density(mcmc.out$phi.save[-(1:1000)],n=1000),col=2,lwd=2)


