##Exercise: Plot your priors

#For K
K <- seq(1,1200,1)
K.prior <- dgamma (K, 0.001, 0.001)
plot(K,K.prior, type="l")
#For r
r <- seq(0, 0.2, 0.0001)
r.prior <- dgamma (r, 0.001, 0.001)
plot(r,r.prior, type="l")
#For tau
tau <- seq(1,2500,1)
tau.prior <- dgamma (tau, 0.001, 0.001)
plot(tau, tau.prior, type="l")


#############################
setwd("C:\\Users\\MBrandt\\Dropbox\\Transfer documents\\conferences\\Bayesian Statistics\\Classwork\\Day4\\JAGS\\Logistic-R-code\\")

rm(list=ls())
pop.data=(read.csv("Logistic-Data-II1.csv"))
names(pop.data)=c("Year","Population Size", "Growth Rate")

inits=list(
list(K=1500, r=.2, tau=2500),# , #chain 1
list(K=1000, r=.15, tau=1000))#, #chain 2
#list(K=900, r=.3, tau=500) #chain 3
#)  

n.xy = nrow(pop.data)

data=list(n=n.xy, x=as.double(pop.data$"Population Size"),y=as.double(pop.data$"Growth Rate"))
library(rjags)

n.adapt = 5000  ##the number of iterations that JAGS will use to choose the sampler
n.update = 10000  ##the number of iterations that will be discarded (i.e., the burn in)
n.iter = 250000 ##the number of iterations that will be stored in the chain

jm=jags.model("Logistic example BUGS.R",data=data, inits,n.chains=length(inits), n.adapt = n.adapt, n.update = n.update)

update(jm, n.iter=10000)
load.module("dic")
zm=coda.samples(jm,variable.names=c("K", "r", "sigma"), n.iter=10000, n.thin=1)
zj=jags.samples(jm,variable.names=c("K", "r", "sigma","mu"), n.iter=5, n.thin=1)


summary(zm) #Gives summary of statistics from the MCMC chain (coda)

summary(zm)$stat[2,1:2]  #Gives you the mean and SD of just r

summary(zm)$stat #Gives just the quantiles

summary(zm)$quantile

##Make a table that contains the mean, SD, median and upper and lower 2.5% CI for parameter estimates


t1 = summary(zm)$stat[1:3,1:2]
t2 = summary(zm)$quantile[1:3,c(1,3,5)]

combined.table <- cbind(t1, t2)

round(combined.table, digits = 3)


b=summary(zj$K,mean)$stat
b=summary(zj$mu,quantile,c(.025,.5,.975))$stat
par(mfrow=c(1,2))
plot(pop.data$"Population Size", pop.data$"Growth Rate", xlab="N", ylab="Per capita growth rate")
lines(pop.data$"Population Size",b[2,])
lines(pop.data$"Population Size",b[1,],lty="dashed")
lines(pop.data$"Population Size",b[3,],lty="dashed")
plot(density(zj$K),xlab="K", main="", xlim=c(800,2500))

##Exercise: understanding coda objects
jm=jags.model("Logistic example BUGS.R",data=data, inits,n.chains=length(inits), n.adapt = 500, n.update = 500)
update(jm, n.iter=20)
load.module("dic")
zm.short=coda.samples(jm,variable.names=c("K", "r", "sigma"), n.iter=20, n.thin=3)

summary(zm.short)


df = as.data.frame(rbind(zm[[1]], zm[[2]], am[[3]]))

max(df$sigma)
mean(df$K[1:1000])
nr=length(df$K)
mean(df$K[(nr-1000):nr])
plot(density(df$K),main="",xlim=c(800,2000),xlab="K")
1-ecdf(df$K)(1600)
ecdf(df$K)(1200)-ecdf(df$K)(800)

summary(zm)
traceplot(zm)
densityplot(zm)




