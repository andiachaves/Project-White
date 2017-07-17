#File Logistic example.R in the folder below

#setwd("/Users/Tom/Documents/Ecological Modeling Course/JAGS Primer")
#r=.2
#K=1200
#sigma=.03
#x=c(70,283,332,459,502,678,690,804,923,1005,1120)
#mu=r-r/K*x
#plot(x,mu)
#
#y=rnorm(x,mu,sigma)
#plot(x,y)
#data = as.matrix(cbind(x,y))
#names(data)=c("Population Size", "Growth Rate")
#write.csv(data, "Logistic Data")

setwd("/Users/Tom/Documents/Ecological Modeling Course/JAGS Primer")
rm(list=ls())
pop.data=(read.csv("Logistic Data II.csv"))
names(pop.data)=c("Year","Population Size", "Growth Rate")

inits=list(
list(K=1500, r=.2, tau=2500),# , #chain 1
list(K=1000, r=.15, tau=1000))#, #chain 2
#list(K=900, r=.3, tau=500) #chain 3
#)  

n.xy = nrow(pop.data)

data=list(n=n.xy, x=as.real(pop.data$"Population Size"),y=as.real(pop.data$"Growth Rate"))
library(rjags)

##call to JAGS
jm=jags.model("Logistic example BUGS.R",data=data, inits,n.chains=length(inits), n.adapt = 5000)

update(jm, n.iter=10000)
load.module("dic")
zm=coda.samples(jm,variable.names=c("K", "r", "sigma"), n.iter=10000, n.thin=1)
zj=jags.samples(jm,variable.names=c("K", "r", "sigma","mu"), n.iter=5, n.thin=1)

dic.j= dic.samples(jm,n.iter=2500, type="pD")


acfplot(zm)


b=summary(zj$K,mean)$stat
b=summary(zj$mu,quantile,c(.025,.5,.975))$stat
par(mfrow=c(1,2))
plot(pop.data$"Population Size", pop.data$"Growth Rate", xlab="N", ylab="Per capita growth rate")
lines(pop.data$"Population Size",b[2,])
lines(pop.data$"Population Size",b[1,],lty="dashed")
lines(pop.data$"Population Size",b[3,],lty="dashed")
plot(density(zj$K),xlab="K", main="", xlim=c(800,2500))


#save(zj, file="/Users/Tom/Documents/Ecological Modeling Course/JAGS Primer/zj_object.Rdata")


#exercises on manipulating objects
df=as.data.frame(rbind(zm[[1]],zm[[2]],zm[[3]]))
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








