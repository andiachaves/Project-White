setwd("/Users/Tom/Documents/Book with Mevin/Serengeti data and modeling/")
rm(list=ls())
library(rjags)
library(boot)

y=read.csv("Wildebeest data from Ray.csv")
#y=read.csv("Serengeti wildebeest population data only.csv")
y=cbind(seq(1,nrow(y)),y)
y=y[!is.na(y$Rain),]
names(y)=c("row", "Year","Rain","N","sd")
T = nrow(y)  #last year has missing data for rain
plot(y$Year,y$N)

N.index=y[!is.na(y$N),1]
N.obs = y[!is.na(y$N),4]
sd.N = y[!is.na(y$N),5]

#Rain.norm= (y$Rain-mean(y$Rain))/mean(y$Rain)
Rain.norm=as.vector(scale(y$Rain, center=TRUE))
N = numeric(T);mu=numeric(T)
N[1] = y$N[1]
sd.proc=.01
n.coef=3
b=numeric(n.coef)
b[1]= .21
b[2] = -.00017
b[3] = .112
b[4] = .0001
for(t in 2:T){
	mu[t] = log(N[t-1]*exp(b[1] + b[2]*N[t-1] + b[3]*Rain.norm[t]))
	N[t] = rlnorm(1,mu[t], sd.proc)
}
#Plot simulated trajectory vs real data
plot(y$Year,N, typ='l', ylim=c(0,2000)) #simulated
points(y$Year,y$N)  #real

N.proc=N[N.index]
sd.obs=numeric(length(N.index))
for(j in 1:length(N.proc)){
	sd.obs[j]=runif(1,.10,.15)*N.proc[j]
	N.obs[j] = rnorm(1,N.proc[j],sd.obs[j])
	}


plot(y$Year[N.index],N.proc,pch=18, typ='l')
points(y$Year[N.index],N.obs)
points(y$Year[N.index], N[N.index], col="red")
	
n.coef=4

######section for forecast
#add data for forecasts, set npred to 0 for no forecast
n.pred = 3
#additional inits
# average
Rain=c(Rain.norm,c(rep(0,n.pred)))
N.init = c(N,rep(N[T], n.pred))

data=list(#N.obs=N.obs,   # use this for simulations
	N.obs=y$N[N.index],   #use this for data
	#Rain=y$Rain, 
	Rain=Rain,
	index=N.index,
	T= T + n.pred,
	n.obs= length(N.index),
	n.coef = n.coef,
	tau.obs = 1/y$sd[N.index]^2  #use this for data
	#tau.obs=1/sd.obs^2 #use this for simulations
	)
	
# inits=list(
	# list(b=b[1:n.coef],	logN=log(N.init)),
	# list(b=b[1:n.coef]*.8, logN=log(N.init*.7))
# )

inits=list(
	list(b=b[1:n.coef],	z=(N.init)),
	list(b=b[1:n.coef]*.8, z=(N.init*.7))
)

#Examine autocorrleation
acf.dat=NULL
for(i in 1:length(acf.dat)){
z1=acf(zj$epsilon[1:19,1,1])	
	
	
}

n.iter=50000

#####JAGS verion
jm=jags.model("Ricker JAGS lognormal 9.R", data=data, inits, n.adapt=3000, n.chains=length(inits))
update(jm,n.iter=30000)
zm=coda.samples(jm,variable.names=c("b", "sigma.p"), n.iter=50000)
zj=jags.samples(jm,variable.names=c("epsilon.proc", "epsilon.obs","z", "fit", "fit.new", "pvalue", "b"), n.iter=n.iter)
#zm.DIC=dic.samples(jm, n.iter=50000, thin=1, type="pD")
end=Sys.time()

plot(zj$fit/100,zj$fit.new/100,cex=.05,xlim=c(0,30000),ylim=c(0,30000),ylab="Simulated discrepancy", xlab="Observed discrepancy" )
abline(a=0,b=1)
p=summary(zj$pvalue,mean)$stat
text(22000,12000, paste("Bayesian P value = ", as.character(signif(p,2))))


acf_test = function(n, jags_object,ts_n){
acf.dat=NULL
for(i in 1:n){
z1=acf(jags_object[1:ts_n,i,1],plot=FALSE)$acf
z2=acf(jags_object[1:ts_n,i,2],plot=FALSE)$acf
acf.dat=rbind(acf.dat,z1,z2)
}

q=matrix(nrow=11, ncol=3)
for(j in 1:11)q[j,] = quantile(acf.dat[,j],c(.025,.5,.975))
#plot(lag,q[,2],lwd=2, col="black", main="Autocorrleation", ylim=c(-1,1), typ="h", ylab="ACF", xlab="Lag")
acf(q[,2], main="Autocorrelation", lwd=3)
lag=seq(0,10)
lines(lag,q[,1], typ="h",lty="dotted", col="red")
lines(lag,q[,3], typ="h",lty="dotted", col="red")
#abline(h=-.5, col="blue")
#abline(h=.5, col="blue")
#abline(h=0,col="blue")
} #end of acf function

acf_test(n=10000, jags_object=zj$epsilon.obs, ts_n=19)


pdf(file="Wildebeest trace plots.pdf")
par(mfrow=c(3,2))
traceplot(zm)
dev.off()



pdf(file="Wildebeest density plots.pdf")

par(mfrow=c(3,2))
plot_prior=function(x){
y.p=dnorm(x,0,1/sqrt(.00001))
lines(x,y.p, lty="dashed")}
plot(density(zj$b[1,,]),xlab=expression(italic(beta[1])),main="")
plot_prior(x=seq(0,.8,.01))

plot(density(zj$b[2,,]),xlab=expression(italic(beta[2])),main="")
abline(h=1,lty="dashed")
plot(density(zj$b[3,,]),xlab=expression(italic(beta[3])),main="")
plot_prior(x=seq(-.5,1,.01))
plot(density(zj$b[4,,]),xlab=expression(italic(beta[4])),main="")
abline(h=1,lty="dashed")

plot(density(zj$sigma.p),xlab=expression(italic(sigma[p])),main="",xlim=c(0,.3))
x=seq(0,.3,.01)
y.p=dgamma(x,.01,.01)
lines(x,y.p,lty="dashed")
dev.off()



plot.BCI=function(x,data.y, data.bar, model.y,y.high,y.low,xlab=years,ylab){
	y.errorbars(x,data.y,data.bar,y.low=y.low,y.high=y.high, 	xlab=xlab, ylab=ylab)

	lines (x,model.y[2,], typ='l', ylim=c(y.low,y.high)) #plots the median
	xx <-c(x, rev(x))
	yy <- c(c(model.y[3,]),rev(c(model.y[1,])))
	polygon(xx,yy,col=topo.colors(5,alpha=.3), border=NA)
}

y.errorbars = function(x,y,ybar, ylab=x, xlab=y, main=" ",y.low,y.high){
	plot(x,y,pch=16,,xlab=xlab, ylab = ylab, main=main, ylim=c(y.low,y.high), family="serif")
	arrows(x,y-ybar,x,y+ybar, code=3, angle=90, length=.01)
}
pdf(file="Wildebeest state plot.pdf", pointsize=14)
par(mfrow=c(1,1))
year=seq(1961,2011)
j=summary(zj$z, quantile, c(.025,.5,.975))
plot.BCI(x=year,data=c(y$N,rep(NA,3)),data.bar=2*y$sd,model.y=j$stat, xlab="Year", ylab="Population size (x1000)", y.low=0, y.high=2000)
par(mfrow=c(1,1))
dev.off()



###old acf code
#Examine autocorrleation
# acf.dat=NULL
# for(i in 1:10000){
# z1=acf(zj$epsilon.proc[1:18,i,1],plot=FALSE)$acf
# z2=acf(zj$epsilon.proc[1:18,i,2],plot=FALSE)$acf
# acf.dat=rbind(acf.dat,z1,z2)
# }

# q=matrix(nrow=11, ncol=3)
# for(j in 1:11) q[j,] = quantile(acf.dat[,j],c(.025,.5,.975))
# lag=seq(0,10)
# #plot(lag,q[,2],lwd=2, col="black", main="Autocorrleation", ylim=c(-1,1), typ="h", ylab="ACF", xlab="Lag")
# acf(q[,2], main="Autocorrelation", lwd=3)
# lines(lag,q[,1], typ="h",lty="dotted", col="red")
# lines(lag,q[,3], typ="h",lty="dotted", col="red")
# #abline(h=-.5, col="blue")
# #abline(h=.5, col="blue")
# #abline(h=0,col="blue")





