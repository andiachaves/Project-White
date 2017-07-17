rm(list=ls())
library(rjags)
setwd("/Users/Tom/Documents/NSF Statistics Workshop/2013 Course/Hobbs course materials/Dynamic models/Lynx lab/")
y=read.csv("Lynx data.csv")

#Function to get beta shape parameters from moments
shape_from_stats <- function(mu = mu.global, sigma = sigma.global){
		 a <-(mu^2-mu^3-mu*sigma^2)/sigma^2
		 b <- (mu-2*mu^2+mu^3-sigma^2+mu*sigma^2)/sigma^2
		shape_ps <- c(a,b)
		return(shape_ps)
}



#get shape parameters for population multiplier, 1/p
shapes=shape_from_stats(.163,.012)
#check prior on p
x = seq(0,1,.001)
p=dbeta(x,shapes[1],shapes[2])
plot(x,p,typ="l",xlim=c(.1,.3))


#visually estimate some data for initial conditions
endyr = nrow(y)
n=numeric(endyr+1)
mu=numeric(endyr+1) #use this for family groups
lambda=1.02
sigma.p=.00001
n[1] = y$census[1]

for(t in 2: (endyr+1)){
	n[t] <- lambda*(y$census[t-1] - .16 * y$harvest[t-1])  #use this for family groups
	}
plot(y$census,ylim=c(0,100))
lines(n)


#Harvest two years out
h=c(0, 10, 25, 50, 75)

#Data for JAGS
data = list(
	y1 = y$census[1],
	y.endyr = endyr,
	y.a=shapes[1],
	y.b=shapes[2],
	y.H=y$harvest,
	y=y$census,
	h=h
)




inits = list(
	list(
	lambda = 1.2,
	sigma.p = .1,
	N=n
	),
	list(
	lambda = 1.01,
	sigma.p = .2,
	N=round(n*1.2)
	),
	list(
	lambda = 1.06,
	sigma.p = .5,
	N=round(n)
	))


model = "Lynx Harvest JAGS.R"

n.update=10000
n.iter=50000
n.adapt=5000
n.thin=1

jm = jags.model(model,data=data,inits=inits, n.adapt=n.adapt, n.chains=length(inits))
update(jm, n.iter=n.update)


z = coda.samples(jm,variable.names=c("lambda","sigma.p","p"), n.iter=n.iter, thin=n.thin)

zj=jags.samples(jm,variable.names=c("N","N.hat", "fg", "fg.hat", "pvalue", "fit", "fit.new"), n.iter=n.iter, thin=n.thin)
summary(z)
plot(z)
heidel.diag(z)
gelman.diag(z)
#look at Bayesian P
zj$pvalue

#Do goodness of fit plot.
par(mfrow=c(1,1))
plot(zj$fit.new,zj$fit, xlab = "Discrepancy observed", ylab= "Discrepancy simulated", cex=.05, xlim=c(0,3000), ylim=c(0,3000))
abline(0,1)
p=summary(zj$pvalue,mean)$stat
text(500,2500, paste("Bayesian P value = ", as.character(signif(p,2))))

par(mfrow=c(1,1))
years=seq(1997,2010)
fg = summary(zj$fg,quantile,c(.025,.5,.975))$stat
y2=c(y$census, NA)
plot(years,y2, ylim=c(0,100), ylab="Number of Lynx Family Groups", xlab="Years")
lines(years,fg[2,])
lines(years,fg[3,], lty="dashed")
lines(years,fg[1,],lty="dashed")

lower = 26
upper = 32
p.in = numeric(length(h))
p.over =numeric(length(h))
p.under = numeric(length(h))

#calculate probability of meeting goals
for(j in 1:length(h)){
	p1 = ecdf(zj$fg.hat[j,,])(upper)
	p.under[j] = ecdf(zj$fg.hat[j,,])(lower)
	p.in[j] = p1 - p.under[j]
	p.over[j] = 1-p1
}


table = as.data.frame(cbind(h,p.under,p.in,p.over))
names(table)=c("Harvest", "P(under)", "P(in)", "P(over)")

#Probabliy population is < 100 next year with no harvest
ecdf(d[,14])(100)



	
	

