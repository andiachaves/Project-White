##Lynx model call

setwd("C:\\Users\\MBrandt\\Dropbox\\Transfer documents\\conferences\\Bayesian Statistics\\Classwork\\Day8\\Lynx-lab-materials\\")
library(rjags)

lynx.data = read.csv("Lynx data.csv")

#Function to get beta shape parameters from moments
shape_from_stats <- function(mu = mu.global, sigma = sigma.global){
		 a <-(mu^2-mu^3-mu*sigma^2)/sigma^2
		 b <- (mu-2*mu^2+mu^3-sigma^2+mu*sigma^2)/sigma^2
		shape_ps <- c(a,b)
		return(shape_ps)
}
#get shape parameters for population multiplier, 1/p
shapes=shape_from_stats(.163,.012)
#check prior on phi
x = seq(0,1,.001)
phi=dbeta(x,shapes[1],shapes[2])
plot(x,phi,typ="l",xlim=c(.1,.3))

##Initial conditions
endyr = nrow(lynx.data) #counts the number of data rows = 13
n = numeric(endyr+1) #creates a vector of 14 zeros
mu = numeric(endyr+1)
lambda = 1.02
sigma.proc = 0.00001
n[1] = lynx.data$census[1] #sets the first value in n = to first value in census data

for(t in 2: (endyr+1)) {
	n[t] <- lambda*(lynx.data$census[t-1] - 0.163 * lynx.data$harvest[t-1])
}
plot(lynx.data$census, ylim = c(0,100)) # plots the census data as points
lines(n) # plots the n values calculated above as a line

#Proposed harvest levels
harv = c(0, 10, 25, 50, 75)


#Data to send to JAGSs

data = list(
	y1 = lynx.data$census[1], #initial value for y
	y.endyr = endyr, 
	y.a = shapes[1],
	y.b=shapes[2],
	y.H=lynx.data$harvest,
	y=lynx.data$census,
	harv=harv
)

inits = list(
	list(
	lambda = 1.2,
	sigma.proc = .1,
	N=n
	),
	list(
	lambda = 1.01,
	sigma.proc = .2,
	N=round(n*1.2)
	),
	list(
	lambda = 1.06,
	sigma.proc = .5,
	N=round(n)
	))


n.update = 10000
n.iter = 50000
n.adapt = 5000
n.thin = 1

jm = jags.model("LynxLab_JAGS-mycode.R", data=data, inits=inits, n.adapt = n.adapt, n.chains = length(inits))
update(jm, n.iter=n.update)

load.module("dic")  #Must load this if you want to monitor for deviance
z = coda.samples(jm,variable.names=c("deviance","lambda","sigma.proc","phi","tau.proc"), n.iter=n.iter, thin=n.thin)

### Plots of traces
st=start(z)
en=end(z) 
th=thin(z)

xyplot(window(z, start=st, end=en, thin=1)[,c("lambda", "sigma.proc", "phi", "tau.proc", "deviance")], layout=c(1,5))
xyplot(window(z, start=st, end=en, thin=1)[,c("deviance")])

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
y2=c(lynx.data$census, NA)
plot(years,y2, ylim=c(0,100), ylab="Number of Lynx Family Groups", xlab="Years")
lines(years,fg[2,])
lines(years,fg[3,], lty="dashed")
lines(years,fg[1,],lty="dashed")


# Environmentalists and hunters have agreed on an acceptable range for lynx abundance
# in the unit, 26-32 family groups
lower = 26 # lower number of family groups
upper = 32 # upper number of family groups
p.in = numeric(length(harv))  #probabilities that the harvest is within the range
p.over =numeric(length(harv))  #probabilities that the harvest is over the range
p.under = numeric(length(harv)) #probabilities that the harvest is less than the range

#calculate probability of meeting goals
for(j in 1:length(harv)){					
	p1 = ecdf(zj$fg.hat[j,,])(upper)
	p.under[j] = ecdf(zj$fg.hat[j,,])(lower)
	p.in[j] = p1 - p.under[j]
	p.over[j] = 1-p1
}


table = as.data.frame(cbind(harv,p.under,p.in,p.over))
names(table)=c("Harvest", "P(under)", "P(in)", "P(over)")

#Probabliy population is < 100 next year with no harvest
p100 = ecdf(zj$fg.hat[1,,])(100)



	
	


