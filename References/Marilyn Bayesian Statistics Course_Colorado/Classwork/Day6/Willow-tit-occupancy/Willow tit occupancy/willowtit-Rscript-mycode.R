setwd("C:\\Users\\MBrandt\\Dropbox\\Transfer documents\\conferences\\Bayesian Statistics\\Classwork\\Day6\\Willow-tit-occupancy\\Willow tit occupancy\\")


#### Load rjags
library(rjags)

####Read in cavitation data
y = read.csv("Swiss BB data.csv")

###Scale elevation data
elev <- as.vector(scale(y[,"elev"],center=TRUE))

###Create data list
bird_data = list(N = 237, elev=elev, forest=y[,5], n=y[,3])

###Set initial conditions to inital MCMC chains. Provide 3 lists of starting
#values for all ROOT nodes.
N = 237
inits = list(z=rep(1,N), b0=0.1, b1=0.1, b2=0.1, b3=0.1, p.detect=0.2)

n.adapt = 1000
n.iter = 25000
n.chains = 1

jm1 = jags.model("willowtit-JAGS-mycode.R", 
 data=bird_data,
 inits=inits,
 n.chains=n.chains,
 n.adapt=n.adapt)


load.module("dic")
zc1 = coda.samples(jm1,variable.names=c("phi","z", "p.detect","b0", "b1", "b2", "b3"),
  n.iter=n.iter)

#### With below command, R will ask you to hit enter to show each page of 
#### output/graphics:
devAskNewPage(ask = TRUE)

#### For ease of coding, create variables that indicate the start and end
#### iterations for the MCMC chains associated with the zc1 coda object:
st=start(zc1)
en=end(zc1) 
th=thin(zc1)

#### View trace plots of each monitored quantity:
xyplot(window(zc1, start=st, end=en, thin=1),layout=c(2,4))
#### View density plots of each monitored quantity:
densityplot(window(zc1, start=st, end=en, thin=1), layout=c(2,4))

#### Gelman-Rubin convergence diagnotics (want r < 1.2)
gelman.plot(zc1, layout=c(2,4))
gelman.diag(zc1)







