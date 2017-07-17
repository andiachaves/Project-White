####
####  Read in Data
####

bird.df=read.csv("birds.csv",header=TRUE)
bird.df
N=dim(bird.df)[1]

####
####  Transform Response and Scale Covariates 
####

bird.df$log.spp=log(bird.df$spp)
bird.df$area.sc=scale(bird.df$area)
bird.df$temp.sc=scale(bird.df$temp)
bird.df$precip.sc=scale(bird.df$precip)

####
####  Look for Collinearity and Outliers  
####

pairs(bird.df[,8:11])  # NOTE: There appear to be two outliers
cor(bird.df[,9:11])

idx.outlier=(1:N)[bird.df$spp==min(bird.df$spp) | bird.df$area==max(bird.df$area)]
bird.sm.df=bird.df[-idx.outlier,]
n=dim(bird.sm.df)[1]

pairs(bird.sm.df[,8:11])  
cor(bird.sm.df[,9:11])  # NOTE: precip and area are highly correlated (r=-0.63)

####
####  Create Design Matrices
####


y=bird.sm.df$log.spp

X.1=model.matrix(log.spp~area.sc+temp.sc,data=bird.sm.df)
X.2=model.matrix(log.spp~precip.sc+temp.sc,data=bird.sm.df)

####
####  Fit Model 1
####

source("norm.reg.mcmc.R")

mcmc.1.out=norm.reg.mcmc(y,X.1,solve(t(X.1)%*%X.1)%*%t(X.1)%*%y,c(0,0,0),100,1,100,5000)

layout(matrix(1:2,1,2))
matplot(t(mcmc.1.out$betasave),type="l",lty=1)
plot(mcmc.1.out$s2save,type="l")
  
####
####  Fit Model 2 
####

mcmc.2.out=norm.reg.mcmc(y,X.2,solve(t(X.2)%*%X.2)%*%t(X.2)%*%y,c(0,0,0),100,1,100,5000)

layout(matrix(1:2,1,2))
matplot(t(mcmc.2.out$betasave),type="l",lty=1)
plot(mcmc.2.out$s2save,type="l")

####
####  Fit Model 1 with Metropolis-Hastings
####

source("norm.reg.MH.mcmc.R")

mcmc.1.MH.out=norm.reg.MH.mcmc(y,X.1,solve(t(X.1)%*%X.1)%*%t(X.1)%*%y,c(0,0,0),100,1,100,.1,5000)

layout(matrix(1:2,1,2))
matplot(t(mcmc.1.MH.out$betasave),type="l",lty=1)
plot(mcmc.1.MH.out$s2save,type="l")




