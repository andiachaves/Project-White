setwd("C:\\Users\\MBrandt\\Dropbox\\Transfer documents\\conferences\\Bayesian Statistics\\Classwork\\Day4\\Morning lab\\Day_4_AM_Lab_1\\Day_4_AM_Lab_1\\")

####
####  Read in Data
####

bird.df=read.csv("birds.csv",header=TRUE)
bird.df
N=dim(bird.df)[1] 

####
####  Transform Response and Scale Covariates 
####

#Transform the response variable using a log 
#We transform the data because the model is normal (not Poisson)
bird.df$log.spp=log(bird.df$spp) #Log transforms species richness data

#Scale the covariates (area, temperature and precipitation)
bird.df$area.sc=scale(bird.df$area) #scales state area data
bird.df$temp.sc=scale(bird.df$temp) #scales temperature data
bird.df$precip.sc=scale(bird.df$precip) #scales precipitation data

####
####  Look for Collinearity and Outliers  
####

#Create a scatterplot matrix (8:11 indicate the columns - which have the
# transformed and scaled data:
pairs(bird.df[,8:11])  # NOTE: There appear to be two outliers (possibly Alaska and Hawaii)
# Calculate correlations between covariates:
cor(bird.df[,9:11])

#This function identifies the outliers
idx.outlier=(1:N)[bird.df$spp==min(bird.df$spp) | bird.df$area==max(bird.df$area)]
#This removes the outliers (Alaska and Hawaii)
bird.sm.df=bird.df[-idx.outlier,]
n=dim(bird.sm.df)[1] #n now indicates the number of states without Alaska and Hawaii


pairs(bird.sm.df[,8:11])  
cor(bird.sm.df[,9:11])  # NOTE: precip and area are highly correlated (r=-0.63)

####
####  Create Design Matrices
####


y=bird.sm.df$log.spp  #y is the new set of log transformed spp data without Alaska and Hawaii

X.1=model.matrix(log.spp~area.sc+temp.sc,data=bird.sm.df)
X.2=model.matrix(log.spp~precip.sc+temp.sc,data=bird.sm.df)

####
####  Fit Model 1
####

#Source the function to do the regression (to fit the model)
source("norm.reg.mcmc.R")

#Arguments for the norm.reg.mcmc function: (y,X,betastrt,betamean,betavar,s2mean,s2sd,ngibbs,no.print=FALSE)
# y - data (log transformed species diversity data)
# X - matrix with predictor variables
# betastrt: solve(t(X.1)%*%X.1)%*%t(X.1)%*%y
# betamean: c(0,0,0)
# betavar: 100
# s2mean: 1
# s2sd: 100 
# ngibbs: 5000

#Fit model where the covariates are the scaled versions of "area" and temp" with a normal prior
# for Beta and inverse gamma prior for sigma^2
mcmc.1.out=norm.reg.mcmc(y,X.1,solve(t(X.1)%*%X.1)%*%t(X.1)%*%y,c(0,0,0),100,1,100,5000)

layout(matrix(1:2,1,2))
matplot(t(mcmc.1.out$betasave),type="l",lty=1)
plot(mcmc.1.out$s2save,type="l")

###CHALLENGE:  Change s2sd to 0.1 from 100
mcmc.1.out=norm.reg.mcmc(y,X.1,solve(t(X.1)%*%X.1)%*%t(X.1)%*%y,c(0,0,0),100,1,0.1,5000)
layout(matrix(1:2,1,2))
matplot(t(mcmc.1.out$betasave),type="l",lty=1)
plot(mcmc.1.out$s2save,type="l")
#Change to 0.01 from 0.1
mcmc.1.out=norm.reg.mcmc(y,X.1,solve(t(X.1)%*%X.1)%*%t(X.1)%*%y,c(0,0,0),100,1,0.01,5000)
layout(matrix(1:2,1,2))
matplot(t(mcmc.1.out$betasave),type="l",lty=1)
plot(mcmc.1.out$s2save,type="l")


  
####
####  Fit Model 2 
####

#Fit model where the covariates are the scaled versions of precip and temp:
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

hist

###*** See notes for outputs:







