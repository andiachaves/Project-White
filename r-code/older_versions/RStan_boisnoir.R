# Authors: Bernd Panassiti, Florian Hartig, Michael Breuer, Robert Biedermann
# Date created: 10.02.2014; last modified 27.01.2015, prepared as supplement: 28.04.2015

# Supplement to "Bayesian inference of environmental and biotic factors determining the occurrence of the grapevine disease ‘bois noir’"

# Description
# RStan code to predict prevalences of a phytoplasma disease in grapevines (bois noir) in the Baden region (SW-Germany) using GLMM evaluated in a Bayesian framework.

#Removes al list objects/cleans the working space
rm(list=ls(all=TRUE))

# load libraries
library(ggplot2)
library(rstan)

# load and define data
setwd("/Users/andia/Documents/ANDIA\ documents/A-\ Postdoc\ UVI\ Brant\ Lab/Projects/EPScor/References/Bayes_diseases_/Pannassiti_etal_Supplemental\ _Info")
read.table("boisnoir.txt")
boisnoir.txt
newdata <- read.table("boisnoir.txt")
grapeselection = c(5:13)   # to select grapevines
enviroselection = c(14:26) # to select remaining environmental variables


### Here begins the RStan model
# Data as a list
dataList = list(
  observedPlantsTotal=newdata$Nsum,
  observedDisease=newdata$illsum,
  nObs=length(newdata$ill),
  
  Grapes = newdata[3:(3-1+length(grapeselection))],
  nGrapes = length(grapeselection),
  
  Environment = newdata[(3+length(grapeselection)):(3-1+length(grapeselection)+length(enviroselection))],
  nEnvironment = length(enviroselection),
  
  nRegions = length(unique(newdata$regionData)),
  nOwners = length(unique(newdata$ownerData)),
  ownerData = newdata$ownerData,
  regionData = newdata$regionData
)


# Define model

model.bn.frequentgrapes <- '
data {
// counters
int<lower=0> nObs;
int<lower=0> nGrapes;
int<lower=0> nEnvironment;  
int<lower=0> nRegions; 
int<lower=0> nOwners;

// predictors
matrix[nObs, nGrapes] Grapes;
matrix[nObs, nEnvironment] Environment;

int<lower=1, upper =253> ownerData[nObs];
int<lower=1, upper=6> regionData[nObs];

// response
int<lower=0> observedPlantsTotal[nObs]; // DISEASE - num cases, total number of plants
int<lower=0> observedDisease[nObs];     // num successes
}


parameters {
// parameters that will be used in the model should be defined here
// these will be saved and printed with the output

real b0;                           // intercept
vector[nGrapes] parGrapes;
vector[nEnvironment] parEnviro;

real <lower=0> tauRegions;
real <lower=0> tauOwners;
vector[nRegions] u_raw;
vector[nOwners] v_raw;
}

transformed parameters {
// generation of derived parameters from parameters above should be defined here
// these will be saved and printed with the output

vector[nRegions] u;
vector[nOwners] v;
real realDisease[nObs];
vector[nObs] mu;

for (j in 1:nRegions) {
u[j] <- u_raw[j]*sqrt(tauRegions);
}

for (h in 1:nOwners){
v[h] <- v_raw[h]*sqrt(tauOwners);
}


for( i in 1:nObs ) {
mu[i] <- u[regionData[i]] + v[ownerData[i]];
realDisease[i]  <-(b0  + Environment[i] * parEnviro  + Grapes[i] * parGrapes + mu[i]);
}
}

model{
// This is area for describing the model including the priors and likelihood
// This is the only place where sampling can be done within Stan
// All data and parameters used in this section must be defined above

// Fixed effects priors
b0 ~ normal(-5.5, 100 ); //
parGrapes ~ normal(0 , 10 );
parEnviro ~ normal(0 , 10 );

//Hyperpriors
tauRegions ~ inv_gamma(0.001, 0.001);
tauOwners  ~ inv_gamma(0.001, 0.001);    

//Random effects priors
u_raw ~ normal(0, 1);
v_raw ~ normal(0, 1);

observedDisease ~ binomial_logit(observedPlantsTotal,realDisease);
}

generated quantities{
// predictive inference
real realDiseasePredictions[nObs]; // vector to store predictions

for( i in 1:nObs ) {
realDiseasePredictions[i]  <-binomial_rng(observedPlantsTotal[i],inv_logit(b0  + Environment[i] * parEnviro  + Grapes[i] * parGrapes + mu[i]));
}
}
'




##### Run model - may take some hours

time <- proc.time()
set.seed(2013)
stan.fit <- stan(model_code = model.bn.frequentgrapes, data = dataList, iter = 10, chains = 4,refresh=1,thin=10)
(proc.time() - time)/60/60

plot(stan.fit)

# to access the samples and plotting them with ggmcmc we can use the extract method, discarding the warm-up period
library(ggmcmc) 
library(coda)
s <- extract(stan.fit, inc_warmup=FALSE) 

#then stan.fit can be converted into a typical coda object by
s <- mcmc.list(mcmc(s[,1,-4]), mcmc(s[,2,-4]), mcmc(s[,3,-4]), mcmc(s[,4,-4]))
# Or a cleaner way thanks to Sebastian Weber:
s <- do.call(mcmc.list, alply(s[,,-4], 2, mcmc))
# Or even a more general way, thanks to Forrest Stevens
s <- do.call(mcmc.list, alply(s[,,-(length(s[1,1,]))], 2, mcmc))

#specific for caterpillars plots per parameter I have to modify for this specific example
#Plots of ggmcmc
library(ggmcmc)
data(radon)
s.radon.short <- radon$s.radon.short

S <- ggs(s.radon.short)
ggmcmc(S)


L.radon.intercepts <- data.frame(
        Parameter=paste("alpha[", radon$counties$id.county, "]", sep=""),
        Label=radon$counties$County)
head(L.radon.intercepts)

S.full <- ggs(radon$s.radon, par_labels=L.radon.intercepts, family="^alpha")

ggs_caterpillar(S.full)
