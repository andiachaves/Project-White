# Author: Andia Chaves Fonnegra. Code based on Panassiti et al. 2015. ECOSPHERE 6(8): 1-13.
# Date created: 08.23.2016; last modified xx.xx.xxxx, prepared as supplement: xx.xx.xxxx.

# Supplement to "Environmental Drivers of White Plague Prevalence in the US Virgin Islands"

# Description
# RStan code to predict prevalences of White PLague disease in corals in the US-Virgin Islands using GLMM evaluated in a Bayesian framework.

#Removes al list objects/cleans the working space
rm(list=ls(all=TRUE))
#-----------------------------------------------------------------------
#This code below  installs rstand 2.60 which finally run the 8schools example. The 2.10 version did not worked, and the 2.11 version did not installed.
        
#devtools::install_url("http://cran.r-project.org/src/contrib/Archive/BH/BH_1.55.0-3.tar.gz")

#devtools::install_url("https://github.com/stan-dev/rstan/releases/download/v2.6.0/rstan_2.6.0.tar.gz",dependencies = FALSE)
#library(rstan)
#set_cppo("fast")
#----------------------------------------------------------------------
#I am using R studio 3.3 and Stan 2.6        

# load libraries
library(StanHeaders)
library(ggplot2)
library(rstan)


# load and define data
setwd("/Users/andia/Documents/ANDIA\ documents/A-\ Postdoc\ UVI\ Brant\ Lab/Projects/EPScor/Rasta\ proj\ analysis/White_Plague_Bayesian-Andia")
White <- read.csv("White.csv", header = TRUE)


White

#Environmental variables include temperature and rainfall
enviroselection = c(5:6) # to select environmental variables


# Data as a list
dataList = list(
        
        
        TotalObservedcolonies = White$Ntot,
        Diseasedcolonies = White$ill,
        nWP = length(White$Plague),
        
        Temperature = White$Temp,
        Rainfall = White$Rain,
        Environment = White[5:(5-1+length(enviroselection))],
        nEnvi = length(enviroselection),
        
        Depth = White$Depth,
        Month = White$Month,
        nDepth = length(unique(White$Depth)),
        nMonth = length(unique(White$Month))
)

# Define model

model.bn.Wplague <- '
data {
// counters
int<lower=0> nWP; //number of White Plague cases, int=integer  and non negative >0 
int<lower=0> nEnvi; //number of environmental variables
int<lower=1, upper=8> nDepth; //number of depths at which was sampled
int<lower=1, upper=10> nMonth; //number of Months sampled


// predictors
matrix[nWP, nEnvi] Environment;

int<lower=1, upper =8> Depth[nWP];
int<lower=1, upper=10> Month[nWP];

// response
int<lower=0> TotalObservedcolonies [nWP]; // Total number of coral colonies (cases)
int<lower=0> Diseasedcolonies [nWP]; // number of diseased coral colonies (successes)
}

parameters {
// parameters that will be used in the model should be defined here
// these will be saved and printed with the output

real b0;                           // intercept
vector[nEnvi] parEnvironment;  //Environmental-level

real <lower=1, upper=8> tauDepth;
real <lower=1> tauMonth;
vector[nDepth] u_raw;
vector[nMonth] v_raw;
}

transformed parameters {
// generation of derived parameters from parameters above should be defined here
// these will be saved and printed with the output

vector[nDepth] u; //depth effects
vector[nMonth] v; //month effects
real RealDisease[nWP];
vector[nWP] mu; //

for (j in 1:nDepth) {
u[j] <- u_raw[j]*(tauDepth);
}

for (h in 1:nMonth){
v[h] <- v_raw[h]*(tauMonth);
}

for( i in 1:nWP ) {
mu[i] <- u[Depth[i]] + v[Month[i]];
RealDisease [i]  <-(b0  + Environment[i] * parEnvironment + mu[i]);
}
}


model{
// This is area for describing the model including the priors and likelihood
// This is the only place where sampling can be done within Stan
// All data and parameters used in this section must be defined above

// Fixed effects priors
b0 ~ normal(-5.5, 100 ); //
parEnvironment ~ normal(0 , 10 );

//Hyperpriors
tauDepth ~ inv_gamma(0.001, 0.001);
tauMonth  ~ inv_gamma(0.001, 0.001);    

//Random effects priors
u_raw ~ normal(0, 1);
v_raw ~ normal(0, 1);

Diseasedcolonies ~ binomial_logit(TotalObservedcolonies,RealDisease);
}


#predictive inference
generated quantities
{
real DiseasePredictions[nWP]; //vector to store predictions

for( i in 1:nWP ) {
DiseasePredictions [i] <- binomial_rng(TotalObservedcolonies[i],inv_logit(b0  + Environment[i] * parEnvironment + mu[i]));
}
}
'

##### Run model - may take some hours

time <- proc.time()
set.seed(2016)
stan.white <- stan(model_code = model.bn.Wplague, data = dataList, iter = 10, chains = 2,refresh=1,thin=1)
(proc.time() - time)/60/60

plot(stan.white)



