# Project-White
Bayesian Hierarchical Model

---
title: White plaque disease - Appendix A
author: 
date: '08.23.2016; last modified: `r format(Sys.Date(), format="%d.%m.%Y")`'
output: 
  pdf_document: 
    number_sections: yes
    toc: yes
header-includes: \usepackage{graphicx}
email: \email{}
---

```{r setup, include=FALSE, warnings=FALSE}
library(knitr)
opts_knit$set(root.dir='../')                     # definining working directory; or normalizePath('../')
opts_chunk$set(fig.align='center',                # aligns all figures
                echo=TRUE,                        # shows r-code
                message=FALSE,                    # suppresses library outputs
                warnings=FALSE,                   # suppresses library outputs
                tidy=TRUE,  # prevents the source code from running off a pdf page
                dev='pdf')                        # pdf device
```


# Introduction
RStan code to predict prevalences of White PLague disease in corals in the US-Virgin Islands using GLMM evaluated in a Bayesian framework.  


# Data



# Joint model
## Loading data and librarie
```{r}
rm(list=ls(all=TRUE))
# Date format - forms part of names of created files or graphs
today<-format(Sys.time(), "%Y.%m.%d")


# load data & libraries
White <- read.csv("data/Deep.csv", header = TRUE)

# load libraries
library(formatR)
library(StanHeaders)
library(rstan)
library(MCMCpack) # rwish function
library(modeest)
library(IDPmisc)
library(ggplot2)
# create directory if not existing
suppressWarnings(dir.create(file.path(getwd(),"results/04_Rstan")))
```

## Data preparation
```{r data_preparation}
#Environmental variables include temperature and rainfall
enviroselection = c(3:4) # to select environmental variables
```


# Data as a list
```{r}
Month=c(sort(rep(1:18,4)))
dataList = list(
        TotalObservedcolonies = White$Ntot,
        Diseasedcolonies = White$ill,
        nWP = length(White$Plague),
        
        Temperature = White$Temp,
        Rainfall = White$Rain,
        Environment = White[3:(3-1+length(enviroselection))],
        nEnvi = length(enviroselection),

        #Month = White$Month,
        Month = Month,
        nMonth = length(unique(Month))
)
```



## Model definition

```{r define_model}
# Define model

model.bn.Wplague <- '
data {
// counters
int<lower=0> nWP;    //number of White Plague cases, int=integer  and non negative >0 
int<lower=0> nEnvi;  //number of environmental variables
int<lower=0> nMonth; //number of Months sampled


// predictors
matrix[nWP, nEnvi] Environment;
int<lower=0> Month[nWP];

// response
int<lower=0> TotalObservedcolonies [nWP]; // Total number of coral colonies (cases)
int<lower=0> Diseasedcolonies [nWP]; // number of diseased coral colonies (successes)
}

parameters {
// parameters that will be used in the model should be defined here
// these will be saved and printed with the output

real b0;                        // intercept
vector[nEnvi] parEnvironment;   //Environmental-level


real <lower=0> tauMonth;
vector[nMonth] v_raw;
}

transformed parameters {
// generation of derived parameters from parameters above should be defined here
// these will be saved and printed with the output


vector[nMonth] v; //month effects
real RealDisease[nWP];
vector[nWP] mu; //


for (h in 1:nMonth){
v[h] = v_raw[h]*(tauMonth);
}

for( i in 1:nWP ) {
mu[i] = v[Month[i]];
RealDisease[i]  = (b0  + Environment[i] * parEnvironment + mu[i]);
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
tauMonth  ~ inv_gamma(0.001, 0.001);    

//Random effects priors
v_raw ~ normal(0, 1);

Diseasedcolonies ~ binomial_logit(TotalObservedcolonies,RealDisease);
}


#predictive inference
#generated quantities
#{
#real DiseasePredictions[nWP]; //vector to store predictions

#for( i in 1:nWP ) {
#DiseasePredictions [i] <- binomial_rng(TotalObservedcolonies[i],inv_logit(b0  + Environment[i] * parEnvironment + mu[i]));
#}
#}
'
```

## Running the model
```{r test_run,echo=FALSE,eval=FALSE}
stan.white <- stan(model_code = model.bn.Wplague, data = dataList, iter = 1000, chains = 1,refresh=1,thin=1)
```


```{r model_run,eval=FALSE}
time <- proc.time()
set.seed(2016)
stan.white <- stan(model_code = model.bn.Wplague, data = dataList, iter = 30000, chains = 4,refresh=1)
(proc.time() - time)/60/60
```


# Session info
```{r}
devtools::session_info()
```


```{r ploting from stan to ggmcmc}
print(stan.white)
plot(stan.white)

library(ggmcmc)
library(coda)

#This is to Install Stan Caterpillar StanCat to build caterpillar graphs from Stan
#devtools::install_github('christophergandrud/StanCat')
StanCat::stan_caterpillar(stan.white, pars = 'v_raw\\[.*\\]')

#To use ggmcmc#the rstan::extract, specify that is extracted from rstan procedure and not from dpylr
library(plyr)
#this one extracts everything but is not plotting histograms properly
s <- rstan::extract(stan.white, inc_warmup=FALSE, permuted=FALSE) 
#this one extracts only some functions

s <- mcmc.list(mcmc(s[,1,-4]), mcmc(s[,2,-4]), mcmc(s[,3,-4]), mcmc(s[,4,-4]))
S <- ggs(s) 
ggmcmc(S)
ggs_traceplots(S)
ggs_caterpillar(S)


#for specific parameters. But is not runing
s <- rstan::extract(stan.white, pars = 'v_raw', permuted = TRUE)$parEnvironment 
s <- do.call(mcmc.list, alply(s[,,-4], 2, mcmc))

```

```{r-To save the r markdown file}

save(dataList, stan.white, file = paste("Figures/", today, "_fit.fullModel.pdf",
sep = ""))
save(dataList, file = paste("Figures/", today, "_fit.fullModel_data.Rdata",
sep = ""))
library(lattice)
trellis.device("pdf", file = paste("Figures/", today, "_Output_rstan.pdf",
sep = ""), width = 14, height = 7)
plot(stan.white)
dev.off()


library(lattice)
trellis.device("pdf")
file= paste("Figures/",today,"_Output_rstan.pdf",sep="2-2016", width=14, height=7)
plot(stan.white)
```

