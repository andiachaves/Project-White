#' # Project-White
#' Bayesian Hierarchical Model
#'
#' ---
#' title: White plaque disease - Appendix A
#' author:
#' date: '08.23.2016; last modified: `r format(Sys.Date(), format="%d.%m.%Y")`'
#' output:
#'   pdf_document:
#'     number_sections: yes
#'     toc: yes
#' header-includes: \usepackage{graphicx}
#' email: \email{}
#' ---
#'
#'
library(knitr)
opts_knit$set(root.dir='../')                     # definining working directory; or normalizePath('../')
opts_chunk$set(fig.align='center',                # aligns all figures
                echo=TRUE,                        # shows r-code
                message=FALSE,                    # suppresses library outputs
                warnings=FALSE,                   # suppresses library outputs
                tidy=TRUE,  # prevents the source code from running off a pdf page
                dev='pdf')                        # pdf device
#'
#'
#'
#' # Introduction
#' RStan code to predict prevalences of White PLague disease in corals in the US-Virgin Islands using GLMM evaluated in a Bayesian framework.
#'
#'
#' # Data
#'
#'
#'
#' # Joint model
#' ## Loading data and librarie
#'
rm(list=ls(all=TRUE))
# Load settings
source('r-code/00_settings.r')


# load data
White <- read.csv("data/Deep.csv", header = TRUE)

# load libraries
if(!require(formatR))
{
  install.packages("formatR", repos = mir)
  require(formatR)
}

if(!require(StanHeaders))
{
  install.packages("StanHeaders", repos = mir)
  require(StanHeaders)
}

if(!require(rstan))
{
  install.packages("rstan", repos = mir)
  require(rstan)
}

if(!require(MCMCpack)) # rwish function
{
  install.packages("MCMCpack", repos = mir)
  require(MCMCpack)
}

if(!require(modeest))
{
  install.packages("modeest", repos = mir)
  require(modeest)
}

if(!require(IDPmisc))
{
  install.packages("IDPmisc", repos = mir)
  require(IDPmisc)
}

if(!require(ggplot2))
{
  install.packages("ggplot2", repos = mir)
  require(ggplot2)
}
#'
#'
#' ## Data preparation
#'
#Environmental variables include temperature and rainfall
enviroselection = c(3:4) # to select environmental variables
#'
#'
#'
#' # Data as a list
#'
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
#'
#'
#'
#'
#' ## Model definition
#'
#'
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
b0 ~ normal(0, 100 ); //
parEnvironment ~ normal(0 , 100 );

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
#'
#'
#' ## Running the model
#'
time <- proc.time()
set.seed(2016)
stan.white <- stan(model_code = model.bn.Wplague, data = dataList, iter = 30000, chains = 4,refresh=1)
(proc.time() - time)/60/60
#'
#'
#'
#'

#'
# note: folder figures synchronized in cloud (dropbox) but not results folder

save(dataList, stan.white, file = paste("results/", today, "_fit.fullModel.Rdata",sep = ""))
save(dataList, file = paste("results/", today, "_fit.fullModel_data.Rdata",sep = ""))

library(lattice)
trellis.device("pdf", file = paste("figures/", today, "_Output_rstan.pdf",
sep = ""), width = 14, height = 7)
plot(stan.white)
dev.off()
#'
#'
#' # Session info
#'
devtools::session_info()
#'
#'
