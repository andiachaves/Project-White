data {
// counters
int<lower=0> nWP;                      // number of White Plague cases, int=integer  and non negative >0
int<lower=1> nEnvi;                    // number of environmental variables

int<lower=1> nYears;               // number of sites-month combinations
int<lower=1,upper=3> Year[nWP];  // sites-month combinations

int<lower=1> nSites;               // number of sites-month combinations
int<lower=1,upper=13> Site[nWP];  // sites-month combinations

matrix[nWP, nEnvi] Environment;

// response
int<lower=0> TotalObservedcolonies [nWP]; // Total number of coral colonies (cases)
int<lower=0> Diseasedcolonies [nWP]; // number of diseased coral colonies (successes)
}
parameters {
real b0; // intercept
// Fixed effects
vector[nEnvi] parEnvironment;   // Environmental-level

// Random effects
//vector[nSites] site_raw;
//real<lower=0> tauSites;
vector[nYears] year_raw;
real<lower=0> tauYears;
}



transformed parameters {
vector[nWP] RealDisease;

RealDisease =  Environment * parEnvironment + year_raw[Year];
//+ site_raw[Site];
} //closes transformed parameter


model{
// Fixed effects priors
b0             ~ normal(0, 10);
parEnvironment ~ normal(0, 10);

//Random effects priors
//site_raw ~ normal(0, tauSites);
year_raw ~ normal(0, tauYears);

//Hyperpriors
//tauSites  ~ inv_gamma(0.001, 0.001);
tauYears  ~ inv_gamma(0.001, 0.001);

// other options
//uniform(0,100);
//student_t(1,0,1);
//tauYears  ~ cauchy(0,2.5);

Diseasedcolonies ~ binomial_logit(TotalObservedcolonies,RealDisease);
}


generated quantities{ //predictive inference
real DiseasePredictions[nWP]; //vector to store predictions
for(i in 1:nWP) DiseasePredictions [i] = binomial_rng(TotalObservedcolonies[i],inv_logit(RealDisease[i]));

}
