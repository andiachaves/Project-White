data {
// counters
int<lower=0> nWP;                      // number of White Plague cases, int=integer  and non negative >0
int<lower=1> nEnvi;                    // number of environmental variables

int<lower=1> nSiteMonth;               // number of sites-month combinations
int<lower=1,upper=130> SiteMonth[nWP];  // sites-month combinations

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
real<lower=0> sigmalev_1;
vector[nSiteMonth] eta_1;
}



transformed parameters {
vector[nWP] RealDisease;
vector[nSiteMonth] sitemonth_raw;

sitemonth_raw  = eta_1 * sigmalev_1;

RealDisease =  Environment * parEnvironment + sitemonth_raw[SiteMonth];
} //closes transformed parameter




model{
// Fixed effects priors
b0             ~ normal(0, 10);
parEnvironment ~ normal(0, 10);

//Random effects priors
eta_1 ~ normal(0, 100);
sigmalev_1 ~ inv_gamma(0.001, 0.001);

Diseasedcolonies ~ binomial_logit(TotalObservedcolonies,RealDisease);
}



generated quantities{ //predictive inference
real DiseasePredictions[nWP]; //vector to store predictions
for(i in 1:nWP) DiseasePredictions [i] = binomial_rng(TotalObservedcolonies[i],inv_logit(RealDisease[i]));

}

