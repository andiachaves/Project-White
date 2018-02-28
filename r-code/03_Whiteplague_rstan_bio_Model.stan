data {
// counters
int<lower=0> nWP;                      // number of White Plague cases, int=integer  and non negative >0
int<lower=1> nEnvi;                    // number of environmental variables

int<lower=1> nYear;               // number of sites-month combinations
int<lower=1,upper=130> Year[nWP];  // sites-month combinations

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
vector[nYear] eta_1;
}



transformed parameters {
vector[nWP] RealDisease;
vector[nYear] year_raw;

year_raw  =eta_1 * sigmalev_1;

RealDisease =  Environment * parEnvironment + year_raw[Year];
} //closes transformed parameter




model{
// Fixed effects priors
b0             ~ normal(0, 10);
parEnvironment ~ normal(0, 10);

//Random effects priors
eta_1 ~ normal(0, 100);

sigmalev_1 ~ uniform(0,100);
//sigmalev_1 ~ student_t(1,0,1);
//sigmalev_1 ~ cauchy(0,2.5);
//sigmalev_1 ~ inv_gamma(0.001, 0.001);

Diseasedcolonies ~ binomial_logit(TotalObservedcolonies,RealDisease);
}


generated quantities{ //predictive inference
real DiseasePredictions[nWP]; //vector to store predictions
for(i in 1:nWP) DiseasePredictions [i] = binomial_rng(TotalObservedcolonies[i],inv_logit(RealDisease[i]));

}
