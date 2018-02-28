data {
// counters
int<lower=0> nWP;                      // number of White Plague cases, int=integer  and non negative >0
int<lower=1> nEnvi;                    // number of environmental variables

int<lower=1> nYear;               // number of years
int<lower=1,upper=4> Year[nWP];  // years

int<lower=1> nDepth;               // number of depths
int<lower=1,upper=10> Depth[nWP];  // depth

//int<lower=1> nSiteYear;               // number of sites-month combinations
//int<lower=1,upper=130> SiteYear[nWP];  // sites-month combinations

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

//real<lower=0> sigmalev_2;
//vector[nSiteYear] eta_2;

real<lower=0> sigmalev_3;
vector[nDepth] eta_3;

}

transformed parameters {
vector[nWP] RealDisease;

vector[nYear] year_raw;
year_raw  = eta_1;

//vector[nSiteYear] siteyear_raw;
//siteyear_raw  = eta_2 * sigmalev_2;

vector[nDepth] depth_raw;
depth_raw  = eta_3 * sigmalev_3;

RealDisease =  Environment * parEnvironment + year_raw[Year]+depth_raw[Depth];
} //closes transformed parameter




model{
// Fixed effects priors
b0             ~ normal(0, 100);
parEnvironment ~ normal(0, 100);


//Random effects priors
eta_1 ~ normal(0, sigmalev_1);
//eta_2 ~ normal(0, sigmalev_2);
eta_3 ~ normal(0, sigmalev_3);

//sigmalev_1 ~ uniform(0,100);
//sigmalev_1 ~ student_t(1,0,1);
sigmalev_1 ~ cauchy(0,2.5);
//sigmalev_1 ~ inv_gamma(0.001, 0.001);
//sigmalev_2 ~ uniform(0,100);
sigmalev_3 ~ uniform(0,100);

Diseasedcolonies ~ binomial_logit(TotalObservedcolonies,RealDisease);
}


generated quantities{ //predictive inference
real DiseasePredictions[nWP]; //vector to store predictions
for(i in 1:nWP) DiseasePredictions [i] = binomial_rng(TotalObservedcolonies[i],inv_logit(RealDisease[i]));

}
