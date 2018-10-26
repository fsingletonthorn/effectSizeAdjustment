library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Run BMWMod_additional.R before running this for the data (which in turn has dependencies)



cat("


", "model.stan")

stanData <- list(orgEffect_FT = jagData$fis.o, 
     n = length(jagData$fis.o),
     repTau = 1/(jagData$seFish.r^2), 
     orgTau = 1/(jagData$seFish.o^2), 
     repEffect_FT = jagData$fis.r,
     nSource = length(unique(jagData$source)), 
     #nStudy = length(unique(jagData$authorsTitle.o)), 
     #study = as.factor(as.character(jagData$authorsTitle.o)), 
     source = as.numeric(as.factor(as.character(jagData$source))))


stanModelCode <- '
data {
  int<lower=0> n; // number of effects included 
  int<lower=0> nSource; // number of replicaiton projects included 
  real y[n]; // estimated treatment effects
  real<lower=0> sigma[n]; // s.e. of effect estimates 
  real orgEffect_FT[n]; // 
  real repEffect_FT[n]; //
  real<lower=0> repTau; //
  real<lower=0> orgTau; //
  int<lower=0> source; //
}
parameters {
  int y<lower=0, upper=1> clust; // binary for replication or not
  real<lower=0> tau; // tau
  real<lower=0, upper=1> phi; // estimated replication prob
  real trueOrgEffect; // est true orig
  real mu_source_ori[nSource]; // est random effect for source
  real<lower=0, upper=1> alpha
}
model {
  for(i in 1:n){
    clust[i] ~ dbern(phi) //# extract errors in variables (FT stands for Fisher-transformed):
    orgEffect_FT[i] ~ normal(trueOrgEffect[i] +  mu_source_ori[source[i]], orgTau[i]) // #  + mu_study_ori[study[i]] 
    repEffect_FT[i] ~ normal(trueRepEffect[i] + mu_source_rep[source[i]] ,  repTau[i]) //# + mu_study_rep[study[i]] 
    trueRepEffect[i] ~ normal(mu[i], tau)
    
    //# if clust[i] = 0 then H0 is true; if clust[i] = 1 then H1 is true and
    //# the replication effect is a function of the original effect, plus a random effect for study:
    mu[i] <- (alpha * (trueOrgEffect[i]) * equals(clust[i], 1))
    //# when clust[i] = 0, then mu[i] = 0;
    //# when clust[i] = 1, then mu[i] = alpha * trueOrgEffect[i]
  }
}
'

stan('Analysis/stanmodel1.stan', data = stanData)



model{
  # Mixture Model Priors:
  alpha ~ dunif(0,1) # flat prior on slope for predicted effect size under H1
  tau ~ dgamma(0.001,0.001) # vague prior on study precision
  phi ~ dbeta(1, 1) # flat prior on the true effect rate
  
  # for(i in 1:nStudy){
  # mu_study_ori[i]~ dnorm(0, 0.000001) # hyperparameter for random intercepts Mean on the study level   ## CHECK 
  # mu_study_rep[i]~ dnorm(0, 0.000001) # hyperparameter for random intercepts Mean on the study level   ## CHECK
  # }
  
  for(i in 1:nSource){
    mu_source_ori[i] ~ dnorm(0, 0.0001) # hyperparameter for random intercepts Mean on the source level (i.e., the OSC project/LOOPR etc.) ## CHECK 
    mu_source_rep[i] ~ dnorm(0, 0.0001) # hyperparameter for random intercepts Mean on the source level (i.e., the OSC project etc.) ## CHECK 
  }
  # prior on true effect size of original studies:
  for (i in 1:n){
    trueOrgEffect[i] ~ dnorm(0, 1)
  }
  
  # Mixture Model Likelihood:
  # Study level
  for(i in 1:n){
    clust[i] ~ dbern(phi)# extract errors in variables (FT stands for Fisher-transformed):
    orgEffect_FT[i] ~ dnorm(trueOrgEffect[i] +  mu_source_ori[source[i]], orgTau[i]) #  + mu_study_ori[study[i]] 
    repEffect_FT[i] ~ dnorm(trueRepEffect[i] + mu_source_rep[source[i]] ,  repTau[i]) # + mu_study_rep[study[i]] 
    trueRepEffect[i] ~ dnorm(mu[i], tau)
    
    # if clust[i] = 0 then H0 is true; if clust[i] = 1 then H1 is true and
    # the replication effect is a function of the original effect, plus a random effect for study:
    mu[i] <- (alpha * (trueOrgEffect[i]) * equals(clust[i], 1))
    # when clust[i] = 0, then mu[i] = 0;
    # when clust[i] = 1, then mu[i] = alpha * trueOrgEffect[i]
  }
}















functions { 
} 
data { 
  int<lower=1> N;  // total number of observations 
  vector[N] Y;  // response variable 
  vector<lower=0>[N] se;  // known sampling error 
  real<lower=0> sigma;  // residual SD 
  // data for group-level effects of ID 1
  int<lower=1> J_1[N];
  int<lower=1> N_1;
  int<lower=1> M_1;
  vector[N] Z_1_1;
  int prior_only;  // should the likelihood be ignored? 
} 
transformed data { 
  vector<lower=0>[N] se2 = square(se); 
} 
parameters { 
  real temp_Intercept;  // temporary intercept 
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // unscaled group-level effects
} 
transformed parameters { 
  // group-level effects 
  vector[N_1] r_1_1 = sd_1[1] * (z_1[1]);
} 
model { 
  vector[N] mu = temp_Intercept + rep_vector(0, N);
  for (n in 1:N) { 
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
  } 
  // priors including all constants 
  target += student_t_lpdf(temp_Intercept | 3, 0, 10); 
  target += student_t_lpdf(sd_1 | 3, 0, 10)
  - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += normal_lpdf(z_1[1] | 0, 1);
  // likelihood including all constants 
  if (!prior_only) { 
    target += normal_lpdf(Y | mu, se);
  } 
} 
generated quantities { 
  // actual population-level intercept 
  real b_Intercept = temp_Intercept; 
} 