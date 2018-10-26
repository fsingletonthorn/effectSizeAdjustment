data {
  int<lower=0> n; // number of effects included 
  int<lower=0> nSource; // number of replicaiton projects included 
  int<lower=1> nStudy; // num groups
  real orgEffect_FT[n]; // 
  real repEffect_FT[n]; //
  real<lower=0> repTau; //
  real<lower=0> orgTau; //
  // int<lower=1,upper=[nSource]> source; //
  // int<lower=1,upper=[nStudy]> study; //

}
parameters {
  // int cluster[n]; // prior binary for replication or not
  real<lower=0> tau; // prior scale, tau
  real<lower=0, upper=1> phi; // estimated replication prob
  real trueOrgEffect[n]; // est true orig
  // real mu_source_ori[nSource]; // est random effect for source
  real<lower=0, upper=1> alpha; // alpha
}
model {
for(i in 1:n) {
  int<lower=0, upper=1> cluster[i]; // ~ binomial(1, phi); // extract errors in variables 
  orgEffect_FT[i] ~ normal(trueOrgEffect[i] , orgTau[i]); //   + mu_study_ori[study[i]] +  mu_source_ori[source[i]]
  repEffect_FT[i] ~ normal(trueRepEffect[i] ,  repTau[i]); // + mu_study_rep[study[i]]  + mu_source_rep[source[i]] 
  trueRepEffect[i] ~ normal(mu[i], tau);
    
    // if clust[i] = 0 then H0 is true; if clust[i] = 1 then H1 is true and
    // the replication effect is a function of the original effect, plus a random effect for study:
    mu[i] <- (alpha * (trueOrgEffect[i]) * equals(clust[i], 1));
    // when clust = 0, then mu[i] = 0;
    // when clust = 1, then mu[i] = alpha * trueOrgEffect
}
}
