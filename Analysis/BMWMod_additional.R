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