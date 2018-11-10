model{
  # Mixture Model Priors:
  tau ~ dgamma(0.001,0.001) # vague prior on study precision
  phi ~ dbeta(1, 1) # flat prior on the true effect rate
  
  for(i in 1:n){
    alpha[i] ~ dunif(0,1) # flat prior on slope for predicted effect size under H1 for each replication project
  }
  
  # prior on true effect size of original studies:
  for (i in 1:n){
    trueOrgEffect[i] ~ dnorm(0, 1) # Normal prior on the original effect size 
  }
  
  # Mixture Model Likelihood:
  # Study level
  for(i in 1:n){
    clust[i] ~ dbern(phi)# extract errors in variables:
    orgEffect[i] ~ dnorm(trueOrgEffect[i] , orgTau[i]) # 
    repEffect[i] ~ dnorm(trueRepEffect[i] ,   repTau[i]) 
    trueRepEffect[i] ~ dnorm(mu[i], tau)
    
    # if clust[i] = 0 then H0 is true; if clust[i] = 1 then H1 is true and
    # the observed replication effect is a function of the original effect:
    mu[i] <- (alpha[i] * trueOrgEffect[i] * equals(clust[i], 1))
    # when clust[i] = 0, then mu[i] = 0;
    # when clust[i] = 1, then mu[i] = alpha * trueOrgEffect[i]
  }
}