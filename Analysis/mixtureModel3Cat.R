model{
  # Mixture Model Priors:
  tau ~ dgamma(0.001,0.001) # vague prior on study precision
  phi ~ ddirch(mPriorProb) # Flat prior on the model priors
  mPriorProb[1] <- 1
  mPriorProb[2] <- 1
  mPriorProb[3] <- 1
    
    alpha ~ dunif(0,1) # flat prior on attenuation factor for each replication project
  
  # prior on true effect size of original studies:
  for (i in 1:n){
    trueOrgEffect[i] ~ dnorm(0, 1) # Normal prior on the original effect size 
  }
  
  # Mixture Model Likelihood:
  # Study level
  for(i in 1:n){
    clust[i] ~ dcat(phi)# extract errors in variables:
    
    orgEffect[i] ~ dnorm(trueOrgEffect[i] , orgTau[i]) # 
    
    # if clust[i] = 0 then H0 is true; if clust[i] = 1 then the true effect size is a function of the original effect size (times alpha),
    # if phi == 2 the the effect is exactly equal to the original effect 
    # the observed replication effect is a function of the original effect:
    mu[i] <- ifelse(clust[i] == 1,  (alpha * trueOrgEffect[i]), ifelse(clust[i] == 2, trueOrgEffect[i], 0))
        
    
    trueRepEffect[i] ~ dnorm(mu[i], tau)    
    repEffect[i] ~ dnorm(trueRepEffect[i] ,   repTau[i]) 
    
  }
}