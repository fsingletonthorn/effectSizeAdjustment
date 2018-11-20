model{
  # Mixture Model Priors:
  tau ~ dgamma(0.001,0.001) # vague prior on study precision
  phi ~ ddirch(mPriorProb) # Flat prior on the model priors
  mPriorProb[1] <- 1 # This sets the priors to be equal
  mPriorProb[2] <- 1 # This sets the priors to be equal
  mPriorProb[3] <- 1 # This sets the priors to be equal
  
  alpha ~ dunif(0,1) # flat prior on attenuation factor for each replication project
  
  # prior on true effect size of original studies:
  for (i in 1:n){
    trueOrgEffect[i] ~ dnorm(0, 1) # Normal prior on the original effect size 
  }
  
  # Mixture Model Likelihood:
  # Study level
  for(i in 1:n){
    clust[i] ~ dcat(phi)# cluster is equal to one of the categories with probability equal to cat 
    orgEffect[i] ~ dnorm(trueOrgEffect[i] , orgTau[i]) # the original effect is from a dist with a mean equal to the true org effect (estimated) w/ a precision equal to the SD of the org
    
    # if clust[i] = 0 then H0 is true; if clust[i] = 1 then the true effect size is a function of the original effect size (times alpha),
    # if phi == 2 the the effect is exactly equal to the original effect 
    # the observed replication effect is a function of the original effect:
    mu[i] <- ifelse(clust[i] == 2,  (alpha * trueOrgEffect[i]), ifelse(clust[i] == 3, trueOrgEffect[i], 0))
    H1original[i] <-  ifelse(clust[i] == 3, 1, 0)  
    H1decrease[i] <-  ifelse(clust[i] == 2, 1, 0)
    H0True[i] <- ifelse(clust[i] == 1, 1, 0)
    
    repEffect[i] ~ dnorm(mu[i],   repTau[i]) 
    
  }
}