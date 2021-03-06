##      model{
##      # Mixture Model Priors:
##      tau ~ dgamma(0.001,0.001) # vague prior on study precision
##      phi ~ dbeta(1, 1) # flat prior on the true effect rate
##      
##      # prior on alpha, the effect size attenuation value 
##      for(i in 1:nSource){
##        alpha[i] ~ dunif(0,1) # flat prior on attenuation factor for each replication project
##      }
##      
##      # prior on true effect size of original studies:
##      for (i in 1:n){
##      trueOrgEffect[i] ~ dnorm(0, 1) # Normal prior on the original effect size 
##      }
##      
##      # Mixture Model Likelihood:
##      # Study level
##      for(i in 1:n){
##      clust[i] ~ dbern(phi)# extract errors in variables:
##      orgEffect[i] ~ dnorm(trueOrgEffect[i] , orgTau[i]) # 
##      repEffect[i] ~ dnorm(trueRepEffect[i] ,   repTau[i]) 
##      trueRepEffect[i] ~ dnorm(mu[i], tau)
##      
##      # if clust[i] = 0 then H0 is true; if clust[i] = 1 then H1 is true and
##      # the observed replication effect is a function of the original effect:
##      mu[i] <- (alpha[source[i]] * trueOrgEffect[i] * equals(clust[i], 1))
##      # when clust[i] = 0, then mu[i] = 0;
##      # when clust[i] = 1, then mu[i] = alpha * trueOrgEffect[i]
##        }
##      }

model{
  # Mixture Model Priors:
  tau ~ dgamma(0.001,0.001) # vague prior on study precision
  alphaTau ~ dgamma(0.01,0.01) # vague prior on the amount of variability in replication project alphas, more likely that higher values hold that tau
  phi ~ dbeta(1, 1) # flat prior on the true effect rate
  
  # prior on alpha, the effect size attenuation value 
  # Prior on meta-alpha, the mean effect size 
  metaAlpha ~ dunif(0,1)
  # prior on alpha, the effect size attenuation value for each study
  for(i in 1:nSource){
    alpha[i] ~ dnorm(metaAlpha,alphaTau) # T(0,1) # flat prior on attenuation factor for each replication project
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
    trueRepEffect[i] ~ dnorm(mu[i], tau) T(0,)
    
    # if clust[i] = 0 then H0 is true; if clust[i] = 1 then H1 is true and
    # the observed replication effect is a function of the original effect:
    mu[i] <- (alpha[source[i]] * trueOrgEffect[i] * equals(clust[i], 1))
    # when clust[i] = 0, then mu[i] = 0;
    # when clust[i] = 1, then mu[i] = alpha * trueOrgEffect[i]
  }
}