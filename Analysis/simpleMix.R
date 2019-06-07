model{
  # Mixture Model Priors:
  phi ~ dbeta(1, 1) # flat prior on the true effect rate
  # Mixture Model Likelihood:
  for(i in 1:n){
    # if clust[i] = 0 then H0 is true; if clust[i] = 1 then H1 is true and
    # the replication effect is sampled from the original effect plus error variance:
    clust[i] ~ dbern(phi)

    altHypothesis[i] ~ dnorm(orgEffect[i], orgTau[i])T(2.220446e-16, )
    
    # when clust[i] = 0, then mu[i] = 0;
    # when clust[i] = 1, then mu[i] = trueOrgEffect[i]
    trueEffect[i] <- (altHypothesis[i]) * clust[i]
    repEffect[i] ~ dnorm(altHypothesis[i], repTau[i])
  }
}