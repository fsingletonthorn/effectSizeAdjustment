library(rstan)
library(brms)

### Need to run data cleaning first 

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# study as factor
allData$study <- as.factor(allData$authorsTitle.o)
# Setting up same MLM as the main analysis but in bayesian framework
mod <- brm(fisherZDiff | se(seDifference.ro) ~ 1 + (1|source/study), data = allData[!is.na(allData$seDifference.ro) & !is.na(allData$fisherZDiff),], cores = parallel::detectCores(), 
           control =list(adapt_delta = .999, max_treedepth = 30), warmup = 1000, iter = 3000)
summary(mod)
REsMod <- ranef(mod)

pairs(mod)

stancode(mod)
 
modFixed <- brm(fisherZDiff | se(seDifference.ro) ~ 1 + source + (1|study), data = allData[!is.na(allData$seDifference.ro) & !is.na(allData$fisherZDiff),], cores = parallel::detectCores(), 
           control =list(adapt_delta = .999, max_treedepth = 30), warmup = 500, iter = 1500)
summary(modFixed)



modSimple <- brm(fisherZDiff | se(seDifference.ro) ~ 1 + (1|source/study), data = allData[!is.na(allData$seDifference.ro) & !is.na(allData$fisherZDiff),], cores = parallel::detectCores(), 
           control =list(adapt_delta = .95, max_treedepth = 15), iter = 1000)





brms::pp_check(mod)
brms::marginal_effects(mod)

mod2 <- brm(fisherZDiff | se(seDifference.ro) ~ 1 + source + (1|study), data = allData[!is.na(allData$seDifference.ro) & !is.na(allData$fisherZDiff),], cores = 4, control =list(adapt_delta = .99, max_treedepth = 15))
REModFixSource
posterior_summary(mod)

summary(mod2)

pairs(mod2)



  
  
  
  
  model{
    # Mixture Model Priors:
     # flat prior on slope for predicted effect size under H1
     # vague prior on study precision
    phi ~ dbeta(1, 1) # flat prior on the true effect rate
    # prior on true effect size of original studies:
    for (i in 1:n){
      trueOrgEffect[i] ~ dnorm(0, 1)
    }
    # Mixture Model Likelihood:
    for(i in 1:n){
      clust[i] ~ dbern(phi)# extract errors in variables (FT stands for Fisher-transformed):
      orgEffect_FT[i] ~ dnorm(trueOrgEffect[i], orgTau[i])
      repEffect_FT[i] ~ dnorm(trueRepEffect[i], repTau[i])
      trueRepEffect[i] ~ dnorm(mu[i], tau)
      # if clust[i] = 0 then H0 is true; if clust[i] = 1 then H1 is true and
      # the replication effect is a function of the original effect:
      mu[i] <- alpha * trueOrgEffect[i] * equals(clust[i], 1)
      # when clust[i] = 0, then mu[i] = 0;
      # when clust[i] = 1, then mu[i] = alpha * trueOrgEffect[i]
    }
  }











