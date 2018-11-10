library(rjags)
library(brms)

# https://osf.io/xhj4d/  - Camerer, C. F., Dreber, A., Holzmeister, F., Ho, T.-H., Huber, J., Johannesson, M., . . . Wu, H. (2018). Evaluating the replicability of social science experiments in Nature and Science between 2010 and 2015. Nature Human Behaviour, 2(9), 637-644. doi:10.1038/s41562-018-0399-z
# Getting rid of missing data
jagData <- allData[ !is.na(allData$fis.o) & !is.na(allData$fis.r)  & !is.na(allData$n.o)  & !is.na(allData$n.r) ,]

## ADDING IN MISSING DATA 
jagData$seFish.o <- sqrt(1/(jagData$n.o-3))
jagData$seFish.r <- sqrt(1/(jagData$n.r-3))

# Setting up the model 
## NOTE: JAGS uses the precision, not the variance, when specifying a normal distribution (precision = 1/variance), i.e., 1/(se^2)
jagMod <- jags.model(file = 'Analysis/BMWMod.R', 
                     data = list(orgEffect_FT = jagData$fis.o, n = length(jagData$fis.o), repTau = 1/(jagData$seFish.r^2), orgTau = 1/(jagData$seFish.o^2), repEffect_FT = jagData$fis.r),
                     n.chains=4)

# Parameters to keep
params <- c("alpha",
            "mu",
            "phi",
            "clust",
            "trueOrgEffect",
            "trueRepEffect")

# Running model and summarising 
samples <- coda.samples(jagMod,params,n.iter = 10000)
samplesSum <- summary(samples)
# Highest prob density interval
HPDinterval(samples)
plot(samples)
# Traceplot - be warned, this takes forever unless you are only looking at one parameter at a time - I suggest alpha, the most important (assuming all else is working)
# traceplot(samples)

# Taking a burn in of 5000 out

View(samplesSum)
# plot(samples)


### additional model - it should also account for random effects for both study and for source

jagMod_additional <- jags.model(file = 'Analysis/BMWMod_additional.R', 
                     data = list(orgEffect = jagData$fis.o, n = length(jagData$fis.o), repTau = 1/(jagData$seFish.r^2), orgTau = 1/(jagData$seFish.o^2), repEffect = jagData$fis.r,
                                 nSource = length(unique(jagData$source)), #nStudy = length(unique(jagData$authorsTitle.o)), 
                                 #study = as.factor(as.character(jagData$authorsTitle.o)), 
                                 source = as.factor(as.character(jagData$source))),
                                 n.chains=4)

# Parameters to keep
params <- c("mu",
            "phi",
  "clust",
  "orgEffect" ,
  "trueOrgEffect",
  "repEffect" ,
  "trueRepEffect",
  "alpha")

# Running model and summarising 
samples2 <- coda.samples(jagMod_additional,params,n.iter = 10000)
samples2Sum<-summary(samples2)
samples2Sum
mean(samples2Sum$statistics[9:313,2])
# plot(samples2)
# View(samples2Sum)

HPDinterval(samples2)
  
## Collecting information from both models
sums <- do.call(cbind.data.frame, samples2Sum)
statsMeans <- sums$statistics.Mean
# par(mar = c(0, 0, 0, 0))
alpha <- statsMeans[grepl(pattern = "alpha", rownames(sums))]
alpha <- statsMeans[grepl(pattern = "phi", rownames(sums))]
jagData$probRealEffect <- statsMeans[grepl(pattern = "clust", rownames(sums))]
jagData$trueRepEffectModel <- statsMeans[grepl(pattern = "trueRepEffect", rownames(sums))]
jagData$trueOrgEffectModel <- statsMeans[grepl(pattern = "trueOrgEffect", rownames(sums))]

plot(jagData$probRealEffect ~ as.factor(jagData$source))
plot(jagData$trueRepEffectModel ~ as.factor(jagData$source))
plot(jagData$probRealEffect ~ jagData$trueRepEffectModel)
plot(jagData$trueRepEffectModel ~ as.factor(jagData$source))
plot(trueOrgEffectModel ~ fis.o, data = jagData)
plot(trueRepEffectModel ~ fis.r, data = jagData)

sum((jagData$probRealEffect > .5)/nrow(jagData))

