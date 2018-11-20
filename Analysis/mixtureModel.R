library(rjags)
# https://osf.io/xhj4d/  - Camerer, C. F., Dreber, A., Holzmeister, F., Ho, T.-H., Huber, J., Johannesson, M., . . . Wu, H. (2018). Evaluating the replicability of social science experiments in Nature and Science between 2010 and 2015. Nature Human Behaviour, 2(9), 637-644. doi:10.1038/s41562-018-0399-z
# Getting rid of missing data

jagData <- allData[ !is.na(allData$fis.o) & !is.na(allData$fis.r)  & !is.na(allData$n.o)  & !is.na(allData$n.r) ,]


# Parameters to keep
params <- c(#"mu",
  "phi",
  "clust",
  "H1original",
  "H0True",
  "H1decrease",
  #  "repEffect" ,
  "trueRepEffect",
  "alpha")

## Filling in any missing SEs 
jagData$seFish.o <- sqrt(1/(jagData$n.o-3))
jagData$seFish.r <- sqrt(1/(jagData$n.r-3))

# Setting up the model 
## NOTE: JAGS uses the precision, not the variance, when specifying a normal distribution (precision = 1/variance), i.e., 1/(se^2)
jagMod <- jags.model(file = 'Analysis/BMWMod.R', 
                     data = list(orgEffect = jagData$fis.o, n = length(jagData$fis.o), repTau = 1/(jagData$seFish.r^2), orgTau = 1/(jagData$seFish.o^2), repEffect = jagData$fis.r),
                     n.chains=4)

# Running model and summarising 
samples <- coda.samples(jagMod,params,n.iter = 100000)
samplesSum <- summary(samples)

## Collecting information from both models
sums <- do.call(cbind.data.frame, samplesSum)
statsMeans <- sums$statistics.Mean
# par(mar = c(0, 0, 0, 0))
alpha <- statsMeans[grepl(pattern = "alpha", rownames(sums))]
phi <- statsMeans[grepl(pattern = "phi", rownames(sums))]
jagData$probRealEffect <- statsMeans[grepl(pattern = "clust", rownames(sums))]
jagData$trueRepEffect <- statsMeans[grepl(pattern = "trueRepEffect", rownames(sums))]

mean((jagData$trueRepEffect-jagData$fis.o)/jagData$fis.o)

# Highest prob density interval on alpha 
intervalSamples <- HPDinterval(samples)
chains <- data.frame(t(sapply(intervalSamples, FUN = function(x) x[1,])))
HDIs <- data.frame(mean(chains$lower), mean(chains$upper))



# Setting up the model 
## NOTE: JAGS uses the precision, not the variance, when specifying a normal distribution (precision = 1/variance), i.e., 1/(se^2)
jagMod <- jags.model(file = 'Analysis/BMWMod_3cat_additional.R', 
                     data = list(orgEffect = jagData$fis.o, n = length(jagData$fis.o), repTau = 1/(jagData$seFish.r^2), orgTau = 1/(jagData$seFish.o^2), repEffect = jagData$fis.r),
                     n.chains=4)

# Running model and summarising 
samples <- coda.samples(jagMod,params,n.iter = 10000)
samplesSum <- summary(samples)

sums <- do.call(cbind.data.frame, samplesSum)
statsMeans <- sums$statistics.Mean
# par(mar = c(0, 0, 0, 0))
alpha <- statsMeans[grepl(pattern = "alpha", rownames(sums))]
# Probability of null, decrease, original effect 
phi <- statsMeans[grepl(pattern = "phi", rownames(sums))]

jagData$probH1Original <- statsMeans[grepl(pattern = "H1original", rownames(sums))]
jagData$probH1decrease <- statsMeans[grepl(pattern = "H1decrease", rownames(sums))]
jagData$probH0True <- statsMeans[grepl(pattern = "H0True", rownames(sums))]

jagData$mostLikelyClass <- apply(jagData[c("probH0True","probH1decrease","probH1Original")], 1,which.max)
jagData$probabilityMostLikelyClass <- apply(jagData[c("probH0True","probH1decrease","probH1Original")], 1,max)
jagData$mostLikelyClass <- factor(jagData$mostLikelyClass, labels = c("H0 true", "H1 attenuated", "H1 original effect"))


# Setting up the model 
## NOTE: JAGS uses the precision, not the variance, when specifying a normal distribution (precision = 1/variance), i.e., 1/(se^2)
jagMod <- jags.model(file = 'Analysis/BMWMod_vanilla_plus.R', 
                     data = list(orgEffect = jagData$fis.o, n = length(jagData$fis.o), repTau = 1/(jagData$seFish.r^2), repEffect = jagData$fis.r),
                     n.chains=4)

# Running model and summarising 
samples <- coda.samples(jagMod,params,n.iter = 10000)
samplesSum <- summary(samples)

## Collecting information from both models
sums <- do.call(cbind.data.frame, samplesSum)
statsMeans <- sums$statistics.Mean
# par(mar = c(0, 0, 0, 0))
alpha <- statsMeans[grepl(pattern = "alpha", rownames(sums))]
phi <- statsMeans[grepl(pattern = "phi", rownames(sums))]
jagData$probRealEffect <- statsMeans[grepl(pattern = "clust", rownames(sums))]

mean((jagData$trueRepEffect-jagData$fis.o)/jagData$fis.o)

# Highest prob density interval on alpha 
intervalSamples <- HPDinterval(samples)
chains <- data.frame(t(sapply(intervalSamples, FUN = function(x) x[1,])))
HDIs <- data.frame(mean(chains$lower), mean(chains$upper))


# Setting up the model 
## NOTE: JAGS uses the precision, not the variance, when specifying a normal distribution (precision = 1/variance), i.e., 1/(se^2)
jagMod_trunc <- jags.model(file = 'Analysis/BMWMod_truncated.R', 
                     data = list(orgEffect = jagData$fis.o, n = length(jagData$fis.o), repTau = 1/(jagData$seFish.r^2), orgTau = 1/(jagData$seFish.o^2), repEffect = jagData$fis.r),
                     n.chains=4)

# Running model and summarising 
samples_trunc <- coda.samples(jagMod_trunc,params,n.iter = 10000)
samplesSum_trunc <- summary(samples_trunc)

## Collecting information from both models
sums_trunc <- do.call(cbind.data.frame, samplesSum_trunc)
statsMeans <- sums_trunc$statistics.Mean
# par(mar = c(0, 0, 0, 0))
alpha <- statsMeans[grepl(pattern = "alpha", rownames(sums_trunc))]
phi <- statsMeans[grepl(pattern = "phi", rownames(sums_trunc))]
jagData$probRealEffect <- statsMeans[grepl(pattern = "clust", rownames(sums_trunc))]

# Highest prob density interval on alpha 
intervalSamples <- HPDinterval(samples)
chains <- data.frame(t(sapply(intervalSamples, FUN = function(x) x[1,])))
HDIs <- data.frame(mean(chains$lower), mean(chains$upper))




### Saving these for future access
# saveRDS(alpha, "Data/mixtureModelOutput/alphaSimple.rds")
# saveRDS(phi, "Data/mixtureModelOutput/phiSimple.rds")
# saveRDS(HDIs, "Data/mixtureModelOutput/HDIsAlphaSimple.rds") # this was used to save the data that was used later for analysis 
# write.csv(jagData, "Data/mixtureModelOutput/jagData.csv", row.names = F)

## Diganostics:
#rejectionRate(samples)
#effectiveSize(samples)
# plot(samples)
# Traceplot - be warned, this takes forever unless you are only looking at one parameter at a time - I suggest alpha, the most important (assuming all else is working)
# traceplot(samples)


### additional model - this accounts for 3 categories - null true effects, attenuated effects, and an additional set of studies which are not attenuated

jagMod_additional3 <- jags.model(file = 'Analysis/mixtureModel3Cat.R', 
                     data = list(orgEffect = jagData$fis.o, n = length(jagData$fis.o), repTau = 1/(jagData$seFish.r^2), orgTau = 1/(jagData$seFish.o^2), repEffect = jagData$fis.r),
                                 n.chains=4)

# Parameters to keep
params <- c("phi",
  # "clust",
  "H1original",
  "H0True",
  "H1decrease",
  "trueRepEffect",
  "alpha")

# Running model and summarising 
samples3 <- coda.samples(jagMod_additional3,params,n.iter = 20000)
samples3Sum<-summary(samples3)
samples3Sum

boundsamples3 <- as.mcmc(do.call(rbind,samples3)) 
HPDSamples3 <- HPDinterval(boundsamples3)
HPDCat <- HPDSamples3[grepl(pattern = "alpha", rownames(HPDSamples3))]


sums <- do.call(cbind.data.frame, samples3Sum)
statsMeans <- sums$statistics.Mean
# par(mar = c(0, 0, 0, 0))
alpha <- statsMeans[grepl(pattern = "alpha", rownames(sums))]
# Probability of null, decrease, original effect 
phi <- statsMeans[grepl(pattern = "phi", rownames(sums))]

jagData$probH1Original <- statsMeans[grepl(pattern = "H1original", rownames(sums))]
jagData$probH1decrease <- statsMeans[grepl(pattern = "H1decrease", rownames(sums))]
jagData$probH0True <- statsMeans[grepl(pattern = "H0True", rownames(sums))]

jagData$trueRepEffect <- statsMeans[grepl(pattern = "trueRepEffect", rownames(sums))]

jagData$mostLikelyClass <- apply(jagData[c("probH0True","probH1decrease","probH1Original")], 1,which.max)
jagData$probabilityMostLikelyClass <- apply(jagData[c("probH0True","probH1decrease","probH1Original")], 1,max)
jagData$mostLikelyClass <- factor(jagData$mostLikelyClass, labels = c("H0 true", "H1 attenuated", "H1 original effect"))










jagMod_additional <- jags.model(file = 'Analysis/BMWMod_additional.R', 
                                data = list(orgEffect = jagData$fis.o, n = length(jagData$fis.o), repTau = 1/(jagData$seFish.r^2), orgTau = 1/(jagData$seFish.o^2), repEffect = jagData$fis.r,
                                            nSource = length(unique(jagData$source)), #nStudy = length(unique(jagData$authorsTitle.o)), 
                                            #study = as.factor(as.character(jagData$authorsTitle.o)), 
                                            source = as.factor(as.character(jagData$source))),
                                n.chains=4)


# Parameters to keep
params <- c("metaAlpha",
            "alphaTau",
  "phi",
  # "clust",
  # "orgEffect" ,
  # "trueOrgEffect",
  # "trueRepEffect",
  "alpha")

# Running model and summarising 
samples2 <- coda.samples(jagMod_additional,params,n.iter = 20000)
samples2Sum<-summary(samples2)
samples2Sum

boundsamples2 <- as.mcmc(do.call(rbind,samples2)) 
# 
hdp <- HPDinterval(boundsamples2)
# mean(samples2Sum$statistics[9:313,2])
# plot(samples2)
# View(samples2Sum)

# HPDinterval(samples2)

#### Plots etc
## Collecting information from both models
sums <- do.call(cbind.data.frame, samples2Sum)
statsMeans <- sums$statistics.Mean
# par(mar = c(0, 0, 0, 0))
alpha <- statsMeans[grepl(pattern = "alpha", rownames(sums))]
metaAlpha <- statsMeans[grepl(pattern = "metaAlpha", rownames(sums))] 
phi <- statsMeans[grepl(pattern = "phi", rownames(sums))]
jagData$probRealEffect <- statsMeans[grepl(pattern = "clust", rownames(sums))]
jagData$trueRepEffect <- statsMeans[grepl(pattern = "trueRepEffect", rownames(sums))]




write.csv(jagData, "Data/mixtureModelOutput/jagData.csv", row.names = F)
