library(rjags)
library(tidyverse)
# Source these two before running to load data:
# source(file = 'Data/data_collection_cleaning.R')
# https://osf.io/xhj4d/  - Camerer, C. F., Dreber, A., Holzmeister, F., Ho, T.-H., Huber, J., Johannesson, M., . . . Wu, H. (2018). Evaluating the replicability of social science experiments in Nature and Science between 2010 and 2015. Nature Human Behaviour, 2(9), 637-644. doi:10.1038/s41562-018-0399-z
# Getting rid of missing data

jagData <- allData[ !is.na(allData$fis.o) & !is.na(allData$fis.r)  & !is.na(allData$n.o)  & !is.na(allData$n.r) ,]

# Parameters to keep
params <- c("phi",
  "clust",
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
# Binding to treat as one MCMC chain 
samplesDechained <- as.mcmc( do.call(rbind, samples))

intervalSamples <- HPDinterval(samplesDechained)

alphaHPD <- intervalSamples[row.names(intervalSamples) == 'alpha']
phiHPD <- intervalSamples[row.names(intervalSamples) == 'phi']

trueRep <- samplesDechained[,which(str_detect(colnames(samplesDechained), "trueRep"))]
clust <- samplesDechained[,which(str_detect(colnames(samplesDechained), "clust"))]

# meltedClust <- as.mcmc(stack(clust))
# meltedtrueRep <- as.mcmc(stack(trueRep))

propBelow.1.BMM <- sum(samplesDechained[,which(str_detect(colnames(samplesDechained), "trueRep"))] < .10036)/ (nrow((samplesDechained[,which(str_detect(colnames(samplesDechained), "trueRep"))] >= .1)) * ncol(samplesDechained[,which(str_detect(colnames(samplesDechained), "trueRep"))]))
propBelow.1THA.BMM <- sum(samplesDechained[,which(str_detect(colnames(samplesDechained), "trueRep"))] < .10036)/ (nrow((samplesDechained[,which(str_detect(colnames(samplesDechained), "trueRep"))] >= .1)) * ncol(samplesDechained[,which(str_detect(colnames(samplesDechained), "trueRep"))]))

### Saving these for future access commented out for the sake of avoiding irritation
# saveRDS(alpha, "Data/mixtureModelOutput/alphaSimple.rds")
# saveRDS(phi, "Data/mixtureModelOutput/phiSimple.rds")
# saveRDS(alphaHPD, "Data/mixtureModelOutput/HDIsAlphaSimple.rds") # this was used to save the data that was used later for analysis 
# saveRDS(phiHPD, "Data/mixtureModelOutput/HPDphiSimple.rds")
#  write.csv(jagData, "Data/mixtureModelOutput/jagData.csv", row.names = F)
#  saveRDS(propBelow.1.BMM, "Data/mixtureModelOutput/ValuesBelow.1.rds")
  