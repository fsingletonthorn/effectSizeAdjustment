library(rjags)
library(tidyverse)
# Source these two before running to load data:
# source(file = 'Data/data_collection_cleaning.R')
# source('simplifiedAnalysisScript.R')
# https://osf.io/xhj4d/  - Camerer, C. F., Dreber, A., Holzmeister, F., Ho, T.-H., Huber, J., Johannesson, M., . . . Wu, H. (2018). Evaluating the replicability of social science experiments in Nature and Science between 2010 and 2015. Nature Human Behaviour, 2(9), 637-644. doi:10.1038/s41562-018-0399-z
# Getting rid of missing data

jagData <- allData[ !is.na(allData$fis.o) & !is.na(allData$fis.r)  & !is.na(allData$n.o)  & !is.na(allData$n.r) ,]

# Parameters to keep
params <- c("phi",
            "clust",
            "altHypothesis",
            "alpha") # Remove alpha 

## Filling in any missing SEs 
jagData$seFish.o <- sqrt(1/(jagData$n.o-3))
jagData$seFish.r <- sqrt(1/(jagData$n.r-3))

# Setting up the model 
## NOTE: JAGS uses the precision, not the variance, when specifying a normal distribution (precision = 1/variance), i.e., 1/(se^2)
jagMod <- jags.model(file = 'Analysis/simpleMix.R', 
                     data = list(orgEffect = jagData$fis.o, n = length(jagData$fis.o), repTau = 1/(jagData$seFish.r^2), orgTau = 1/(jagData$seFish.o^2), repEffect = jagData$fis.r),
                     n.chains=4)

# Running model and summarising 
samples <- coda.samples(jagMod,params,n.iter = 50000)
samplesSum <- summary(samples)

## Collecting information from both models
sums <- do.call(cbind.data.frame, samplesSum)
statsMeans <- sums$statistics.Mean

phi <- statsMeans[grepl(pattern = "phi", rownames(sums))]
alpha <- statsMeans[grepl(pattern = "alpha", rownames(sums))] # Remove 
jagData$probRealEffect <- statsMeans[grepl(pattern = "clust", rownames(sums))]
jagData$trueRepEffect <- statsMeans[grepl(pattern = "altHypothesis", rownames(sums))]

mean((jagData$trueRepEffect-jagData$fis.o)/jagData$fis.o)

# Highest prob density interval on alpha 
# Binding to treat as one MCMC chain 
samplesDechained <- as.mcmc( do.call(rbind, samples))

intervalSamples <- HPDinterval(samplesDechained)

# alphaHPD <- intervalSamples[row.names(intervalSamples) == 'alpha']
phiHPD <- intervalSamples[row.names(intervalSamples) == 'phi']

trueRep <- samplesDechained[,which(str_detect(colnames(samplesDechained), "trueRep"))]
clust <- samplesDechained[,which(str_detect(colnames(samplesDechained), "clust"))]

propBelow.1 <- sum(samplesDechained[,which(str_detect(colnames(samplesDechained), "trueRep"))] < .10036)/ (nrow((samplesDechained[,which(str_detect(colnames(samplesDechained), "trueRep"))] >= .1)) * ncol(samplesDechained[,which(str_detect(colnames(samplesDechained), "trueRep"))]))

### Saving these for future access commented out for the sake of avoiding irritation
# saveRDS(alpha, "Data/mixtureModelOutput/alphaSimple.rds")
# saveRDS(phi, "Data/mixtureModelOutput/phiSimple.rds")
# saveRDS(alphaHPD, "Data/mixtureModelOutput/HDIsAlphaSimple.rds") # this was used to save the data that was used later for analysis 
# saveRDS(phiHPD, "Data/mixtureModelOutput/HPDphiSimple.rds")
# write.csv(jagData, "Data/mixtureModelOutput/jagData.csv", row.names = F)
# saveRDS(propBelow.1.BMM, "Data/mixtureModelOutput/ValuesBelow.1.rds")

ggplot(jagData, aes(x = fis.o, y = fis.r,  color = probRealEffect, size = n.r)) +
    geom_abline( slope = 1, intercept = 0) + geom_point(alpha = .8, na.rm = T)+ theme_classic() +
    guides(color = guide_legend(title = "Posterior\nassignment\nrate"),
                       shape = guide_legend(title = "True effect size < 0.1"),
                       size = guide_legend(title = "Replication\nSample size",
                                          values= trans_format("identity", function(x) round(exp(x),0)), order = 2)) +
    scale_size(trans = "log", breaks = c(10, 100, 1000, 10000)) + geom_point(colour = "black", na.rm = T, size = .5, shape = 3) +
    xlab("Original correlation")+ ylab("Replication correlation") + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) 
