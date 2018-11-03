# source(file = 'Data/data_collection_cleaning.R')
source(file = "R code effect size estimation/appendixCodeFunctionsJeffreys.R")
source(file = "Analysis/Wilson score interval.R")
# source(file = 'Data/data_collection_cleaning.R') # this calls the data
library(readr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

sumStats <- function(x, na.rm = T){
return(list(mean = mean(x, na.rm = na.rm), sd = sd(x, na.rm = na.rm), median = median(x, na.rm = na.rm) , quantile = quantile(x, na.rm = na.rm), n = sum(!is.na(x)), nNA = sum(is.na(x))))
}

tableAverageDecrease <- allData %>%
  group_by(source) %>%
  dplyr::summarise(mean=mean(fisherZDiff, na.rm=T), sd=sd(fisherZDiff, na.rm=T))

# Could approximate SEs using: # Not valid for F stats w/ df1 > 2 / chi square tests
allData$seFishAprox.r  <- ifelse( is.na(allData$seFish.r ), sqrt(1/(allData$n.r-3)), allData$seFish.r)
allData$seFishAprox.o  <- ifelse( is.na(allData$seFish.o ), sqrt(1/(allData$n.o-3)), allData$seFish.o)

# first off simple calculation of the proportion of studies conducted per published result given a significnat main effect
# How to get to .9 or .75 of the literature being sig, assuming all studies are statistically significant
studiesPerPublishedPaper <- paste(round(.75/.44,2), "to", round(.9/.44,2))

# converting to cohen's d
effsizes.r <- compute.es::res(allData$correlation.r, n = allData$n.r, verbose = F)
allData$cohensD.r <- effsizes.r$d

# converting to cohen's d 
effsizes.o <- compute.es::res(allData$correlation.o, n = allData$n.o, verbose = F)
allData$cohensD.o <- effsizes.o$d

# Sample characteristics
# Calculating the number included in the meta-analysis (also the number included in the questions asking about change)
nMeta <- sum((!is.na(allData$fisherZDiff)))
nMetaIndividualPapers <- length(unique(allData$authorsTitle.o[!is.na(allData$fisherZDiff)]))

nEffectSizeFell <- sum(allData$correlationDifference.ro < 0, na.rm = T)
propEffectSizeFell <- mean(allData$correlationDifference.ro < 0, na.rm = T)

CIonProp <- WilsonBinCI(sum(!is.na(allData$correlationDifference.ro < 0)), p = propEffectSizeFell)

# N with valid SEs
nSEsValid <- sum(!is.na(allData$seFish.o) & !is.na(allData$seFish.r))

# sum stats effect size (corr and d)
summaryCor.r <- sumStats(allData$correlation.o, na.rm = T)
summaryCor.o <- sumStats(allData$correlation.o, na.rm = T)
summaryD.r <-   sumStats(allData$cohensD.r, na.rm = T)
summaryD.o <-   sumStats(allData$cohensD.o, na.rm = T)

#### Default bayes factors ####
# Setting up functions for BFs for apply function 
# From Wagenmakers, E. J., Verhagen, J., & Ly, A. (2016). How to quantify the evidence for the absence of a correlation. Behav Res Methods, 48(2), 413-426. doi:10.3758/s13428-015-0593-0

bfapply10 <- function(x){
  if(sum(is.na(x)>0)) {
    return(NA)
    }
     else{
    bf10JeffreysIntegrate( r = as.numeric(x[1]),n = as.numeric(x[2]))
     }
  }

bfapplyplus0 <- function(x){
  if(sum(is.na(x)>0)) {
    return(NA)
  }
  else{
    bfPlus0JeffreysIntegrate( r = as.numeric(x[1]),n = as.numeric(x[2]))
  }
}

bfapplyRep0 <- function(x){
  if(sum(is.na(x))>0) {
    return(NA)
  }
  else{
    tryCatch(repBfR0(rOri=as.numeric(x[1]), nOri = as.numeric(x[2]),  rRep = as.numeric(x[3]), nRep = as.numeric(x[4])), error=function(err) NA)
  }
}

# Two sided default Bayes factor
allData$bf10 <- apply(data.frame(allData$correlation.r, allData$n.r), MARGIN = 1, FUN = bfapply10) 
allData$bf01 <- 1/allData$bf10 

# one sided ~ probably makes more sense as they have all been shifted to positive direction
allData$bfplus0 <- apply(data.frame(allData$correlation.r, allData$n.r), MARGIN = 1, FUN = bfapplyplus0) 
allData$bf0plus <- 1/allData$bfplus0 

# Replication Bayes factor
allData$bfRep0 <- apply(data.frame(allData$correlation.o, allData$n.o, allData$correlation.r, allData$n.r), MARGIN = 1, FUN = bfapplyRep0) 
allData$bf0Rep <- 1/allData$bfRep0

allData$bf10.o <-  apply(data.frame(allData$correlation.o, allData$n.o), MARGIN = 1, FUN = bfapply10)
allData$bfRep0_updated <-  allData$bf10.o * allData$bf10

# mean(allData$bfRep0_updated - allData$bfRep0, na.rm =T)
# plot(log(allData$bfRep0_updated), log(allData$bfRep0))
# ## Finding minimium effect that would have been significant in the original study - 
minimumEffectDetectableZ <- qnorm(.05, mean = 0, sd = allData$seFishAprox.o, lower.tail = FALSE)
CIs90UB <-   allData$fis.r + (qnorm(0.95)*allData$seFishAprox.r)
CIs90LB <-   allData$fis.r - (qnorm(0.95)*allData$seFishAprox.r)

upper95.r <- allData$fis.r + ( 1.96 * allData$seFish.r )
lower95.r <- allData$fis.r - ( 1.96 * allData$seFish.r )

# Extracting all studies 
# View(allData[lower95.r > 0 & !allData$significantSameDirection.r,]);  View(allData[lower95.r < 0 & allData$significantSameDirection.r,]); 
nonMatchesStatisticalSig <- sum( c(lower95.r > 0 & !allData$significantSameDirection.r,lower95.r < 0 & allData$significantSameDirection.r), na.rm = T )

tableBayesFactors <-data_frame(Article = natSci$Authors, Original_r = allData$correlation.o[allData$abrev=="natSci"], Original_N =allData$n.o[allData$abrev=="natSci"], Replication_r = allData$correlation.r[allData$abrev=="natSci"], Replication_n = allData$n.r[allData$abrev=="natSci"], Camerer_et_al._BFP0 = natSci$bfP0_rp, Camerer_et_al._BFRep0 = natSci$bfR0_rp, BFrep0 =  1/allData$bf0Rep[allData$abrev=="natSci"],  BF0plus = 1/allData$bf0plus[allData$abrev=="natSci"], BF01 = 1/allData$bf01[allData$abrev=="natSci"])

minimumEffectDetectableR <- ztor(minimumEffectDetectableZ)

# calculating the number of equiv studies
proportionEquiv <- mean(CIs90UB < minimumEffectDetectableZ, na.rm =T)
nEquiv <- sum(CIs90UB < minimumEffectDetectableZ, na.rm =T)

allData$statisticallyEquiv.ro <- CIs90UB < minimumEffectDetectableZ

#Tested with:
#allData$statisticallyEquiv.ro <- (CIs90UB < minimumEffectDetectableZ)
#i = 8
#TOSTER::TOSTr(n = allData$n.r[i], r = allData$correlation.r[i], low_eqbound_r = -1, high_eqbound_r = minimumEffectDetectableR[i])

######################################
#### Amount of change in subsets #####
######################################
# Reduction in all studies: 
reductionAllData <- allData %>% 
  filter(!is.na(fisherZDiff)) %>% 
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T), 
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),    
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(!is.na(fisherZDiff)), nValid = sum(!is.na(allData$fis.r - allData$fis.o)))
# extracting amount of change in subsets
# Sig results replication only (same direction)
reductionSignifcantR <- allData %>% 
  filter(significantSameDirection.r & !is.na(significantSameDirection.r)) %>% 
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T), 
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),    
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(!is.na(significantSameDirection.r & !is.na(significantSameDirection.r))), nValid = sum(!is.na(allData$significantSameDirection.r) & !is.na(allData$fis.r - allData$fis.o)))
# BF01 < 3 (exclude those with moderate or greater evidence for the null)
reductionSbf01LessThan3 <- allData %>% 
  filter(bf01<3) %>% 
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),   
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(!is.na(bf01<3)), nValid = sum(!is.na(allData$bf01) & !is.na(allData$fis.r - allData$fis.o)))
# BF10 > 3 (Only those with evidence for the alternative)
reductionSbf10MoreThan3 <- allData %>% 
  filter(bf10>3) %>% 
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) , 
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),    
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(!is.na(bf10>3)), nValid = sum(!is.na(allData$bf10) & !is.na(allData$fis.r - allData$fis.o)))
# BF0Plus < 3 (exclude those with moderate evidence for the null, one sided alternative)
reductionSbf0PlusLessThan3 <- allData %>% 
  filter(bf0plus<3) %>% 
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),   
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(!is.na(bf0plus<3)), nValid = sum(!is.na(allData$bf0plus)& !is.na(allData$fis.r - allData$fis.o)))
# BFPlus0 (Exclude those without moderate evidence for the one sided alternative)
reductionSbfPlus0MoreThan3 <- allData %>% 
  filter(bfplus0>3) %>%
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),    
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(!is.na(bfplus0>3)), nValid = sum(!is.na(allData$bfplus0) & !is.na(allData$fis.r - allData$fis.o)))
# BF0Rep < 3 (exclude those with moderate evidence for the null, replication alternative)
reductionSbf0RepLessThan3 <- allData %>% 
  filter(bf0Rep<3) %>% 
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),  
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(!is.na(bf0Rep<3)), nValid = sum(!is.na(allData$bf0Rep) & !is.na(allData$fis.r - allData$fis.o)))
# BFRep0 (Exclude those without moderate evidence for the replication alternative)
reductionSbfRep0MoreThan3 <- allData %>% 
  filter(bfRep0>3) %>%
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),  
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(!is.na(bfRep0>3)), nValid = sum(!is.na(allData$bfRep0) & !is.na(allData$fis.r - allData$fis.o)))
# Significant TOST (Exclude those without moderate evidence for the one sided alternative)
reductionEquiv <- allData %>% 
  filter(!statisticallyEquiv.ro) %>%
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),   
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(!is.na(statisticallyEquiv.ro)), nValid = sum(!is.na(allData$statisticallyEquiv.ro)& !is.na(allData$fis.r - allData$fis.o)))


# Bringing all of the above together:
tableReductions <- rbind(Overall = reductionAllData, StatisticalSignificance = reductionSignifcantR, Nonequivalence = reductionEquiv, BF0RepBelow3 = reductionSbf0RepLessThan3 , BFRep0Above3 = reductionSbfRep0MoreThan3,  BF01Below3 = reductionSbf01LessThan3 , BF10Above3 = reductionSbf10MoreThan3, BF0PBelow3 = reductionSbf0PlusLessThan3, BFP0Above3 = reductionSbfPlus0MoreThan3)

# Calculating naive CIs around mean differences 
boundDist <- qt(0.975, df = tableReductions$nTrue - 1)*tableReductions$sdDiff/sqrt(tableReductions$nTrue)
tableReductions$lbMeanDiff <- tableReductions$meanDiff - boundDist
tableReductions$ubMeanDiff <- tableReductions$meanDiff + boundDist

# Converting into Rs 
tableReductions[,-which(str_detect(names(tableReductions), "Prop|nTrue|nValid"))] <- ztor(tableReductions[,-which(str_detect(names(tableReductions), "Prop|nTrue|nValid"))])
# Renaming 
names(tableReductions) <- c("Mean proportion change", "Mean replication ES", "Mean original ES", "Mean ES difference", 
  "SD difference", "Median proportion change", "Median replicaiton ES", "Median original ES", "Median ES difference",
  "n included", "n criteria calculable for", "95% CI LB Mean ES Change", "95% CI UB Mean ES Change")

# Reordering cols
tableReductions <- tableReductions[c("n included", "n criteria calculable for",  "Mean original ES", "Median original ES",  "Mean replication ES", "Median replicaiton ES", 
                  "Mean ES difference",  "95% CI LB Mean ES Change", "95% CI UB Mean ES Change",  "Median ES difference", "SD difference", "Mean proportion change", "Median proportion change")]

kableReductions <- kable(tableReductions, digits = 2)

#########################################################
#### Plots of the effects plotted against each other ####
#########################################################

# Plotting effects on each other 
plotAllData <- ggplot(allData, aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) + geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = .5)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("All data")

# Plotting only non-equivalent studies  
plotNonequiv <- ggplot(allData[!allData$statisticallyEquiv.ro,], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) + geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = .5)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Non-equivalent studies")

# Plotting just significant studies (sig and in same direction)
plotSigR <- ggplot(allData[allData$significantSameDirection.r==TRUE,], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source)))+ geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = .5)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Significant replications in the same direction only")

# two sided test, Exluding studies with BF10 lower than 3, i.e., without evidence moderate or greater for the alternative
plotBF10Greater3 <- ggplot(allData[(allData$bf10>3),], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source)))+ geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = .5)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Exluding studies with BF10 lower than 3")

# two sided test, excluding those with evidence for null > 3
plotBF01Lesser3 <- ggplot(allData[(allData$bf01<3) & (sign(allData$correlation.o) == sign(allData$correlation.r)),], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) + geom_abline( slope = 1, intercept = 0)+  geom_point(na.rm = T, alpha = .5)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Exluding studies with BF01 greater than 3")

# One sided default bayes factor, excluding studies with ev for null > 3
plotBF0plusLesser3 <- ggplot(allData[allData$bf0plus<3,], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) + geom_abline( slope = 1, intercept = 0)+  geom_point(na.rm = T, alpha = .5)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Exluding studies with BF0+ (one sided) greater than 3")

# One sided default bayes factor, including only studies with ev for alternative > 3
plotBFPlus0Greater3 <- ggplot(allData[allData$bfplus0>3,], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) + geom_abline( slope = 1, intercept = 0)+  geom_point(na.rm = T, alpha = .5)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Exluding studies with BF+0 less than than 3")

# Replication bayes factor, including only studies *without* evidence for the null vs the original effect over 3
plotBF0RepLesser3 <- ggplot(allData[allData$bf0Rep <3,], aes(correlation.o, correlation.r, size = log(n.r), colour = as.factor(source))) + geom_abline( slope = 1, intercept = 0)+  geom_point(na.rm = T, alpha = .5)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17))+ ggtitle("Exluding studies with BF0Rep greater than than 3")

# Replication bayes factor, including only studies with evidence for the alternative of greater than 3
plotBFRep0Lesser3 <- ggplot(allData[allData$bfRep0>3,], aes(correlation.o, correlation.r, size = log(n.r), colour = as.factor(source))) + geom_abline( slope = 1, intercept = 0)+  geom_point(na.rm = T, alpha = .5)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Exluding studies with BFRep0 less than than 3")



### Descriptives
nNotSig.r <- sum(allData$significant.r)


######################################
###### Multilevel meta-analysis ######
######################################

# Random effects model with random effects for authors nested within source
REMod <- rma.mv(yi = allData$fisherZDiff, V = allData$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o,  data = allData)
 # cdREMOD <- cooks.distance(REMod, parallel = "snow", ncpus = parallel::detectCores())

# Empirical Bayes estimates for random Effects 


# Model output summary function plus converting results into an understandable format function. 
modelOutputSummary <- function(REMod) {
zDecrease <- as.numeric(REMod[1])
zDecreaseCILB <- as.numeric(REMod[6])
zDecreaseCIUB <- as.numeric(REMod[7])
rDecrease <- ztor(as.numeric(REMod[1]))
rDecreaseCILB <- ztor(as.numeric(REMod[6]))
rDecreaseCIUB <- ztor(as.numeric(REMod[7]))
cDcreaseCIs <- c(rDecreaseCILB, rDecreaseCIUB)
# And in cohen's d: 
dDecrease <- (2*rDecrease) / sqrt(1-rDecrease^2)
dDecreaseCIs <- (2*cDcreaseCIs) / sqrt(1-cDcreaseCIs^2)
randEffREModSource <- ranef(REMod)[1]
list(estimate_Z = zDecrease, CI95_Z = c(zDecreaseCILB, zDecreaseCIUB) , 
     estimate_Cor = rDecrease, CI95_Cor = c(rDecreaseCILB, rDecreaseCIUB), 
     estimate_D = dDecrease, CI95_d = dDecreaseCIs, 
     sourceEstimate = randEffREModSource)
}

REModSum <- modelOutputSummary(REMod)
  
# The first model but with only significant replications
REModOnlySigR <- rma.mv(yi = fisherZDiff, V = allData[allData$significantSameDirection.r==TRUE,]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = allData[allData$significantSameDirection.r==TRUE,])

# the first model removing all studies with BFs0Plus > 3 (i.e., moderate evidence for a null)
REModBF0PlusGreaterThan3Excluded <- rma.mv(yi = fisherZDiff, V = allData[allData$bf0plus < 3 & !is.na(allData$bf0plus < 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = allData[allData$bf0plus < 3 & !is.na(allData$bf0plus > 3),])

# the first model removing all studies with BFsPlus0 < 3 (i.e., those without evidence for the alternative)
REModBFPlus0LessThan3Excluded <- rma.mv(yi = fisherZDiff, V = allData[allData$bfplus0 > 3 & !is.na(allData$bfplus0 > 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = allData[allData$bfplus0 > 3 & !is.na(allData$bfplus0 > 3),])

# the first model removing all studies with BFs01 > 3 (i.e., moderate evidence for a null)
REModBF01GreaterThan3Excluded <- rma.mv(yi = fisherZDiff, V = allData[allData$bf01 < 3 & !is.na(allData$bf01 < 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = allData[allData$bf01 < 3 & !is.na(allData$bf01 > 3),])

# the first model removing all studies with BFs10 < 3 (i.e., those without evidence for the alternative)
REModBF10LessThan3Excluded <- rma.mv(yi = fisherZDiff, V = allData[allData$bf10 > 3 & !is.na(allData$bf10 > 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = allData[allData$bf10 > 3 & !is.na(allData$bf10 > 3),])

# the first model removing all studies with BFs0Rep < 3 (i.e., those without evidence for the for the null vs. origianl greater than 3)
REModBF0RepGreaterThan3Excluded <-  rma.mv(yi = fisherZDiff, V = allData[allData$bf0Rep < 3 & !is.na(allData$bf0Rep < 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = allData[allData$bf0Rep < 3 & !is.na(allData$bf0Rep < 3),])

# the first model removing all studies with BFs0Rep < 3 (i.e., those without evidence for the for the null vs. origianl greater than 3)
REModBFRep0LessThan3Excluded <- rma.mv(yi = fisherZDiff, V = allData[allData$bfRep0 > 3 & !is.na(allData$bfRep0 > 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data =  allData[allData$bfRep0 > 3 & !is.na(allData$bfRep0 > 3),])

# The first model but with only non-equiv 
REModNonequiv <- rma.mv(yi = fisherZDiff, V = allData[allData$statisticallyEquiv.ro==FALSE & !is.na(allData$statisticallyEquiv.ro==FALSE),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = allData[allData$statisticallyEquiv.ro==FALSE & !is.na(allData$statisticallyEquiv.ro==FALSE),])

modRes <- data.frame( modelN = REMod$k, modelEstimate = REMod$b, MLM95lb = REMod$ci.lb, MLM95ub = REMod$ci.ub, row.names = "Overall")
modRes1 <- data.frame(modelN = REModOnlySigR$k, modelEstimate = REModOnlySigR$b, MLM95lb = REModOnlySigR$ci.lb, MLM95ub = REModOnlySigR$ci.ub, row.names = "StatisticalSignificance")
modRes2 <- data.frame(modelN = REModBF0PlusGreaterThan3Excluded$k, modelEstimate = REModBF0PlusGreaterThan3Excluded$b, MLM95lb = REModBF0PlusGreaterThan3Excluded$ci.lb, MLM95ub = REModBF0PlusGreaterThan3Excluded$ci.ub, row.names = "BF0PBelow3")
modRes3 <- data.frame(modelN = REModBFPlus0LessThan3Excluded$k, modelEstimate = REModBFPlus0LessThan3Excluded$b, MLM95lb = REModBFPlus0LessThan3Excluded$ci.lb, MLM95ub = REModBFPlus0LessThan3Excluded$ci.ub, row.names = "BFP0Above3")
modRes4 <- data.frame(modelN = REModBF01GreaterThan3Excluded$k, modelEstimate = REModBF01GreaterThan3Excluded$b, MLM95lb = REModBF01GreaterThan3Excluded$ci.lb, MLM95ub = REModBF01GreaterThan3Excluded$ci.ub, row.names = "BF01Below3")
modRes5 <- data.frame(modelN = REModBF10LessThan3Excluded$k, modelEstimate = REModBF10LessThan3Excluded$b, MLM95lb = REModBF10LessThan3Excluded$ci.lb, MLM95ub = REModBF10LessThan3Excluded$ci.ub, row.names = "BF10Above3")
modRes6 <- data.frame(modelN = REModBF0RepGreaterThan3Excluded$k, modelEstimate = REModBF0RepGreaterThan3Excluded$b, MLM95lb = REModBF0RepGreaterThan3Excluded$ci.lb, MLM95ub = REModBF0RepGreaterThan3Excluded$ci.ub, row.names = "BF0RepBelow3")
modRes7 <- data.frame(modelN = REModBFRep0LessThan3Excluded$k, modelEstimate = REModBFRep0LessThan3Excluded$b, MLM95lb = REModBFRep0LessThan3Excluded$ci.lb, MLM95ub = REModBFRep0LessThan3Excluded$ci.ub, row.names = "BFRep0Above3")
modRes8 <- data.frame(modelN = REModNonequiv$k, modelEstimate = REModNonequiv$b, MLM95lb = REModNonequiv$ci.lb, MLM95ub = REModNonequiv$ci.ub, row.names = "Nonequivalence")

modSumaries <- rbind(modRes, modRes1, modRes2, modRes3, modRes4, modRes5, modRes6, modRes7, modRes8)

tableAllEstimates <- merge.data.frame(tableReductions, modSumaries, by = "row.names", sort = F)



#### Leave one out cross validation for main model, excluding first each source
# # # Commented out to avoid having to run this each time I knit
# LOOTracking <- list()
# # First initial models
# tempData <- allData
# LOOTracking[[8]] <-  rma.mv(yi = tempData$fisherZDiff, V = tempData$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o,  data = tempData)
# LOOTracking[[1]] <-    rma.mv(yi = fisherZDiff, V = tempData[tempData$significantSameDirection.r==TRUE,]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$significantSameDirection.r==TRUE & !is.na(tempData$significantSameDirection.r==TRUE),])
# LOOTracking[[2]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bf0plus < 3 & !is.na(tempData$bf0plus < 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bf0plus < 3 & !is.na(tempData$bf0plus > 3),])
# LOOTracking[[3]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bfplus0 > 3 & !is.na(tempData$bfplus0 > 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bfplus0 > 3 & !is.na(tempData$bfplus0 > 3),])
# LOOTracking[[4]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bf01 < 3 & !is.na(tempData$bf01 < 3),]$seDifference.ro^2,       random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bf01 < 3 & !is.na(tempData$bf01 > 3),])
# LOOTracking[[5]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bf10 > 3 & !is.na(tempData$bf10 > 3),]$seDifference.ro^2,       random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bf10 > 3 & !is.na(tempData$bf10 > 3),])
# LOOTracking[[6]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bf0Rep < 3 & !is.na(tempData$bf0Rep < 3),]$seDifference.ro^2,   random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bf0Rep < 3 & !is.na(tempData$bf0Rep < 3),])
# LOOTracking[[7]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bfRep0 > 3 & !is.na(tempData$bfRep0 > 3),]$seDifference.ro^2,   random =  ~ 1|source/authorsTitle.o, data =  tempData[tempData$bfRep0 > 3 & !is.na(tempData$bfRep0 > 3),])
# 
# 
#  for(i in 1:length(unique(allData$source))) {
#    # exlude <- unique(allData$source)[i]
#    tempData <- allData[-which(allData$source == unique(allData$source)[i]),]
#    LOOTracking[[8+i]] <-  rma.mv(yi = tempData$fisherZDiff, V = tempData$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o,  data = tempData)
#    LOOTracking[[8+i+length(unique(allData$source))]] <-    rma.mv(yi = fisherZDiff, V = tempData[tempData$significantSameDirection.r==TRUE,]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$significantSameDirection.r==TRUE & !is.na(tempData$significantSameDirection.r==TRUE),])
#    LOOTracking[[8+i+2*length(unique(allData$source))]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bf0plus < 3 & !is.na(tempData$bf0plus < 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bf0plus < 3 & !is.na(tempData$bf0plus > 3),])
#    LOOTracking[[8+i+3*length(unique(allData$source))]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bfplus0 > 3 & !is.na(tempData$bfplus0 > 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bfplus0 > 3 & !is.na(tempData$bfplus0 > 3),])
#    LOOTracking[[8+i+4*length(unique(allData$source))]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bf01 < 3 & !is.na(tempData$bf01 < 3),]$seDifference.ro^2,       random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bf01 < 3 & !is.na(tempData$bf01 > 3),])
#    LOOTracking[[8+i+5*length(unique(allData$source))]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bf10 > 3 & !is.na(tempData$bf10 > 3),]$seDifference.ro^2,       random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bf10 > 3 & !is.na(tempData$bf10 > 3),])
#    LOOTracking[[8+i+6*length(unique(allData$source))]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bf0Rep < 3 & !is.na(tempData$bf0Rep < 3),]$seDifference.ro^2,   random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bf0Rep < 3 & !is.na(tempData$bf0Rep < 3),])
#    LOOTracking[[8+i+7*length(unique(allData$source))]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bfRep0 > 3 & !is.na(tempData$bfRep0 > 3),]$seDifference.ro^2,   random =  ~ 1|source/authorsTitle.o, data =  tempData[tempData$bfRep0 > 3 & !is.na(tempData$bfRep0 > 3),])
# }
#  # extrating estiamtes from list
# LOOProject <- sapply(LOOTracking, function(m) m[c(1,5, 105)], simplify = T)
# 
# LOOProjectDF <- data.frame(as.data.frame(t(LOOProject)))
# LOOProjectDF <- data.frame(b = as.numeric(LOOProjectDF$b), p = as.numeric(LOOProjectDF$pval), call = as.character(LOOProjectDF$call))

#  write.csv( LOOProjectDF, "Data/LOOProject.csv", row.names = F)
LOOProjectDF <- read_csv("Data/LOOProject.csv")

LOOProjectDFSum <- LOOProjectDF %>%
  group_by(call) %>%
  dplyr::summarise('proporiton significant'= mean(p<.05), 'Minimum b' = quantile(b)[1], 'Quantile b 25%' = quantile(b)[2], 'Median b' = quantile(b)[3], 'Quantile b 75%' = quantile(b)[4], 'Maximum b' = quantile(b)[5])

maxDiffLOOProject <- max(LOOProjectDFSum[,7] - LOOProjectDFSum[,3])




#### LOO removing studies - this takes ~ 20 mins to run, it was run once and then saved ####
#  for(i in 1:length(unique(allData$authorsTitle.o))) {
#    # exlude <- unique(allData$source)[i]
#    tempData <- allData[-which(allData$authorsTitle.o == unique(allData$authorsTitle.o)[i]),]
#    LOOTracking[[i]] <- rma.mv(yi = tempData$fisherZDiff, V = tempData$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o,  data = tempData)
#    LOOTracking[[i+length(unique(allData$authorsTitle.o))]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$significantSameDirection.r==TRUE,]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$significantSameDirection.r==TRUE,])
#    LOOTracking[[i+ (2*length(unique(allData$authorsTitle.o)))]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bf0plus < 3 & !is.na(tempData$bf0plus < 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bf0plus < 3 & !is.na(tempData$bf0plus > 3),])
#    LOOTracking[[i+ (3*length(unique(allData$authorsTitle.o)))]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bfplus0 > 3 & !is.na(tempData$bfplus0 > 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bfplus0 > 3 & !is.na(tempData$bfplus0 > 3),])
#    LOOTracking[[i+ (4*length(unique(allData$authorsTitle.o)))]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bf01 < 3 & !is.na(tempData$bf01 < 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bf01 < 3 & !is.na(tempData$bf01 > 3),])
#    LOOTracking[[i+ (5*length(unique(allData$authorsTitle.o)))]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bf10 > 3 & !is.na(tempData$bf10 > 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bf10 > 3 & !is.na(tempData$bf10 > 3),])
#    LOOTracking[[i+ (6*length(unique(allData$authorsTitle.o)))]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bf0Rep < 3 & !is.na(tempData$bf0Rep < 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = tempData[tempData$bf0Rep < 3 & !is.na(tempData$bf0Rep < 3),])
#    LOOTracking[[i+ (7*length(unique(allData$authorsTitle.o)))]] <-  rma.mv(yi = fisherZDiff, V = tempData[tempData$bfRep0 > 3 & !is.na(tempData$bfRep0 > 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data =  tempData[tempData$bfRep0 > 3 & !is.na(tempData$bfRep0 > 3),])
#  }

#LOOStudy <- sapply(LOOTracking, function(m) m[c(1,5, 105)], simplify = T)
#LOOStudyDF <- data.frame(as.data.frame(t(LOOStudy)))
#LOOStudyDF <- data.frame(b = as.numeric(LOOStudyDF$b), p = as.numeric(LOOStudyDF$pval), call = as.character(LOOStudyDF$call))

 LOOStudyDF <- read.csv("Data/LOOStudyDF.csv")
 
LOOStudyDFSum <- LOOStudyDF %>%
  group_by(call) %>%
  dplyr::summarise('proporiton significant'= mean(p<.05), 'Minimum b' = quantile(b)[1], 'Quantile b 25%' = quantile(b)[2], 'Median b' = quantile(b)[3], 'Quantile b 75%' = quantile(b)[4], 'Maximum b' = quantile(b)[5])

maxDiffLOOStudy <- max(LOOStudyDFSum[,7] - LOOStudyDFSum[,3])

  maximumStudyDif <- max(as.numeric(REMod[1]) - LOOStudyDF[1])
  
# inf <- cooks.distance.rma.mv(REModBF01GreaterThan3Excluded)

#####################
###### Figures ######
#####################

# functions used in plots 
lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)

# Fishers Z change 
pdf(file = "Figures/ViolinPlotZscoreChange.pdf", width = 5)
ggplot(allData,aes(y = fisherZDiff,x= source, 
                   fill = source, colour = source)) + geom_hline(yintercept = 0, alpha = .1) + 
  geom_flat_violin(width=1.35,
                   position = position_nudge(x = .2, y = 0), alpha = .8) + 
  geom_boxplot(width = .25, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  geom_point(position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  coord_flip() + 
  theme_classic() +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") + 
  theme(legend.position =  'none', axis.text.y = element_text(size=8), 
        axis.title.y = element_blank()) + ylab("Differences in Fisher's Z transformed effect sizes")
dev.off()

pdf(file = "Figures/ViolinPlotPercentageChange.pdf")

# percentage change
ggplot(allData, aes(y = percentageChangeES.ro,x= source, 
                    fill = source, colour = source))+ geom_hline(yintercept = 0, alpha = .1) + 
  geom_flat_violin(width=1.45,
    position = position_nudge(x = .2, y = 0), alpha = .8) + 
  geom_boxplot(width = .3, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  geom_point(position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  coord_flip() + 
  theme_classic() +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") + 
  theme(legend.position =  'none', axis.text.y = element_text(size=8), 
        axis.title.y = element_blank()) + ylab("Percentage change in effect sizes from original to replicaiton")
dev.off()


# Correlation change 
pdf(file = "Figures/ViolinPlotCorrelationDifference.pdf", width = 5)
ggplot(allData,aes(y = correlationDifference.ro,x= source, 
                   fill = source, colour = source)) + geom_hline(yintercept = 0, alpha = .1) + 
  geom_flat_violin(width=1.35,
                   position = position_nudge(x = .2, y = 0), alpha = .8) + 
  geom_boxplot(width = .25, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  geom_point(position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  coord_flip() + 
  theme_classic() +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") + 
  theme(legend.position =  'none', axis.text.y = element_text(size=8), 
        axis.title.y = element_blank()) + ylab("Differences in effect size (correlations)") 
dev.off()



###### Simulation data analysis #######
tableAllEstimatesSim <- read.csv(file = "Data/SimulationModelOutput.csv")

tableAllEstimatesSim$squaredError <- (-tableAllEstimatesSim$propAttenuation-tableAllEstimatesSim$meanPropChange)^2
tableAllEstimatesSim$absoluteError <- abs(-tableAllEstimatesSim$propAttenuation-tableAllEstimatesSim$meanPropChange)

simulationSum <-as.tibble(tableAllEstimatesSim) %>%
  group_by(propNull, propAttenuation, Row.names) %>%
  dplyr::summarise(meanMeanPropChange = mean(meanPropChange), sdMeanPropChange = sd(meanPropChange), meanChange = mean(meanDiff), 
                   meanModelEstimate = mean(modelEstimate), sdModelEstimate = sd(modelEstimate),  nSims = n(), MSE = mean(squaredError), 
                   RMSE = sqrt(mean(squaredError)),  MAE = mean(absoluteError), meanTrueDiff = mean(trueMeanDifference, na.rm = T),
                   meanTrueDiff = mean(trueMeanDifferenceNo0s, na.rm = T), errorMeanPropReduction = mean(-propAttenuation - meanPropChange))

vis <- simulationSum 

# # Plotting distance between estimated mean change and true prop
# ggplot(vis, aes(x = propAttenuation, y = propNull)) +
#   geom_raster(aes(fill = errorMeanPropReduction), interpolate=F)  +
#   scale_fill_gradient2(low="navy", mid="white", high="red", 
#                        midpoint=0)+ facet_wrap(~ Row.names) + theme_bw()
# 
# # Plotting SD by method at different values 
# ggplot(vis, aes(x = propAttenuation, y = propNull)) +
#   geom_raster(aes(fill = sdMeanPropChange), interpolate=F)  +
#   scale_fill_gradient2(low="navy", mid="white", high="red", 
#                        midpoint=0)+ facet_wrap(~ Row.names) + theme_bw()
# 
# # Plotting MSE by method at different values 
# ggplot(vis, aes(x = propAttenuation, y = propNull)) +
#   geom_raster(aes(fill = MSE), interpolate=F)  +
#   scale_fill_gradient2(low="navy", mid="white", high="red", 
#                        midpoint=0)+ facet_wrap(~ Row.names) + theme_bw()
# 
# # Plotting MAE by method at different values 
# ggplot(vis, aes(x = propAttenuation, y = propNull)) +
#   geom_raster(aes(fill = MAE), interpolate=F)  +
#   scale_fill_gradient2(low="navy", mid="white", high="red", 
#                        midpoint=0)+ facet_wrap(~ Row.names) + theme_bw()
#
simulationSumByType <- as.tibble(tableAllEstimatesSim) %>%
  group_by(Row.names) %>%
  dplyr::summarise(meanMeanPropChange = mean(meanPropChange, na.rm = T), sdMeanPropChange = sd(meanPropChange, na.rm = T), 
                   meanModelEstimate = mean(modelEstimate, na.rm = T), sdModelEstimate = sd(modelEstimate, na.rm = T),  nSims = n(), MSE = mean(squaredError, na.rm = T), 
                   RMSE = sqrt(mean(squaredError, na.rm = T)),  MAE = mean(absoluteError, na.rm = T))
#View(simulationSumByType)
#
#
#
simulationSumByTypeLessThan75 <- as.tibble(tableAllEstimatesSim) %>%
  group_by(Row.names, below.8s = propNull < .8 & propAttenuation < .8) %>%
  dplyr::summarise(meanMeanPropChange = mean(meanPropChange, na.rm = T), sdMeanPropChange = sd(meanPropChange, na.rm = T), 
                   meanModelEstimate = mean(modelEstimate, na.rm = T), sdModelEstimate = sd(modelEstimate, na.rm = T),  nSims = n(), MSE = mean(squaredError, na.rm = T), 
                   RMSE = sqrt(mean(squaredError, na.rm = T)),  MAE = mean(absoluteError, na.rm = T))
#View(simulationSumByTypeLessThan75[simulationSumByTypeLessThan75$`propNull < 0.8 & propAttenuation < 0.8`==T,])
#View(simulationSumByTypeLessThan75)
#
#
## Plotting mean distance between estimated mean change and true prop
#ggplot(vis, aes(x = propAttenuation, y = propNull)) +
#  geom_raster(aes(fill = meanModelEstimate/mean(allData$fis.o, na.rm =T)), interpolate=F)  +
#  scale_fill_gradient2(low="navy", mid="white", high="red", 
#                       midpoint=0)+ facet_wrap(~ Row.names) + theme_bw()
#
#
## Plotting mean distance between estimated mean change and true prop
#ggplot(vis, aes(x = propAttenuation, y = propNull)) +
#  geom_raster(aes(fill = meanMeanPropChange), interpolate=F)  +
#  scale_fill_gradient2(low="navy", mid="white", high="red", 
#                       midpoint=0, limits = c(-1,0))+ facet_wrap(~ Row.names) + theme_bw()
#
#

