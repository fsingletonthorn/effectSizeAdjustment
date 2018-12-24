# source(file = 'Data/data_collection_cleaning.R')
source(file = "R code effect size estimation/appendixCodeFunctionsJeffreys.R")
source(file = "Analysis/Wilson score interval.R")
# source(file = 'Data/data_collection_cleaning.R') # this sources the data
library(readr)
library(MASS)
library(tidyverse)
library(ggplot2)
library(scales)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(bestNormalize)
library(reshape2)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

# Setting up summary functions:
sumStats <- function(x, na.rm = T){
return(list(mean = mean(x, na.rm = na.rm), sd = sd(x, na.rm = na.rm), median = median(x, na.rm = na.rm) , quantile = quantile(x, na.rm = na.rm), n = sum(!is.na(x)), nNA = sum(is.na(x))))
}

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

niceMLMESum <- function(REMod) {
  data_frame(Estimate = c(REMod$b, rep(NA, 3)), "95% CI LB" = c(REMod$ci.lb, rep(NA, 3)), "95% CI UB" = c(REMod$ci.ub, rep(NA, 3)), SE = c(REMod$se,  rep(NA, 3)), p = c(ifelse(REMod$pval<.001, "< .001", REMod$pval),  rep(NA, 3)), 
             "Random effects" = c(NA, paste0("Project variance = ", round(REMod$sigma2[1], 3), ", n = ", REMod$s.nlevels[1]), paste0("Article variance = ", round(REMod$sigma2[2], 3), ", n = ", REMod$s.nlevels[2]),  paste0("QE(",REMod$k-1, ") = ", round(REMod$QE, 2),  ", p ", ifelse(REMod$QEp <.001, "< .001", paste("=" , round(REMod$QEp, 2))))))
}


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

# Could approximate SEs using: # Not valid for F stats w/ df1 > 2 / chi square tests, used as an approximation where necessary
allData$seFishAprox.r  <- ifelse( is.na(allData$seFish.r ), sqrt(1/(allData$n.r-3)), allData$seFish.r)
allData$seFishAprox.o  <- ifelse( is.na(allData$seFish.o ), sqrt(1/(allData$n.o-3)), allData$seFish.o)

# mean raw decrease
meanDecrease <- mean(allData$percentageChangeES.ro, na.rm = T)

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

# Two sided default Bayes factor
allData$bf10 <- apply(data.frame(allData$correlation.r, allData$n.r), MARGIN = 1, FUN = bfapply10) 
allData$bf01 <- 1/allData$bf10 

# one sided ~ probably makes more sense as they have all been shifted to positive direction
allData$bfplus0 <- apply(data.frame(allData$correlation.r, allData$n.r), MARGIN = 1, FUN = bfapplyplus0) 
allData$bf0plus <- 1/allData$bfplus0 

# Replication Bayes factor
allData$bfRep0 <- apply(data.frame(allData$correlation.o, allData$n.o, allData$correlation.r, allData$n.r), MARGIN = 1, FUN = bfapplyRep0) 
allData$bf0Rep <- 1/allData$bfRep0

# some of these fail because the BF is ~ infinity, setting them to that here
allData$bfRep0[is.na(allData$bfRep0) & (!is.na(allData$correlation.o) & !is.na(allData$n.o) & !is.na(allData$correlation.r) & !is.na(allData$n.r))] <- Inf
allData$bf0Rep[is.na(allData$bf0Rep) & (!is.na(allData$correlation.o) & !is.na(allData$n.o) & !is.na(allData$correlation.r) & !is.na(allData$n.r))] <- 1/Inf

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
            nTrue = sum((allData$significantSameDirection.r & !is.na(allData$fisherZDiff))), nValid = sum(!is.na(allData$significantSameDirection.r) & !is.na(allData$fis.r - allData$fis.o)))
# BF01 < 3 (exclude those with moderate or greater evidence for the null)
reductionSbf01LessThan3 <- allData %>% 
  filter(bf01<3) %>% 
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),   
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(bf01<3), nValid = sum(!is.na(allData$bf01) & !is.na(allData$fis.r - allData$fis.o)))
# BF10 > 3 (Only those with evidence for the alternative)
reductionSbf10MoreThan3 <- allData %>% 
  filter(bf10>3) %>% 
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) , 
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),    
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(bf10>3), nValid = sum(!is.na(allData$bf10) & !is.na(allData$fis.r - allData$fis.o)))
# BF0Plus < 3 (exclude those with moderate evidence for the null, one sided alternative)
reductionSbf0PlusLessThan3 <- allData %>% 
  filter(bf0plus<3) %>% 
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),   
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(bf0plus<3), nValid = sum(!is.na(allData$bf0plus)& !is.na(allData$fis.r - allData$fis.o)))
# BFPlus0 (Exclude those without moderate evidence for the one sided alternative)
reductionSbfPlus0MoreThan3 <- allData %>% 
  filter(bfplus0>3) %>%
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),    
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(bfplus0>3), nValid = sum(!is.na(allData$bfplus0) & !is.na(allData$fis.r - allData$fis.o)))
# BF0Rep < 3 (exclude those with moderate evidence for the null, replication alternative)
reductionSbf0RepLessThan3 <- allData %>% 
  filter(bf0Rep<3) %>% 
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),  
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(bf0Rep<3), nValid = sum(!is.na(allData$bf0Rep) & !is.na(allData$fis.r - allData$fis.o)))
# BFRep0 (Exclude those without moderate evidence for the replication alternative)
reductionSbfRep0MoreThan3 <- allData %>% 
  filter(bfRep0>3) %>%
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),  
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(bfRep0>3), nValid = sum(!is.na(allData$bfRep0) & !is.na(allData$fis.r - allData$fis.o)))
# Significant TOST (Exclude those without moderate evidence for the one sided alternative)
reductionEquiv <- allData %>% 
  filter(!statisticallyEquiv.ro) %>%
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),   
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(!allData$statisticallyEquiv.ro, na.rm = T), nValid = sum(!is.na(allData$statisticallyEquiv.ro)& !is.na(allData$fis.r - allData$fis.o)))

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
plotAllData <- ggplot(allData, aes(correlation.o, correlation.r, size = n.r, colour = as.factor(source))) +
  geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = .5) +
  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + 
  guides( size = guide_legend(title = "Replication\nSample size", values= trans_format("identity", function(x) round(exp(x),0)), order = 2),
          colour = guide_legend(title = "Replication projects"), override.aes = list(alpha = 1), order = 1) +
 xlab("Original correlation")+ ylab("Replication correlation") + scale_size(trans = "log", breaks = c(150, 3000, 60000))

# Plotting only non-equivalent studies  
plotNonequiv <- ggplot(allData[!allData$statisticallyEquiv.ro & !is.na(allData$statisticallyEquiv.ro),], aes(correlation.o, correlation.r, size = n.r, colour = as.factor(source))) +
  geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = .5) +
  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + 
  guides( size = guide_legend(title = "Replication\nSample size", values= trans_format("identity", function(x) round(exp(x),0)), order = 2),
          colour = guide_legend(title = "Replication projects"), override.aes = list(alpha = 1), order = 1) +
  xlab("Original correlation")+ ylab("Replication correlation") + scale_size(trans = "log", breaks = c(150, 3000, 60000))

# Plotting just significant studies (sig and in same direction)
plotSigR <- ggplot(allData[allData$significantSameDirection.r==TRUE,], aes(correlation.o, correlation.r, size = n.r, colour = as.factor(source))) +
  geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = .5) +
  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + 
  guides( size = guide_legend(title = "Replication\nSample size", values= trans_format("identity", function(x) round(exp(x),0)), order = 2),
          colour = guide_legend(title = "Replication projects"), override.aes = list(alpha = 1), order = 1) +
  xlab("Original correlation")+ ylab("Replication correlation") + scale_size(trans = "log", breaks = c(150, 3000, 60000))

# two sided test, Exluding studies with BF10 lower than 3, i.e., without evidence moderate or greater for the alternative
plotBF10Greater3 <- ggplot(allData[(allData$bf10>3) & !is.na(allData$bf10<3),],aes(correlation.o, correlation.r, size = n.r, colour = as.factor(source))) +
  geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = .5) +
  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + 
  guides( size = guide_legend(title = "Replication\nSample size", values= trans_format("identity", function(x) round(exp(x),0)), order = 2),
          colour = guide_legend(title = "Replication projects"), override.aes = list(alpha = 1), order = 1) +
  xlab("Original correlation")+ ylab("Replication correlation") + scale_size(trans = "log", breaks = c(150, 3000, 60000))

# two sided test, excluding those with evidence for null > 3
plotBF01Lesser3 <- ggplot(allData[(allData$bf01<3) & !is.na(allData$bf01<3),], aes(correlation.o, correlation.r, size = n.r, colour = as.factor(source))) +
  geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = .5) +
  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + 
  guides( size = guide_legend(title = "Replication\nSample size", values= trans_format("identity", function(x) round(exp(x),0)), order = 2),
          colour = guide_legend(title = "Replication projects"), override.aes = list(alpha = 1), order = 1) +
  xlab("Original correlation")+ ylab("Replication correlation") + scale_size(trans = "log", breaks = c(150, 3000, 60000))

# One sided default bayes factor, excluding studies with ev for null > 3
plotBF0plusLesser3 <- ggplot(allData[allData$bf0plus<3 & !is.na(allData$bf0plus<3),], aes(correlation.o, correlation.r, size = n.r, colour = as.factor(source))) +
  geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = .5) +
  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + 
  guides( size = guide_legend(title = "Replication\nSample size", values= trans_format("identity", function(x) round(exp(x),0)), order = 2),
          colour = guide_legend(title = "Replication projects"), override.aes = list(alpha = 1), order = 1) +
  xlab("Original correlation")+ ylab("Replication correlation") + scale_size(trans = "log", breaks = c(150, 3000, 60000))

# One sided default bayes factor, including only studies with ev for alternative > 3
plotBFPlus0Greater3 <- ggplot(allData[allData$bfplus0>3 & !is.na(allData$bfplus0>3),], aes(correlation.o, correlation.r, size = n.r, colour = as.factor(source))) +
  geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = .5) +
  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + 
  guides( size = guide_legend(title = "Replication\nSample size", values= trans_format("identity", function(x) round(exp(x),0)), order = 2),
          colour = guide_legend(title = "Replication projects"), override.aes = list(alpha = 1), order = 1) +
  xlab("Original correlation")+ ylab("Replication correlation") + scale_size(trans = "log", breaks = c(150, 3000, 60000))

# Replication bayes factor, including only studies *without* evidence for the null vs the original effect over 3
plotBF0RepLesser3 <- ggplot(allData[allData$bf0Rep <3 & !is.na(allData$bf0Rep <3),], aes(correlation.o, correlation.r, size = n.r, colour = as.factor(source))) +
  geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = .5) +
  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + 
  guides( size = guide_legend(title = "Replication\nSample size", values= trans_format("identity", function(x) round(exp(x),0)), order = 2),
          colour = guide_legend(title = "Replication projects"), override.aes = list(alpha = 1), order = 1) +
  xlab("Original correlation")+ ylab("Replication correlation") + scale_size(trans = "log", breaks = c(150, 3000, 60000))

# Replication bayes factor, including only studies with evidence for the alternative of greater than 3
plotBFRep0Greater3 <- ggplot(allData[allData$bfRep0>3 & !is.na(allData$bfRep0>3),], aes(correlation.o, correlation.r, size = n.r, colour = as.factor(source))) +
  geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = .5) +
  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + 
  guides( size = guide_legend(title = "Replication\nSample size", values= trans_format("identity", function(x) round(exp(x),0)), order = 2),
          colour = guide_legend(title = "Replication projects"), override.aes = list(alpha = 1), order = 1) +
  xlab("Original correlation")+ ylab("Replication correlation") + scale_size(trans = "log", breaks = c(150, 3000, 60000))

### Descriptives
nNotSig.r <- sum(allData$significant.r)

######################################
###### Multilevel meta-analysis ######
######################################

# Random effects model with random effects for authors nested within source
REMod <- rma.mv(yi = allData$fisherZDiff, V = allData$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o,  data = allData)
 # cdREMOD <- cooks.distance(REMod, parallel = "snow", ncpus = parallel::detectCores())

# Normalising replication p values for inclusion in the model
# bestNormalize(allData$cleanedpVal.r) - selected normalising transform, this version used as it is faster:
temp <- orderNorm(allData$cleanedpVal.r, na.rm = T)
allData$normalisedpVal.r  <- transf.arcsin(allData$cleanedpVal.r)


REMod.p.val.tukey <- rma.mv(yi = fisherZDiff, V = seDifference.ro^2 , mod = normalisedpVal.r, random =  ~ 1|source/authorsTitle.o,  data = allData)
REMod.p.val.norm <- rma.mv(yi = fisherZDiff, V = seDifference.ro^2 , mod = temp$x.t, random =  ~ 1|source/authorsTitle.o,  data = allData)
REMod.p.val.cleaned <- rma.mv(yi = fisherZDiff, V = seDifference.ro^2 , mod = cleanedpVal.r, random =  ~ 1|source/authorsTitle.o,  data = allData)

REMod.p.val.cleaned.sum <-   data_frame(" " = c("Estimate", "p value", NA,NA,NA), Estimate = c(REMod.p.val.cleaned$b, rep(NA, 3)), "95% CI LB" = c(REMod.p.val.cleaned$ci.lb, rep(NA, 3)), "95% CI UB" = c(REMod.p.val.cleaned$ci.ub, rep(NA, 3)), SE = c(REMod.p.val.cleaned$se,  rep(NA, 3)), p = c(ifelse(REMod.p.val.cleaned$pval<.001, "< .001", round(REMod.p.val.cleaned$pval, 3)),  rep(NA, 3)), 
                                        "Random effects" = c(NA, NA, paste0("Project variance = ", round(REMod.p.val.cleaned$sigma2[1], 3), ", n = ", REMod.p.val.cleaned$s.nlevels[1]), paste0("Article variance = ", round(REMod.p.val.cleaned$sigma2[2], 3), ", n = ", REMod.p.val.cleaned$s.nlevels[2]),  paste0("QE(",REMod.p.val.cleaned$k-1, ") = ", round(REMod.p.val.cleaned$QE, 2),  ", p ", ifelse(REMod.p.val.cleaned$QEp <.001, "< .001", paste("=" , round(REMod.p.val.cleaned$QEp, 2))))))

# Empirical Bayes estimates for random Effects 
BLUPsSource <- ranef(REMod)[1]


niceREModSum<- niceMLMESum(REMod)

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

# brining these together
modSumaries <- rbind(modRes, modRes1, modRes8, modRes2, modRes3, modRes4, modRes5, modRes6, modRes7)
# Estiating the degree of effect size change as a proportion of the average effect size in psychology 
modSumaries$`Estimated % attenuation` <- (modSumaries$modelEstimate/mean(allData$fis.o, na.rm = T))*100
modSumaries$`LB % attenuation` <- (modSumaries$MLM95lb/mean(allData$fis.o, na.rm = T)*100)
modSumaries$`UB % attenuation` <- (modSumaries$MLM95ub/mean(allData$fis.o, na.rm = T)*100)
modSumariesR <- modSumaries
# converting to z
modSumariesR[2:4] <- ztor(modSumaries[2:4])

niceModelSums <- lapply(X = list("All Data" = REMod, "Only Signficant replications" = REModOnlySigR,  "BF0Plus < 3" = REModBF0PlusGreaterThan3Excluded, 
                                 "BFplus0 > 3" = REModBFPlus0LessThan3Excluded, "BF01 < 3" = REModBF01GreaterThan3Excluded, 
                                 "BF10 > 3" = REModBF10LessThan3Excluded, "BF0rep < 3" = REModBF0RepGreaterThan3Excluded, 
                                 "BFrep0 > 3" = REModBFRep0LessThan3Excluded, "Non-equivalent studies" = REModNonequiv), niceMLMESum)
  
  

tableAllEstimates <- merge.data.frame(tableReductions, modSumaries, by = "row.names", sort = F)


tableAverageDecrease <- allData %>%
  group_by(source) %>%
  dplyr::summarise(mean=mean(fisherZDiff, na.rm=T), sd=sd(fisherZDiff, na.rm=T))



# estimating the probability of obtaining results as extreme given that x studies were removed randomly

#applicData <- allData[!is.na(allData$bfRep0),]
#simulationOutput <- 1000 %>% rerun(
#  applicData[sample( x = 1:nrow(applicData), size = modRes1$modelN,  replace = F),]  %>%
#  rma.mv(yi = fisherZDiff, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = .) %>% 
#  .$b)
#
#mean(unlist(simulationOutput) > as.numeric(REModBFRep0LessThan3Excluded$b))
#

#applicData <- allData[!is.na(allData$bfRep0),]
#simulationOutput <- 1000 %>% rerun(
#  applicData[sample( x = 1:nrow(applicData), size = modRes1$modelN,  replace = F),]  %>%
#  rma.mv(yi = fisherZDiff, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = .) %>% 
#  .$b)
#
#mean(unlist(simulationOutput) > as.numeric(REModBFRep0LessThan3Excluded$b))
#


#### Leave one out cross validation for main model, excluding first each source ####
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
# This reads in the data from the LOO analysis above 
LOOProjectDF <- read_csv("Data/LOOProject.csv")

LOOProjectDFSum <- LOOProjectDF %>%
  group_by(Subsample = call) %>%
  dplyr::summarise('Proportion significant'= mean(p<.05), 'Minimum estimate' = quantile(b)[1],  '25th percentile' = quantile(b)[2], 'Median' = quantile(b)[3], '75th percentile' = quantile(b)[4], 'Maximum estimate' = quantile(b)[5])

maxDiffLOOProject <- max(LOOProjectDFSum[,7] - LOOProjectDFSum[,3])

# Renaming files
namesR <- str_split(str_split( LOOProjectDFSum$Subsample, "tempData\\[tempData\\$", simplify = T)[,2], " ", simplify = T)[,1:3]
namesR[!str_detect(namesR[,2], ">|<"),c(2,3)] <- ""
namesR[str_detect(namesR[,1], "signif"),1] <- "Significant in same direction"
namesR[namesR[,1] == "",1] <- "All studies included"
namesR <- apply(namesR, 1, paste, collapse="")

LOOProjectDFSum$Subsample <- namesR

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
# This reads in the data produced by the above text 
 LOOStudyDF <- read.csv("Data/LOOStudyDF.csv")
 
LOOStudyDFSum <- LOOStudyDF %>%
  group_by(Subsample = call) %>%
  dplyr::summarise('Proportion significant'= mean(p<.05), 'Minimum estimate' = quantile(b)[1],  '25th percentile' = quantile(b)[2], 'Median' = quantile(b)[3], '75th percentile' = quantile(b)[4], 'Maximum estimate' = quantile(b)[5])

# This is just relabeling things for easy presentation 
namesR <- str_split(str_split( LOOStudyDFSum$Subsample, "tempData\\[tempData\\$", simplify = T)[,2], " ", simplify = T)[,1:3]
namesR[!str_detect(namesR[,2], ">|<"),c(2,3)] <- ""
namesR[str_detect(namesR[,1], "signif"),1] <- "Significant in same direction"
namesR[namesR[,1] == "",1] <- "All studies included"
namesR <- apply(namesR, 1, paste, collapse="")

LOOStudyDFSum$Subsample <- namesR

maxDiffLOOStudy <- max(LOOStudyDFSum[,7] - LOOStudyDFSum[,3])

maximumStudyDif <- max(as.numeric(REMod[1]) - LOOStudyDF[1])
  
#####################
###### Figures ######
#####################

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
tableAllEstimatesSim$modelAbsoluteError <- abs(tableAllEstimatesSim$trueMeanDifferenceNo0s - tableAllEstimatesSim$modelEstimate)

simulationSum <-as.tibble(tableAllEstimatesSim) %>%
  group_by(propNull, propAttenuation, Row.names) %>%
  dplyr::summarise(meanChange = mean(meanDiff, na.rm = T), 
                   meanModelEstimate = mean(modelEstimate, na.rm = T), sdModelEstimate = sd(modelEstimate, na.rm = T),  nSims = n(), MSE = mean(squaredError, na.rm = T), 
                   RMSE = sqrt(mean(squaredError, na.rm = T)),  MAE = mean(absoluteError, na.rm = T), meanTrueDiff = mean(trueMeanDifference, na.rm = T),
                   meanTrueDiff = mean(trueMeanDifferenceNo0s, na.rm = T), errorMeanPropReduction = mean(-propAttenuation - meanPropChange, na.rm = T), Error_SD = sd(-propAttenuation - meanPropChange, na.rm = T))

vis <- simulationSum 

# Plotting distance between estimated mean change and true prop
meanErrorPlot <- ggplot(vis, aes(x = propAttenuation, y = propNull)) +
  geom_raster(aes(fill = errorMeanPropReduction), interpolate=F)  +
  scale_fill_gradient2("Mean error \n", low="navy", mid="white", high="red", 
                       midpoint=0)+ facet_wrap(~ Row.names) + theme_bw() +
  labs(x = "Proportion true attenuation", y = "Proportion true null effects")

# Plotting SD by method at different values 
sdErrorPlot <- ggplot(vis, aes(x = propAttenuation, y = propNull)) +
  geom_raster(aes(fill = Error_SD), interpolate=F)  +
  scale_fill_gradient2("Error SD \n", low="navy", mid="white", high="red", 
                       midpoint=0)+ facet_wrap(~ Row.names) + theme_bw() +
  labs(x = "Proportion true attenuation", y = "Proportion true null effects")

# Plotting MAE by method at different values 
MAEPlot <- ggplot(vis, aes(x = propAttenuation, y = propNull)) +
  geom_raster(aes(fill = MAE), interpolate=F)  +
  scale_fill_gradient2("MAE \n", low="navy", mid="white", high="red", 
                       midpoint=0)+ facet_wrap(~ Row.names) + theme_bw() +
  labs(x = "Proportion true attenuation", y = "Proportion true null effects")

simulationSumByType <- as.tibble(tableAllEstimatesSim) %>%
  group_by(Subsample = Row.names) %>%
  dplyr::summarise( # meanMeanPropChange = mean(meanPropChange, na.rm = T), sdMeanPropChange = sd(meanPropChange, na.rm = T), 
    # meanModelEstimate = mean(modelEstimate, na.rm = T), sdModelEstimate = sd(modelEstimate, na.rm = T), 
    nSims = n(), MSE = mean(squaredError, na.rm = T), Mean_Error = mean(-propAttenuation-meanPropChange, na.rm = T),
    RMSE = sqrt(mean(squaredError, na.rm = T)),  MAE = mean(absoluteError, na.rm = T), # AE_SD = sd(absoluteError, na.rm = T), 
    Error_SD = sd(-propAttenuation - meanPropChange, na.rm = T))

simulationSumByTypeLessThan75 <- as.tibble(tableAllEstimatesSim) %>%
  group_by(Subsample = Row.names, below.8s = propNull < .8 & propAttenuation < .8) %>%
  dplyr::summarise(# meanMeanPropChange = mean(meanPropChange, na.rm = T), sdMeanPropChange = sd(meanPropChange, na.rm = T), 
    # meanModelEstimate = mean(modelEstimate, na.rm = T), sdModelEstimate = sd(modelEstimate, na.rm = T), 
    nSims = n(), MSE = mean(squaredError, na.rm = T), 
    RMSE = sqrt(mean(squaredError, na.rm = T)),  MAE = mean(absoluteError, na.rm = T), AE_SD = sd(absoluteError, na.rm = T), 
    Error_SD = sd(-propAttenuation - meanPropChange, na.rm = T)) 

accuracyOutput <- read.csv(file = "Data/SimulationAccuracyCriteria.csv")

simulationAccuracy  <- as.tibble(accuracyOutput) %>%
  group_by(propNull, propAttenuation) %>%
  dplyr::summarise( 
    mean_StatisticallyEqui = mean(StatisticallyEquiv, na.rm = T),
    sd_StatisticallyEquiv = sd(StatisticallyEquiv, na.rm = T), 
    mean_StatisticalSig = mean(StatisticalSig, na.rm = T),    
    sd_StatisticalSig = sd(StatisticalSig, na.rm = T), 
    mean_BF01 = mean(BF01, na.rm = T),              
    sd_BF0 = sd(BF01, na.rm = T),
    mean_BF10 = mean(BF10, na.rm = T),            
    sd_BF10 = sd(BF01, na.rm = T),
    mean_BFplus0 = mean(BFplus0, na.rm = T),           
    sd_BFplus0 = sd(BFplus0, na.rm = T), 
    mean_BF0plus = mean(BF0plus, na.rm = T),           
    sd_BF0plus = sd(BF0plus, na.rm = T), 
    mean_BFrep0 = mean(BFrep0, na.rm = T),            
    sd_BFrep0 = sd(BFrep0, na.rm = T), 
    mean_BF0rep = mean(BF0rep, na.rm = T),            
    sd_BF0rep = sd(BF0rep, na.rm = T), n())

vis <- simulationAccuracy %>% 
  gather(key = "Subsample", value = "Accuracy", names(simulationAccuracy)[str_which(names(simulationAccuracy), "mean")])
vis$Subsample <- str_remove_all(vis$Subsample, "mean_")
vis$Subsample <- ifelse(str_detect(vis$Subsample, "BFrep0|BFplus0|BF10"), paste(vis$Subsample,"> 3"), ifelse(str_detect(vis$Subsample, "BF0rep|BF0plus|BF01"), paste(vis$Subsample,"< 3"), vis$Subsample))
vis$Subsample <- ifelse(vis$Subsample == "StatisticalSig", "Statistical significance", ifelse(vis$Subsample == "StatisticallyEqui", "Not statistically equivalent", vis$Subsample) )


# Plotting mean accuracy 
accuracyPlot<-   ggplot(vis, aes(x = propAttenuation, y = propNull)) +
  geom_raster(aes(fill = Accuracy), interpolate=F)  +
  scale_fill_gradient2("Accuracy \n", low="red", mid="white", high="navy", 
                       midpoint=.5, limits = c(0,1))+ facet_wrap(~ Subsample) + theme_bw() +
  labs(x = "Proportion true attenuation", y = "Proportion true null effects") + theme(legend.position = c(0.825, .11), legend.direction="horizontal")

simulationAccuracyByType  <- as.tibble(accuracyOutput) %>%
  dplyr::summarise( 
    mean_StatisticallyEqui = mean(StatisticallyEquiv, na.rm = T),
    sd_StatisticallyEquiv = sd(StatisticallyEquiv, na.rm = T), 
    mean_StatisticalSig = mean(StatisticalSig, na.rm = T),    
    sd_StatisticalSig = sd(StatisticalSig, na.rm = T), 
    mean_BF01 = mean(BF01, na.rm = T),              
    sd_BF0 = sd(BF01, na.rm = T),
    mean_BF10 = mean(BF10, na.rm = T),            
    sd_BF10 = sd(BF01, na.rm = T),
    mean_BFplus0 = mean(BFplus0, na.rm = T),           
    sd_BFplus0 = sd(BFplus0, na.rm = T), 
    mean_BF0plus = mean(BF0plus, na.rm = T),           
    sd_BF0plus = sd(BF0plus, na.rm = T), 
    mean_BFrep0 = mean(BFrep0, na.rm = T),            
    sd_BFrep0 = sd(BFrep0, na.rm = T), 
    mean_BF0rep = mean(BF0rep, na.rm = T),            
    sd_BF0rep = sd(BF0rep, na.rm = T), n())

names(simulationAccuracyByType) <- str_remove_all(names(simulationAccuracyByType), "mean_")
names(simulationAccuracyByType) <- str_replace_all(names(simulationAccuracyByType), "sd_", "SD ")
names(simulationAccuracyByType) <- ifelse(str_detect(names(simulationAccuracyByType), "BFrep0|BFplus0|BF10"), paste(names(simulationAccuracyByType),"> 3"), ifelse(str_detect(names(simulationAccuracyByType), "BF0rep|BF0plus|BF01"), paste(names(simulationAccuracyByType),"< 3"), names(simulationAccuracyByType)))
names(simulationAccuracyByType) <- ifelse(names(simulationAccuracyByType)== "StatisticalSig", "Statistical significance", ifelse(names(simulationAccuracyByType) == "StatisticallyEqui", "Not statistically equivalent", names(simulationAccuracyByType)) )

simulationAccuracyByTypeDF <- data_frame("Data inclusion rule" = names(simulationAccuracyByType)[-str_which(names(simulationAccuracyByType), 'SD|n\\(\\)')])
simulationAccuracyByTypeDF$Accuracy <- t(simulationAccuracyByType[-str_which(names(simulationAccuracyByType), 'SD|n\\(\\)')])
simulationAccuracyByTypeDF$`Accuracy SD` <- t(simulationAccuracyByType[str_which(names(simulationAccuracyByType), 'SD')])
nSimsimulationAccuracyByTypeDF <- simulationAccuracyByType[(names(simulationAccuracyByType) == 'n()')]


## more plots - catapillar plot of correlation differences 
plotDat <- allData[!is.na(allData$fisherZDiff),]
plotDat <- plotDat[order(plotDat$correlationDifference.ro),]
names(plotDat)[names(plotDat)=="source"]  <- "Project"

catPlot <- ggplot(plotDat, aes(1:nrow(plotDat), 
                               correlationDifference.ro, 
                               ymin = correlationDifference.ro + ztor(1.96*seDifference.ro), 
                               ymax = correlationDifference.ro -  ztor(1.96*seDifference.ro), 
                               color = Project)) + 
  geom_point(na.rm = T) + geom_errorbar(na.rm = T) + theme_classic() + 
  xlab(NULL) + ylab("Correlation difference") + 
  theme(axis.title.x=element_blank(),                                            
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank()) +  ochRe::scale_colour_ochre(palette = "tasmania") 
                                
## Loading data from mixture model to avoid having to rerun the whole model
HDISimple <- readRDS("Data/mixtureModelOutput/HDIsAlphaSimple.rds")
alphaSimple <- readRDS("Data/mixtureModelOutput/alphaSimple.rds")
phiSimple <- readRDS("Data/mixtureModelOutput/phiSimple.rds")
phiSimpleHDI <- readRDS("Data/mixtureModelOutput/HPDphiSimple.rds")
propBelow.1.BMM <- readRDS("Data/mixtureModelOutput/ValuesBelow.1.rds")
jagData <- read_csv("Data/mixtureModelOutput/jagData.csv")

# Ploting the results
mixtureModelPlot <- ggplot(jagData, aes(x = correlation.o, y = correlation.r,  color = probRealEffect, size = n.r)) +
  geom_abline( slope = 1, intercept = 0) + geom_point(alpha = .8, na.rm = T)+ theme_classic() +
  guides(color = guide_legend(title = "Posterior\nassignment\nrate"),
         shape = guide_legend(title = "True effect size < 0.1"),
         size = guide_legend(title = "Replication\nSample size",
                             values= trans_format("identity", function(x) round(exp(x),0)), order = 2)) +
  scale_size(trans = "log", breaks = c(150, 3000, 60000)) + geom_point(colour = "black", na.rm = T, size = .5, shape = 3) +
  xlab("Original correlation")+ ylab("Replication correlation") + ylim(c(-.5, 1))+ xlim(c(-0, 1)) 