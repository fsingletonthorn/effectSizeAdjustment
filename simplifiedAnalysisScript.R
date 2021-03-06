source(file = "Analysis/Wilson score interval and flat violin.R")
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
library(car)
library(broom)
library(gridExtra)
library(grid)

# Excluding missind data
allData <- allData[!is.na(allData$fis.o) & !is.na(allData$fis.r) & !is.na(allData$n.o) & !is.na(allData$n.r),]

allData$id <- 1:nrow(allData)

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

# Function to create OK looking talbe of RMA output
niceMLMESum <- function(REMod) {
  data_frame(Estimate = c(REMod$b, rep(NA, 4)), "95% CI LB" = c(REMod$ci.lb, rep(NA, 4)), "95% CI UB" = c(REMod$ci.ub, rep(NA, 4)), SE = c(REMod$se,  rep(NA, 4)), p = c(ifelse(REMod$pval<.001, "< .001", REMod$pval),  rep(NA, 4)), 
             "Random effects" = c(NA, paste0("Project variance = ", round(REMod$sigma2[1], 3), ", n = ", 
                                             REMod$s.nlevels[1]),
                                  paste0("Article variance = ", round(REMod$sigma2[2], 3), ", n = ", REMod$s.nlevels[2]), 
                                  paste0("Effect variance = ", round(REMod$sigma2[3], 3), ", n = ", REMod$s.nlevels[3]),
                                  paste0("QE(",REMod$k-1, ") = ", round(REMod$QE, 2),  ", p ", ifelse(REMod$QEp <.001, "< .001", paste("=" , round(REMod$QEp, 2))))))
}

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

# approximate SEs
allData$seFishAprox.r  <- ifelse( is.na(allData$seFish.r ), sqrt(1/(allData$n.r-3)), allData$seFish.r)
allData$seFishAprox.o  <- ifelse( is.na(allData$seFish.o ), sqrt(1/(allData$n.o-3)), allData$seFish.o)

# ## Finding minimium effect that would have been significant in the original study - 
minimumEffectDetectableZ <- qnorm(.05, mean = 0, sd = allData$seFishAprox.o, lower.tail = FALSE)
CIs90UB <-   allData$fis.r + (qnorm(0.95)*allData$seFishAprox.r)
CIs90LB <-   allData$fis.r - (qnorm(0.95)*allData$seFishAprox.r)

upper95.r <- allData$fis.r + ( 1.96 * allData$seFishAprox.r )
lower95.r <- allData$fis.r - ( 1.96 * allData$seFishAprox.r )

# Extracting all studies 
nonMatchesStatisticalSig <- sum( c(lower95.r > 0 & !allData$significantSameDirection.r,lower95.r < 0 & allData$significantSameDirection.r), na.rm = T )

minimumEffectDetectableR <- ztor(minimumEffectDetectableZ)

# calculating the number of equiv studies
proportionEquiv <- mean(CIs90UB < minimumEffectDetectableZ, na.rm =T)
nEquiv <- sum(CIs90UB < minimumEffectDetectableZ, na.rm =T)

allData$statisticallyEquiv.ro <- CIs90UB < minimumEffectDetectableZ

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
            nTrue = sum(!is.na(fisherZDiff)))
# extracting amount of change in subsets
# Sig results replication only (same direction)
reductionSignifcantR <- allData %>% 
  filter(significantSameDirection.r & !is.na(significantSameDirection.r)) %>% 
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T), 
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),    
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum((allData$significantSameDirection.r & !is.na(allData$fisherZDiff))))
# Significant TOST (Exclude those without moderate evidence for the one sided alternative)
reductionEquiv <- allData %>% 
  filter(!statisticallyEquiv.ro) %>%
  summarise(meanPropChange = mean((fisherZDiff)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiff, na.rm = T), sdDiff = sd(fisherZDiff, na.rm = T),   
            medianPropChange = median((fisherZDiff)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiff, na.rm = T),
            nTrue = sum(!allData$statisticallyEquiv.ro, na.rm = T))

# Bringing all of the above together:
tableReductions <- rbind(Overall = reductionAllData, StatisticalSignificance = reductionSignifcantR, Nonequivalence = reductionEquiv)

# Converting into Rs 
tableReductions[,-which(str_detect(names(tableReductions), "Prop|nTrue|nValid"))] <- ztor(tableReductions[,-which(str_detect(names(tableReductions), "Prop|nTrue|nValid"))])
# Renaming 
names(tableReductions) <- c("Mean proportion change", "Mean replication ES", "Mean original ES", "Mean ES difference", 
                            "SD difference", "Median proportion change", "Median replicaiton ES", "Median original ES", "Median ES difference",
                            "n included")

# Reordering cols
tableReductions <- tableReductions[c("n included", "Mean original ES", "Median original ES",  "Mean replication ES", "Median replicaiton ES", 
                                     "Mean ES difference",  "Median ES difference", "SD difference", "Mean proportion change", "Median proportion change")]

kableReductions <- kable(tableReductions, digits = 2)

#########################################################
#### Plots of the effects plotted against each other ####
#########################################################

# Plotting effects on each other 
plotAllData <- ggplot(allData, aes(correlation.o, correlation.r, size = n.r, colour = as.factor(source))) +
  geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = 1) +
  scale_colour_viridis_d() + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + 
  guides( size = guide_legend(title = "Replication\nSample size", values= trans_format("identity", function(x) round(exp(x),0)), order = 2),
          colour = guide_legend(title = "Replication projects"), override.aes = list(alpha = 1), order = 1) +
  xlab("Original correlation")+ ylab("Replication correlation") + scale_size(trans = "log", breaks = c(10, 100, 1000, 10000))

# Plotting only non-equivalent studies  
plotNonequiv <- ggplot(allData[!allData$statisticallyEquiv.ro & !is.na(allData$statisticallyEquiv.ro),], aes(correlation.o, correlation.r, size = n.r, colour = as.factor(source))) +
  geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = 1) +
  scale_colour_viridis_d() + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + 
  guides( size = guide_legend(title = "Replication\nSample size", values= trans_format("identity", function(x) round(exp(x),0)), order = 2),
          colour = guide_legend(title = "Replication projects"), override.aes = list(alpha = 1), order = 1) +
  xlab("Original correlation")+ ylab("Replication correlation") + scale_size(trans = "log", breaks = c(10, 100, 1000, 10000))

# Plotting just significant studies (sig and in same direction)
plotSigR <- ggplot(allData[allData$significantSameDirection.r==TRUE,], aes(correlation.o, correlation.r, size = n.r, colour = as.factor(source))) +
  geom_abline( slope = 1, intercept = 0) +  geom_point(na.rm = T, alpha = 1) +
  scale_colour_viridis_d() + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + 
  guides( size = guide_legend(title = "Replication\nSample size", values= trans_format("identity", function(x) round(exp(x),0)), order = 2),
          colour = guide_legend(title = "Replication projects"), override.aes = list(alpha = 1), order = 1) +
  xlab("Original correlation")+ ylab("Replication correlation") + scale_size(trans = "log", breaks = c(10, 100, 1000, 10000))

### Descriptives
nNotSig.r <- sum(allData$significant.r)

######################################
###### Multilevel meta-analysis ######
######################################

###### Analysis 1
# Random effects model with random effects for authors nested within source
REMod <- rma.mv(yi = fisherZDiff, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id,  data = allData)

# Checking if anything changes excluding those where the SE is approximated
REModValidSE <- rma.mv(yi = fisherZDiff, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id,  data = allData[!is.na(allData$seFish.o) & !is.na(allData$seFish.r),])
validSEDiff <- REMod$b - REModValidSE$b
validSEDiff <- abs(REMod$sigma2 - REModValidSE$sigma2)

# Empirical Bayes estimates for random Effects 
BLUPsSource <- ranef(REMod)[1]

niceREModSum<- niceMLMESum(REMod)

REModSum <- modelOutputSummary(REMod)

## Analysis 2
# The first model but with only significant replications
REModOnlySigR <- rma.mv(yi = fisherZDiff, V = allData[allData$significantSameDirection.r==TRUE,]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id, data = allData[allData$significantSameDirection.r==TRUE,])
# Checking if anything changes excluding those where the SE is approximated
REModOnlySigRValidSE <- rma.mv(yi = fisherZDiff, V = allData[allData$significantSameDirection.r==TRUE,]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id, data = allData[allData$significantSameDirection.r==TRUE & !is.na(allData$seFish.o) & !is.na(allData$seFish.r),])
validSEDiff[2] <- REModOnlySigR$b - REModOnlySigRValidSE$b

## Analysis 3
# The first model but with only non-equiv 
REModNonequiv <- rma.mv(yi = fisherZDiff, V = allData[allData$statisticallyEquiv.ro==FALSE & !is.na(allData$statisticallyEquiv.ro==FALSE),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id, data = allData[allData$statisticallyEquiv.ro==FALSE & !is.na(allData$statisticallyEquiv.ro==FALSE),])
# Checking if anything changes excluding those where the SE is approximated
REModNonequivValidSE <- rma.mv(yi = fisherZDiff, V = allData[allData$statisticallyEquiv.ro==FALSE & !is.na(allData$statisticallyEquiv.ro==FALSE),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id, data = allData[allData$statisticallyEquiv.ro==FALSE & !is.na(allData$statisticallyEquiv.ro==FALSE) & !is.na(allData$seFish.o) & !is.na(allData$seFish.r),])
validSEDiff[3] <- REModNonequiv$b - REModNonequivValidSE$b

# Collecting things for easy reporting
modRes <- data.frame( modelN = REMod$k, modelEstimate = REMod$b, MLM95lb = REMod$ci.lb, MLM95ub = REMod$ci.ub, row.names = "Overall")
modRes1 <- data.frame(modelN = REModOnlySigR$k, modelEstimate = REModOnlySigR$b, MLM95lb = REModOnlySigR$ci.lb, MLM95ub = REModOnlySigR$ci.ub, row.names = "StatisticalSignificance")
modRes2 <- data.frame(modelN = REModNonequiv$k, modelEstimate = REModNonequiv$b, MLM95lb = REModNonequiv$ci.lb, MLM95ub = REModNonequiv$ci.ub, row.names = "Nonequivalence")

modSumaries <- rbind(modRes, modRes1, modRes2)

# Estiating the degree of effect size change as a proportion of the average effect size in psychology 
modSumaries$`Estimated % attenuation` <- (modSumaries$modelEstimate/mean(allData$fis.o, na.rm = T))*100
modSumaries$`LB % attenuation` <- (modSumaries$MLM95lb/mean(allData$fis.o, na.rm = T)*100)
modSumaries$`UB % attenuation` <- (modSumaries$MLM95ub/mean(allData$fis.o, na.rm = T)*100)
modSumariesR <- modSumaries
# converting from r to z
modSumariesR[2:4] <- ztor(modSumaries[2:4])
niceModelSums <- lapply(X = list("All Data" = REMod, "Non-equivalent studies" = REModNonequiv, "Only Signficant replications" = REModOnlySigR), niceMLMESum)

tableAllEstimates <- merge.data.frame(tableReductions, modSumaries, by = "row.names", sort = F)

tableAverageDecrease <- allData %>%
  group_by(source) %>%
  dplyr::summarise(mean=mean(fisherZDiff, na.rm=T), sd=sd(fisherZDiff, na.rm=T))


### Leave one out cross validation for main model, excluding first each source ####
# Commented out to avoid having to run this each time I knit this document
#  LOOTracking <- list()
#  
#  for(i in 1:length(unique(allData$source))) {
#    # exlude <- unique(allData$source)[i]
#    tempData <- allData[-which(allData$source == unique(allData$source)[i]),]
#    LOOTracking[[i]]                                  <-  rma.mv(yi = tempData$fisherZDiff, V = tempData$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id,  data = tempData)
#    LOOTracking[[i+1*length(unique(allData$source))]] <-  rma.mv(yi = fisherZDiff, V = tempData$seDifference.ro^2, mod = cleanedpVal.r, random =  ~ 1|source/authorsTitle.o/id,  data = tempData)
#    LOOTracking[[i+ 2*length(unique(allData$source))]] <-  rma.mv(yi = fisherZDiff, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id,  data = (tempData[tempData$statisticallyEquiv.ro==FALSE,]))
#    LOOTracking[[i+ 3*length(unique(allData$source))]] <-  rma.mv(yi = fisherZDiff, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id,  data = tempData[tempData$significantSameDirection.r==TRUE,])
#    
#  }
#  
#  # extrating estiamtes from list
#   LOOProject <- sapply(LOOTracking, function(m) {m[c(1,5, 108)]}, simplify = T)
#   LOOProjectDF <- data.frame(as.data.frame(t(LOOProject)))
#   LOOProjectDF[,1] <- str_remove_all(LOOProjectDF[,1], ",.+|c\\(")
#   LOOProjectDF[,2] <- str_remove_all(LOOProjectDF[,2], ",.+|c\\(")
#   LOOProjectDF <- data.frame(b = as.numeric(LOOProjectDF$b), p = as.numeric(LOOProjectDF$pval), call = as.character(LOOProjectDF$call))
 
# write.csv( LOOProjectDF, "Data/LOOProjectSimp.csv", row.names = F); beepr::beep(8)
 
# This reads in the data from the LOO analysis above 
LOOProjectDF <- read_csv("Data/LOOProjectSimp.csv")

# Summarising (excluding models with all data)
LOOProjectDFSum <- LOOProjectDF %>%
  group_by(Subsample = call) %>%
  dplyr::summarise('Proportion significant'= mean(p<.05), 'Minimum estimate' = quantile(b)[1],  '25th percentile' = quantile(b)[2], 'Median' = quantile(b)[3], '75th percentile' = quantile(b)[4], 'Maximum estimate' = quantile(b)[5])

# Renaming files
namesR <- as.character(LOOProjectDFSum$Subsample)
namesR[str_detect(LOOProjectDFSum$Subsample, "mods \\=")] <- "P value as Moderator"
namesR[str_detect(LOOProjectDFSum$Subsample, "signif")] <- "Only significant replications"
namesR[str_detect(LOOProjectDFSum$Subsample, "Equiv")] <- "Only non-equivalent replications"
namesR[str_detect(LOOProjectDFSum$Subsample, "^((?!mods|statistic|signif).)*$")] <- "All data"

LOOProjectDFSum$Subsample <- namesR

# Again commented out to avoid knitting this every time
#LOO removing studies - this takes ~ 20 mins to run, it was run once and then saved ####
#   for(i in 1:length(unique(allData$authorsTitle.o))) {
#     tempData <- allData[-which(allData$authorsTitle.o == unique(allData$authorsTitle.o)[i]),]
#     LOOTracking[[i]] <- rma.mv(yi = fisherZDiff, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id,  data = tempData)
#     LOOTracking[[i+ (length(unique(allData$authorsTitle.o)))]] <-  rma.mv(yi = fisherZDiff, V = tempData$seDifference.ro^2, mod = cleanedpVal.r, random =  ~ 1|source/authorsTitle.o/id,  data = tempData)
#     LOOTracking[[i+ 2*(length(unique(allData$authorsTitle.o)))]] <-  rma.mv(yi = fisherZDiff, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id,  data = tempData[tempData$statisticallyEquiv.ro==FALSE,])
#     LOOTracking[[i+ 3*(length(unique(allData$authorsTitle.o)))]] <-  rma.mv(yi = fisherZDiff, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id,  data = tempData[tempData$significantSameDirection.r==TRUE,])
#   }
#
# 
# LOOStudy <- sapply(LOOTracking, function(m) {m[c(1,5, 108)]}, simplify = T)
# LOOStudyDF <- data.frame(as.data.frame(t(LOOStudy)))
# LOOStudyDF[,1] <- str_remove_all(LOOStudyDF[,1], ",.+|c\\(")
# LOOStudyDF[,2] <- str_remove_all(LOOStudyDF[,2], ",.+|c\\(")
# LOOStudyDF <- data.frame(b = as.numeric(LOOStudyDF$b), p = as.numeric(LOOStudyDF$pval), call = as.character(LOOStudyDF$call))
##  
#   write.csv(LOOStudyDF ,"Data/LOOStudyDFSimp.csv"); beepr::beep(8)
# This reads in the data produced by the above  
LOOStudyDF <- read.csv("Data/LOOStudyDFSimp.csv")
LOOStudyDF <-LOOStudyDF[,c(2,3,4)]

LOOStudyDFSum <- LOOStudyDF %>%
  group_by(Subsample = call) %>%
  dplyr::summarise('Proportion significant'= mean(p<.05), 'Minimum estimate' = quantile(b)[1],  '25th percentile' = quantile(b)[2], 'Median' = quantile(b)[3], '75th percentile' = quantile(b)[4], 'Maximum estimate' = quantile(b)[5])

# This is just relabeling things for easy presentation 
namesR <- as.character(LOOStudyDFSum$Subsample)
namesR[str_detect(LOOStudyDFSum$Subsample, "mods \\=")] <- "P value as Moderator"
namesR[str_detect(LOOStudyDFSum$Subsample, "signif")] <- "Only significant replications"
namesR[str_detect(LOOStudyDFSum$Subsample, "Equiv")] <- "Only non-equivalent replications"
namesR[str_detect(LOOStudyDFSum$Subsample, "^((?!mods|statistic|signif).)*$")] <- "All data"
LOOStudyDFSum$Subsample <- namesR

# Loo on the effect level
# Again commented out to avoid knitting this every time
# LOO removing studies - this takes ~ 20 mins to run, it was run once and then saved ####
#  LOOTracking <- list()
#   for(i in 1:nrow(allData)) {
#        tempData <- allData[-i,]
#        LOOTracking[[i]] <- rma.mv(yi = tempData$fisherZDiff, V = tempData$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id,  data = tempData)
#        LOOTracking[[i+ nrow(allData)]] <-  rma.mv(yi = tempData$fisherZDiff, V = tempData$seDifference.ro^2, mod = cleanedpVal.r, random =  ~ 1|source/authorsTitle.o/id,  data = tempData)
#        LOOTracking[[i+ 2*nrow(allData)]] <-  rma.mv(yi = fisherZDiff, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id,  data = tempData[tempData$statisticallyEquiv.ro==FALSE,])
#        LOOTracking[[i+ 3*nrow(allData)]] <-  rma.mv(yi = fisherZDiff, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o/id,  data = tempData[tempData$significantSameDirection.r==TRUE,])
#   }
#    
#    LOOEffect <- sapply(LOOTracking, function(m) {m[c(1,5, 108)]}, simplify = T)
#    LOOEffectDF <- data.frame(as.data.frame(t(LOOEffect)))
#    LOOEffectDF[,1] <- str_remove_all(LOOEffectDF[,1], ",.+|c\\(")
#    LOOEffectDF[,2] <- str_remove_all(LOOEffectDF[,2], ",.+|c\\(")
#    LOOEffectDF <- data.frame(b = as.numeric(LOOEffectDF$b), p = as.numeric(LOOEffectDF$pval), call = as.character(LOOEffectDF$call))
# 
#  write.csv(LOOEffectDF ,"Data/LOOSEffectSimp.csv"); beep(8)
# This reads in the data produced by the above  
LOOEffectDF <- read.csv("Data/LOOSEffectSimp.csv")
LOOEffectDF <- LOOEffectDF[,c(2,3,4)]

LOOEffectDFSum <- LOOEffectDF %>%
  group_by(Subsample = call) %>%
  dplyr::summarise('Proportion significant'= mean(p<.05), 'Minimum estimate' = quantile(b)[1],  '25th percentile' = quantile(b)[2], 'Median' = quantile(b)[3], '75th percentile' = quantile(b)[4], 'Maximum estimate' = quantile(b)[5])

# This is just relabeling things for easy presentation 
namesR <- as.character(LOOEffectDFSum$Subsample)
namesR[str_detect(LOOEffectDFSum$Subsample, "mods \\=")] <- "P value as Moderator"
namesR[str_detect(LOOEffectDFSum$Subsample, "signif")] <- "Only significant replications"
namesR[str_detect(LOOEffectDFSum$Subsample, "Equiv")] <- "Only non-equivalent replications"
namesR[str_detect(LOOEffectDFSum$Subsample, "^((?!mods|statistic|signif).)*$")] <- "All data"
LOOEffectDFSum$Subsample <- namesR

LOOoutput <- bind_rows(LOOProjectDFSum,LOOStudyDFSum,LOOEffectDFSum, .id = "LOO exclusions")
LOOoutput$`LOO exclusions` <- c(rep("Replication project", 4),rep( "Study",4), rep( "Effect",4))

LooMaxDiff <- max(c(abs(REMod$b - LOOoutput$`Minimum estimate`[LOOoutput$Subsample == "All data"]), 
  abs(REMod$b - LOOoutput$`Maximum estimate`[LOOoutput$Subsample == "All data"]),
  # abs(REMod.p.val.cleaned$b[1] - LOOoutput$`Minimum estimate`[LOOoutput$Subsample == "P value as Moderator"]),
  # abs(REMod.p.val.cleaned$b[1] - LOOoutput$`Maximum estimate`[LOOoutput$Subsample == "P value as Moderator"]),
  abs(REModOnlySigR$b - LOOoutput$`Minimum estimate`[LOOoutput$Subsample == "Only significant replications"]),
  abs(REModOnlySigR$b - LOOoutput$`Maximum estimate`[LOOoutput$Subsample == "Only significant replications"]),
  abs(REModNonequiv$b - LOOoutput$`Minimum estimate`[LOOoutput$Subsample == "Only non-equivalent replications"]),
  abs(REModNonequiv$b - LOOoutput$`Maximum estimate`[LOOoutput$Subsample == "Only non-equivalent replications"])
))


# Adding something to adjust for the hetrogeneity of multistudy papers -
adjustedDat <- allData
adjustedDat$adjusted_n.r <- adjustedDat$n.r/adjustedDat$n_studies
adjustedDat$adjustedSE.ro <-  sqrt(1/(adjustedDat$n.o-3) + 1/(adjustedDat$adjusted_n.r-3))
adjustedDat$id2 <- 1:nrow(adjustedDat)

# Analysis 1
REModAdjusted <- rma.mv(yi = fisherZDiff, V = adjustedSE.ro^2, random =  ~ 1|source/authorsTitle.o/id2,  data = adjustedDat)
# Analysis 2
REModOnlySigAdjusted <- rma.mv(yi = fisherZDiff, V = adjustedSE.ro^2, random =  ~ 1|source/authorsTitle.o/id2,  data = adjustedDat[allData$significantSameDirection.r==TRUE,])
# Analysis 3
REModNonequivAdjusted <- rma.mv(yi = fisherZDiff, V = adjustedSE.ro^2, random =  ~ 1|source/authorsTitle.o/id2,  data = adjustedDat[adjustedDat$statisticallyEquiv.ro==FALSE & !is.na(adjustedDat$statisticallyEquiv.ro==FALSE),])

#####################
###### Figures ######
#####################

# Fishers Z change 
fisherZChangePlot <- ggplot(allData,aes(y = fisherZDiff,x= source, 
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
        axis.title.y = element_blank()) + 
  ylab("Differences in Fisher's Z transformed effect sizes")

ggsave(plot =  fisherZChangePlot, 
       file = "Figures/ViolinPlotZscoreChange.jpg", 
       width = 6, height = 8, device = "jpeg", dpi = 320)


## more plots - catapillar plot of correlation differences 
plotDat <- allData[!is.na(allData$fisherZDiff),]
plotDat <- plotDat[order(plotDat$fisherZDiff),]
names(plotDat)[names(plotDat)=="source"]  <- "Project"

catPlot <- ggplot(plotDat, aes(1:nrow(plotDat), 
                               ztor(fisherZDiff), 
                               ymin = ztor( fisherZDiff + 1.96*seDifference.ro), 
                               ymax =  ztor( fisherZDiff - 1.96*seDifference.ro), 
                               color = Project)) + 
  geom_point(na.rm = T) + geom_errorbar(na.rm = T) + theme_classic() + 
  xlab(NULL) + ylab("Correlation difference") + 
  theme(axis.title.x=element_blank(),                                            
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank()) + scale_colour_viridis_d()

## Loading data from mixture model to avoid having to rerun the whole model
HDISimple <- readRDS("Data/mixtureModelOutput/HDIsAlphaSimple.rds")
alphaSimple <- readRDS("Data/mixtureModelOutput/alphaSimple.rds")
phiSimple <- readRDS("Data/mixtureModelOutput/phiSimple.rds")
phiSimpleHDI <- readRDS("Data/mixtureModelOutput/HPDphiSimple.rds")
propBelow.1.BMM <- readRDS("Data/mixtureModelOutput/ValuesBelow.1.rds")
jagData <- read_csv("Data/mixtureModelOutput/jagData.csv")
rhatAlpha <- readRDS("Data/mixtureModelOutput/rhatAlpha.rds")
rhatPhi <- readRDS("Data/mixtureModelOutput/rhatPhi.rds")

# Ploting the results
mixtureModelPlot <- ggplot(jagData, aes(x = correlation.o, y = correlation.r,  color = probRealEffect, size = n.r)) + # scale_colour_continuous(low = "#ece2f0", high = "#1c9099")+
  geom_abline( slope = 1, intercept = 0) + geom_point(alpha = 1, na.rm = T)+ scale_fill_brewer()+theme_classic() +
  guides(color = guide_legend(title = "Posterior\nassignment\nrate"),
         shape = guide_legend(title = "True effect size < 0.1"),
         size = guide_legend(title = "Replication\nSample size",
                             values= trans_format("identity", function(x) round(exp(x),0)), order = 2)) +
  # scale_color_continuous(breaks = c(.5,.75, .999)) +
  scale_size(trans = "log", breaks = c(10, 100, 1000, 10000)) + geom_point(colour = "black", na.rm = T, size = .1, shape = 3) +
  xlab("Original correlation")+ ylab("Replication correlation") + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) 


# Summary table of included data


incStud <- allData %>% 
  group_by(source) %>%
  dplyr::summarise(sum(!is.na(fisherZDiff)))

studies <- c("Camerer, C. F., Dreber, A., Forsell, E., Ho, T.-H., Huber, J., Johannesson, M., . . . Wu, H. (2016). Evaluating replicability of laboratory experiments in economics. Science, 351(6280), 1433. DOI: 10.1126/science.aaf0918",
             "Camerer, C. F., Dreber, A., Holzmeister, F., Ho, T.-H., Huber, J., Johannesson, M., . . . Wu, H. (2018). Evaluating the replicability of social science experiments in Nature and Science between 2010 and 2015. Nature Human Behaviour, 2(9), 637-644. doi:10.1038/s41562-018-0399-z",
             "Cova, F., Strickland, B., Abatista, A., Allard, A., Andow, J., Attie, M., . . . Colombo, M. (2018). Estimating the reproducibility of experimental philosophy. Review of Philosophy and Psychology, 1-36. doi: 10.1007/s13164-018-0407-2.",
             "Ebersole, C. R., Atherton, O. E., Belanger, A. L., Skulborstad, H. M., Allen, J. M., Banks, J. B., . . . Nosek, B. A. (2016). Many Labs 3: Evaluating participant pool quality across the academic semester via replication. Journal of Experimental Social Psychology, 67, 68-82. doi:10.1016/j.jesp.2015.10.012",
             "Klein, R. A., Ratliff, K. A., Vianello, M., Adams, R. B., Bahník, Š., Bernstein, M. J., . . . Nosek, B. A. (2014). Investigating Variation in Replicability. Social Psychology, 45(3), 142-152. doi:10.1027/1864-9335/a000178 $_b$",
             "Klein, R. A., Vianello, M., Hasselman, F., Adams, B. G., Adams, R. B., Alper, S., ... Nosek, B. A. (2018). Many Labs 2: Investigating Variation in Replicability Across Samples and Settings. Advances In Methods and Practices in Psychological Science, 1(4), 443-490. doi:10.1177/2515245918810225",
             "Open Science Collaboration. (2015). Estimating the reproducibility of psychological science. Science, 349(6251), aac4716. doi:10.1126/science.aac4716",
             "Soto et al (2018) $_a$")

nPossibleStudies <- sum(c(18, 21, 37, 9, 16, 28, 97, 121) )

replicationProjects <- data_frame("Replication projects" = studies,  "Number of replication studies performed" = c(18, 21, 37, 9, "16 (13 effects)", 28, 97, 121), "Reported replication rate (statistically significant results in the same direction)" = c("61%","62%","78%","33%","88% (85%)", "54%","36%","86%")
                                  , "Included studies" = as.numeric(unlist(incStud[,2])))
weights <- c(as.numeric(replicationProjects$`Number of replication studies performed`))
weights[5] <- 16
weights2 <- c(18, 21, 37, 9, 13, 28, 97, 121)

replicationProjects[9,] <- c("All projects", sum(c(as.numeric(replicationProjects$`Number of replication studies performed`),16), na.rm = T), 
                             paste0(round(weighted.mean(as.numeric(str_split(replicationProjects$`Reported replication rate (statistically significant results in the same direction)`, "%", simplify = T)[,1], na.rm = T), weights))
                                    #,"% (",
                                    # round(weighted.mean(as.numeric(str_split(replicationProjects$`Reported replication rate (statistically significant results in the same direction)`, "%", simplify = T)[,1], na.rm = T), weights2)), "%)"),
                             , "%"), sum(replicationProjects$`Included studies`))

replicationProjects$`Percent of performed replication studies included` <-paste0( round( (as.numeric(str_split(replicationProjects$`Included studies`, " ", simplify = T)) / as.numeric(str_split(replicationProjects$`Number of replication studies performed`, " ", simplify = T)[,1]))*100), "%")

replicationProjects$`Percent of performed replication studies included`[5] <- paste0(replicationProjects$`Percent of performed replication studies included`[5], " (", round((as.numeric(replicationProjects$`Included studies`[5])-3) /13 *100), "%)")

#  I^2 for ML metas =  following Nakagawa, 2012 #1023}
W <- diag(1/allData$seDifference.ro^2)
X <- model.matrix(REMod)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
Io2 <- 100 * sum(REMod$sigma2) / (sum(REMod$sigma2) + (REMod$k-REMod$p)/sum(diag(P)))

# Estimates of the variability 
confIntsREMod <- confint(REMod)
colourPal <- c("#440154FF", "#46337EFF", "#365C8DFF", "#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
 # c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00',"#8dd3c7",'#a65628','#f781bf')

#### MLM variability plots ####
# # MLMplot 
projectEst <-  ranef(REMod)[[1]] +  REMod$b[[1]]
projectEst$projectName <- row.names(projectEst) %>%
  str_remove( "(\\s\\n.*)|(\\n.*)")
projectEst$names <- projectEst$projectName 
projectEst$colour <- colourPal

# Calculating Estimated study effects
studyEst <-  ranef(REMod)[[2]]
studyEst$names <- row.names(studyEst)
studyEst$projectName <- row.names(studyEst) %>%
  str_remove( "(\\s\\n.*)|(\\n.*)")
temporaryOut <- left_join(studyEst, projectEst, by = c("projectName"))
studyEst$intrcpt <- studyEst$intrcpt + temporaryOut$intrcpt.y
studyEst$colour <- temporaryOut$colour

# Calcualting estiamted effect - effects
effectEst <-  ranef(REMod)[[3]]
effectEst$names <- row.names(effectEst) %>%
  str_remove( "/\\d*$")
effectEst$projectName <- row.names(effectEst) %>%
  str_remove( "(\\s\\n.*)|(\\n.*)")
effectEst$names <- row.names(effectEst) %>%
  str_remove( "/\\d*$")
temporaryOut <- left_join(effectEst, studyEst, by = c("names" = "names"))
effectEst$intrcpt <- effectEst$intrcpt  + temporaryOut$intrcpt.y
effectEst$colour <- temporaryOut$colour

# Allests
allEsts <- bind_rows(projectEst, studyEst, effectEst, .id = "level")
allEsts$level <- factor(allEsts$level, labels = c("Project", "Study", "Effect"))
allEsts$project <- factor(str_remove( allEsts$names, "(\\s\\n.*)|(\\n.*)") %>% str_remove(",")) 

project_title <- grobTree(rectGrob(gp=gpar(fill="lightgrey")),
                 textGrob("Project", x=0.5, hjust=0.5,
                          gp=gpar(col="black", cex=1.3)))

article_title <- grobTree(rectGrob(gp=gpar(fill="lightgrey")),
                 textGrob("Article", x=0.5, hjust=0.5,
                          gp=gpar(col="black", cex=1.3)))

effect_title <- grobTree(rectGrob(gp=gpar(fill="lightgrey")),
                 textGrob("Effect", x=0.5, hjust=0.5,
                          gp=gpar(col="black", cex=1.3)))


######## Plots of BLUPs at each level ##### 

pannel1 <-  allEsts %>%
  filter(level == "Project") %>%
  ggplot(aes(x = intrcpt,  fill = projectName)) +
  geom_histogram( bins = 40) +
#  geom_jitter(aes(x = intrcpt, y = 6), height = .8, width = 0, size = 2.5, alpha = .7, colour = "lightgrey") +
  scale_fill_manual(values = colourPal)  +
  scale_colour_manual(values = colourPal) +
  ylab("") + xlab("") + scale_x_continuous(limits = c(-.7,.7), breaks = round(seq(from = -.6,to =  .6, by = .2),2 )) +
  theme_classic() +
  theme(#axis.title.y=element_blank(),
    #axis.text.y=element_blank(),
    #axis.line.y = element_blank(),
    #axis.ticks.y=element_blank(),
    legend.position = "none")

pannel2 <-  allEsts %>%
  filter(level == "Study") %>% 
  ggplot(aes(x = intrcpt,  fill = project) ) + # , colour = projectName, fill = projectName)) +
  geom_histogram( bins = 40)  + ylim(0,40) + 
#  geom_jitter(aes(x = intrcpt, y = 6), height = .8, width = 0, size = 2.5, alpha = .7, colour = "lightgrey") +
  scale_fill_manual(values = colourPal)  +
  scale_colour_manual(values = colourPal) +
  ylab("") + xlab("") +
  scale_x_continuous(limits = c(-.7,.7),
                     breaks = round(seq(from = -.6,to =  .6, by = .2),2 )) +
  theme_classic() +
  theme(#axis.title.y=element_blank(),
    #axis.text.y=element_blank(),
    #axis.line.y = element_blank(),
    #axis.ticks.y=element_blank(),
    legend.position = "none")


pannel3 <-  allEsts %>%
  filter(level == "Effect") %>%
  ggplot(aes(x = intrcpt,  fill = project) ) + # , colour = projectName, fill = projectName)) +
  geom_histogram( bins = 40) + ylim(0,40) + 
 # geom_jitter(aes(x = intrcpt, y = 6), height = .8, width = 0, size = 2.5, alpha = .5, colour = "lightgrey") +
  scale_fill_manual(values = colourPal)  +
  scale_colour_manual(values = colourPal) +
  ylab("") + xlab("")  + 
  scale_x_continuous(limits = c(-.7,.7), breaks = round(seq(from = -.6,to =  .6, by = .2),2 )) +
  theme_classic() +
  theme(#axis.title.y=element_blank(),
    #axis.text.y=element_blank(),
    #axis.line.y = element_blank(),
    #axis.ticks.y=element_blank(),
    legend.position = "bottom",
    legend.title = element_blank())

