# source(file = 'Data/data_collection_cleaning.R')
source(file = "R code effect size estimation/appendixCodeFunctionsJeffreys.R")
source(file = "Analysis/Wilson score interval.R")
library(readr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

# install.packages("BayesFactor")
# library(BayesFactor)


# first off simple calculation of the proportion of studies conducted per published result given a significnat main effect

# How to get to .9 or .75 of the literature being sig, assuming all studies are statistically signifincat

studiesPerPublishedPaper <- paste(round(.75/.44,2), "to", round(.9/.44,2))


# Sample characteristics
# Calculating the number included in the meta-analysis (also the number included in the questions asking about change)
nMeta <- sum((!is.na(allData$fisherZDiff)))
nMetaIndividualPapers <- length(unique(allData$authorsTitle.o[!is.na(allData$fisherZDiff)]))

nEffectSizeFell <- sum(allData$correlationDifference.ro < 0, na.rm = T)
propEffectSizeFell <- mean(allData$correlationDifference.ro < 0, na.rm = T)

CIonProp <- WilsonBinCI(sum(!is.na(allData$correlationDifference.ro < 0)), p = propEffectSizeFell)

# N with valid SEs
nSEsValid <- sum(!is.na(allData$seFish.o) & !is.na(allData$seFish.r))


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


z <- 1:nrow(allData)
x <- data.frame(allData$correlation.o, allData$n.o, allData$correlation.r, allData$n.r)
for(i in 1:nrow(allData)) {
y <- i
z[i] <- bfapplyRep0(x[y,])
}
z <- x[y,]

# Two sided default Bayes factor
allData$bf10 <- apply(data.frame(allData$correlation.r, allData$n.r), MARGIN = 1, FUN = bfapply10) 
allData$bf01 <- 1/allData$bf10 

# one sided ~ probably makes more sense as they have all been shifted to positive direction
allData$bfplus0 <- apply(data.frame(allData$correlation.r, allData$n.r), MARGIN = 1, FUN = bfapplyplus0) 
allData$bf0plus <- 1/allData$bfplus0 

# Replication Bayes factor
allData$bfRep0 <- apply(data.frame(allData$correlation.o, allData$n.o, allData$correlation.r, allData$n.r), MARGIN = 1, FUN = bfapplyRep0) 
allData$bf0Rep <- 1/allData$bfRep0 


sum(allData$bf0Rep < 1, na.rm = T)/sum(!is.na(allData$bf0Rep < 3))


#### Amount of change in subsets ####

# extracting amount of change in subsetss 
# Sig results replication only (same direction)
reductionSignifcantR <- allData %>% 
  filter(significantSameDirection.r) %>% 
  summarise(arr = mean((fis.r - fis.o)/fis.o, na.rm = TRUE), nSig = sum(!is.na(pVal.r)), n = length(fis.o))
# BF01 < 3 (exclude those with moderate or greater evidence for the null)
reductionSbf01LessThan3<- allData %>% 
  filter(bf01<3) %>% 
  summarise(arr = mean((fis.r - fis.o)/fis.o, na.rm = TRUE), nApplicable = sum(!is.na(bf01)), n = length(fis.o))
# BF10 > 3 (Only those with evidence for the alternative)
reductionSbf10MoreThan3 <- allData %>% 
  filter(bf10>3) %>% 
  summarise(arr = mean((fis.r - fis.o)/fis.o, na.rm = TRUE), nApplicable = sum(!is.na(bf10)), n = length(fis.o))
# BF0Plus < 3 (exclude those with moderate evidence for the null, one sided alternative)
reductionSbf0PlusLessThan3 <- allData %>% 
  filter(bf0plus<3) %>% 
  summarise(arr = mean((fis.r - fis.o)/fis.o, na.rm = TRUE), nApplicable = sum(!is.na(bf0plus)), n = length(fis.o))
# BFPlus0 (Exclude those without moderate evidence for the one sided alternative)
reductionSbfPlus0MoreThan3 <- allData %>% 
  filter(bfplus0>3) %>%
  summarise(arr = mean((fis.r - fis.o)/fis.o, na.rm = TRUE), nApplicable = sum(!is.na(bfplus0)), n = length(fis.o))

#### Plots of the effects plotted against each other ####

# Plotting just significant studies (sig and in same direction)
plotSigR <- ggplot(allData[allData$significantSameDirection.r==TRUE,], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) +  geom_point(na.rm = T)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Significant replications in the same direction only")

# two sided prior, Exluding studies with BF10 lower than 3, i.e., without evidence moderate or greater for the alternative
plotBF10Greater3 <- ggplot(allData[(allData$bf10>3),], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) +  geom_point(na.rm = T)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Exluding studies with BF10 lower than 3")

# two sided prior, excluding those with evidence for null > 3
plotBF01Lesser3 <- ggplot(allData[(allData$bf01<3) & (sign(allData$correlation.o) == sign(allData$correlation.r)),], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) +  geom_point(na.rm = T)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Exluding studies with BF01 greater than 3")

# One sided default bayes factor, excluding studies with ev for null > 3
plotBF0plusLesser3 <- ggplot(allData[allData$bf0plus<3,], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) +  geom_point(na.rm = T)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Exluding studies with BF0+ (one sided) greater than 3")

# One sided default bayes factor, including only studies with ev for alternative > 3
plotBFPlus0Greater3 <- ggplot(allData[allData$bfplus0>3,], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) +  geom_point(na.rm = T)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Exluding studies with BF+0 less than than 3")

# Replication bayes factor, including only studies *without* evidence for the null vs the original effect over 3
plotBFRep0Lesser3 <- ggplot(allData[allData$bf0Rep <3,], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) +  geom_point(na.rm = T)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Exluding studies with BF0Rep greater than than 3")

# average effect size decrease, exlcuding studies with BF0Reps greater than 3 



### Descriptives
nNotSig.r <- sum(allData$significant.r)

sum(allData$bf0plus<1, na.rm=T)/sum(!is.na(allData$bf0plus))

#### Meta-analyses ####
# Random effects model with random effects for authors nested within source
REMod <- rma.mv(yi = allData$fisherZDiff, V = allData$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o,  data = allData)
summary(REMod)

# Empirical Bayes estimates for random Effects 
ranef(REMod)

# Converting results into an understandable format:
rDecrease <- ztor(as.numeric(REMod[1]))
rDecreaseCILB <- ztor(as.numeric(REMod[6]))
rDecreaseCIUB <- ztor(as.numeric(REMod[7]))
cDcreaseCIs <- c(rDecrease, rDecreaseCILB, rDecreaseCIUB)
# And in cohen's d: 
dDecrease <- (2*cDcreaseCIs) / sqrt(1-cDcreaseCIs^2)


# The first model but with only significant replications
REModOnlySigR <- rma.mv(yi = fisherZDiff, V = allData[allData$significant.r==TRUE,]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = allData[allData$significant.r==TRUE,])
summary(REModOnlySigR)

correlationDecreaseOnlySig <- ztor(as.numeric(REModOnlySigR[1]))
dDecreaseOnlySig <- (2*correlationDecreaseOnlySig) / sqrt(1-correlationDecreaseOnlySig^2)

# the first model removing all studies with BFs0Plus > 3 (i.e., moderate evidence for a null)
REModBF0PlusGreaterThan3Excluded <- rma.mv(yi = fisherZDiff, V = allData[allData$bf0plus > 3 & !is.na(allData$bf0plus > 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = allData[allData$bf0plus > 3 & !is.na(allData$bf0plus > 3),])
summary(REModBF0PlusGreaterThan3Excluded)

# the first model removing all studies with BFsPlus0 < 3 (i.e., those without evidence for the alternative)
REModBFPlus0LessThan3Excluded <- rma.mv(yi = fisherZDiff, V = allData[allData$bfplus0 < 3 & !is.na(allData$bfplus0 < 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = allData[allData$bfplus0 < 3 & !is.na(allData$bfplus0 < 3),])
summary(REModBFPlus0LessThan3Excluded)

# the first model removing all studies with BFs0Rep < 3 (i.e., those without evidence for the for the null vs. origianl greater than 3)
REModBF0RepGreaterThan3Excluded <- rma.mv(yi = fisherZDiff, V = allData[allData$bfplus0 < 3 & !is.na(allData$bfplus0 < 3),]$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = allData[allData$bfplus0 < 3 & !is.na(allData$bfplus0 < 3),])
summary(REModBFPlus0LessThan3Excluded)


#### Leave one out cross validation for main model, excluding first each source
LOOTracking <- list()

for(i in 1:length(unique(allData$source))) {
  # exlude <- unique(allData$source)[i]
  tempData <- allData[-which(allData$source == unique(allData$source)[i]),]
  LOOTracking[[i]] <- rma.mv(yi = tempData$fisherZDiff, V = tempData$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o,  data = tempData)
}


estimatesLOOSource <- sapply(LOOTracking, function(m) m[1], simplify = T)
pvaluesLOOSource <- sapply(LOOTracking, function(m) m[5], simplify = T)

# Checking that none go below sig 
# pvaluesLOOSource < .05

maximumSourceDif <-  (as.numeric(estimatesLOOSource)- as.numeric(REMod[1]))
 unique(allData$source)[which(as.numeric(REMod[1]) - as.numeric(estimatesLOOSource) ==(max(as.numeric(REMod[1]) - as.numeric(estimatesLOOSource))))]


## LOO removing studies
for(i in 1:length(unique(allData$authorsTitle.o))) {
  # exlude <- unique(allData$source)[i]
  tempData <- allData[-which(allData$authorsTitle.o == unique(allData$authorsTitle.o)[i]),]
  LOOTracking[[i]] <- rma.mv(yi = tempData$fisherZDiff, V = tempData$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o,  data = tempData)
}

estimatesLOOStudy <- sapply(LOOTracking, function(m) m[1], simplify = T)
pvaluesLOOStudy <- sapply(LOOTracking, function(m) m[5], simplify = T)

estimatesLOOStudy 

# No changes in significance
# sum(as.numeric(pvaluesLOOStudy) > .05)
# And no large changes in intercept 
maximumStudyDif <- max(as.numeric(REMod[1]) - as.numeric(estimatesLOOStudy))


LOOTracking <- list()




for(i in 1:length(unique(allData$source))) {
  # exlude <- unique(allData$source)[i]
  tempData <- allData[-which(allData$source == unique(allData$source)[i]),]
  LOOTracking[[i]] <- rma.mv(yi = tempData$fisherZDiff, V = tempData$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o,  data = tempData)
}





View(LOOTracking)
###### Figures ######

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
        axis.title.y = element_blank()) + ylab("Percentage change in effect sizes from original to replicaiton") +
  ylab()
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




