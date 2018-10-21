# source(file = 'Data/data_collection_cleaning.R')
source(file = "R code effect size estimation/appendixCodeFunctionsJeffreys.R")



# first off simple calculation of the proportion of studies conducted per published result given a significnat main effect

# How to get to .9 or .75 of the literature being sig, assuming all studies are statistically signifincat

studiesPerPublishedPaper <- paste(round(.75/.44,2), "to", round(.9/.44,2))


# Sample characteristics
# Calculating the number included in the meta-analysis (also the number included in the study)
nMeta <- sum((!is.na(allData$fisherZDiff)))

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

# Two sided default Bayes factor
allData$bf10 <- apply(data.frame(allData$correlation.r, allData$n.r), MARGIN = 1, FUN = bfapply10) 
allData$bf01 <- 1/bf10 

# one sided ~ probably makes more sense as they have all been shifted to positive direction
allData$bfplus0 <- apply(data.frame(allData$correlation.r, allData$n.r), MARGIN = 1, FUN = bfapplyplus0) 
allData$bf0plus <- 1/bfplus0 

#### Amount of change in subsets ####

# extracting amount of change in subsetss 
# Sig results replication only (same direction)
reductionSignifcantR<- allData %>% 
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

# Plotting just significant studies
plotSigR <- ggplot(allData[allData$significant.r==TRUE,], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) +  geom_point(na.rm = T)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Significant replications in the same direction only")

# two sided prior, Exluding studies with BF10 lower than 3, i.e., without evidence moderate or greater for the alternative
plotBF10Greater3 <- ggplot(allData[(bf10>3),], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) +  geom_point(na.rm = T)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Exluding studies with BF10 lower than 3")

# two sided prior, excluding those with evidence for null > 3
plotBF01Lesser3 <- ggplot(allData[(bf01<3) & (sign(allData$correlation.o) == sign(allData$correlation.r)),], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) +  geom_point(na.rm = T)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Exluding studies with BF01 greater than 3")

# One sided default bayes factor, excluding studies with ev for null > 3
plotBF0plusLesser3 <- ggplot(allData[bf0plus<3,], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) +  geom_point(na.rm = T)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Exluding studies with BF0+ (one sided) greater than 3")

# One sided default bayes factor, including only studies with ev for alternative > 3
plotBFPlus0Greater3 <- ggplot(allData[bfplus0>3,], aes(correlation.o, correlation.r,size = log(n.r), colour = as.factor(source))) +  geom_point(na.rm = T)+  ochRe::scale_colour_ochre(palette = "tasmania") + theme_classic() + ylim(c(-.5, 1))+ xlim(c(-.5, 1)) + scale_shape_manual(values = c(0,1,2,15,16,17)) + ggtitle("Exluding studies with BF+0 less than than 3")

#### Meta-analyses ####
# Random effects model with random effects for authors nested within source
REMod <- rma.mv(yi = fisherZDiff, V = allData$seDifference.ro^2, random =  ~ 1|source/authorsTitle.o,  data = allData)
summary(REMod)
summary(oldREMod)
# oldREMod <- REMod
#

# Random effects model with random effects for article of effect and fixed effect for Source
REModFixSource <- rma.mv(yi = fisherZDiff, V = allData$seDifference.ro^2, random =  ~ authorsTitle.o, mods = ~source,  data = allData)
summary(REModFixSource) 

# The first model with only significant replications
REModOnlySigR <- rma.mv(yi = fisherZDiff, V = allData[allData$significant.r==TRUE,]$seDifference.ro^2, random =  ~ 1|authorsTitle.o/source,  data = allData[allData$significant.r==TRUE,])
summary(REModOnlySigR)

REAuthorMod <- rma(yi = fisherZDiff, V = allData$seDifference.ro^2, random = ~ authorsTitle.o, mods = ~ source,  data = allData)






