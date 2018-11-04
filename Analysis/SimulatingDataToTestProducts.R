simData <- allData[c("correlation.o", "fis.o", "n.o", "n.r","seFishAprox.o" ,"seFish.o", "seFish.r" ,"seFishAprox.r", "seDifference.ro", "pVal.o", "source", "authorsTitle.o")]

# Must run data cleaning and analysis scripts before this one 
nSim <- 10000

# Tracking vectors
simEquivProp  <- rep(NA, length(nSim))
simPropTrueEffect <- rep(NA, length(nSim))
descriptivesTrackingSim <- data.frame(means=rep(NA, length(tableReductions$`n included`)*nSim), model=rep(NA, length(tableReductions$`n included`)*nSim), 
                                      nTrue=rep(NA, length(tableReductions$`n included`)*nSim), nValid=rep(NA, length(tableReductions$`n included`)*nSim),
                                      meanPropChange = rep(NA, length(tableReductions$`n included`)*nSim), modelEstimate = rep(NA, length(tableReductions$`n included`)*nSim),
                                      simNum =  rep(NA, length(tableReductions$`n included`)*nSim), propNull =  rep(NA, length(tableReductions$`n included`)*nSim),
                                      propAttenuation = rep(NA, length(tableReductions$`n included`)*nSim), trueMeanDifference = rep(NA, length(tableReductions$`n included`)*nSim))

# This controls whether you are adding or starting over again in the simulation - THIS HAS TO BE FALSE FOR THE FIRST ROUND 
add <- TRUE

for(i in 1:nSim) {
  propNull <- sample(x = seq(0,1,by=.1), 1)
  propAttenuation <- sample(x = seq(0,1,by=.1), 1)
# Simulating selection from the true effect - i.e., adding random variability as the se of the original (i.e. finding the "true" effect) and replciation studies, 
# this value assumes that if there were no Ns reported, the SE was actually the mean of all of the standard errors

simData$simulatedFish.true <- (rnorm(nrow(simData), simData$fis.o, ifelse(!is.na(simData$seFish.o),simData$seFish.o, mean(simData$seFish.o, na.rm = T) ))) * (1- propAttenuation)

simData$simulatedFish.true[!is.na(simData$simulatedFish.true)][sample(x = 1:sum(!is.na(simData$simulatedFish.true)), 
                                                                      size =   round(sum(!is.na(simData$simulatedFish.true))*propNull), replace = F)] <- 0
simData$fis.r <- simData$simulatedFish.true + (rnorm(nrow(simData), 0, ifelse(!is.na(simData$seFish.r),simData$seFish.r, mean(simData$seFish.r, na.rm = T)) ))

simData$correlation.r <- ztor(simData$fis.r)

# calculating differences
simData$fisherZDiffSim <- simData$fis.r - simData$fis.o

trueMeanDiff <- mean(simData$simulatedFish.true - simData$fis.r, na.rm = T)

trueMeanDiffNo0s <- mean(simData$simulatedFish.true[simData$simulatedFish.true!=0&!is.na(simData$simulatedFish.true)] - simData$fis.r[simData$simulatedFish.true!=0&!is.na(simData$simulatedFish.true)] )

simData$p.r <- pnorm(q = simData$fis.r, mean = 0, sd = simData$seFishAprox.r, lower.tail = F)

simData$significantSameDirection.r <- (simData$p.r<.05 & simData$fis.r > 0)

# Two sided default Bayes factor
simData$bf10Sim <- apply(data.frame(simData$correlation.r, simData$n.r), MARGIN = 1, FUN = bfapply10) 
simData$bf01Sim <- 1/simData$bf10Sim 

# one sided ~ probably makes more sense as they have all been shifted to positive direction
simData$bfplus0Sim <- apply(data.frame(simData$correlation.r, simData$n.r), MARGIN = 1, FUN = bfapplyplus0) 
simData$bf0plusSim <- 1/simData$bfplus0Sim 

# Replication Bayes factor
simData$bfrep0Sim <- apply(data.frame(simData$correlation.o, simData$n.o, simData$correlation.r, simData$n.r), MARGIN = 1, FUN = bfapplyRep0) 
simData$bf0repSim <- 1/simData$bfrep0Sim

simData$bf10Sim.o <-  apply(data.frame(simData$correlation.o, simData$n.o), MARGIN = 1, FUN = bfapply10)

# mean(simData$bfrep0Sim_updated - simData$bfrep0Sim, na.rm =T)
# plot(log(simData$bfrep0Sim_updated), log(simData$bfrep0Sim))
# ## Finding minimium effect that would have been significant in the original study - 
minimumEffectDetectableZ <- qnorm(.05, mean = 0, sd = simData$seFishAprox.o, lower.tail = FALSE)
CIs90UB <-   simData$fis.r + (qnorm(0.95)*simData$seFishAprox.r)
CIs90LB <-   simData$fis.r - (qnorm(0.95)*simData$seFishAprox.r)

# calculating the number of equiv studies
simData$simStatisticallyEquiv.ro <- CIs90UB < minimumEffectDetectableZ

######################################
#### Amount of change in subsets #####
######################################
# Reduction in all studies: 
reductionsimData <- simData %>% 
  filter(!is.na(fisherZDiffSim)) %>% 
  summarise(meanPropChange = mean((fisherZDiffSim)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T), 
            meanDiff = mean(fisherZDiffSim, na.rm = T), sdDiff = sd(fisherZDiffSim, na.rm = T),    
            medianPropChange = median((fisherZDiffSim)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiffSim, na.rm = T),
            nTrue = sum(!is.na(fisherZDiffSim)), nValid = sum(!is.na(simData$fis.r - simData$fis.o)))
# extracting amount of change in subsets
# Sig results replication only (same direction)
reductionSignifcantR <- simData %>% 
  filter(significantSameDirection.r & !is.na(significantSameDirection.r)) %>% 
  summarise(meanPropChange = mean((fisherZDiffSim)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T), 
            meanDiff = mean(fisherZDiffSim, na.rm = T), sdDiff = sd(fisherZDiffSim, na.rm = T),    
            medianPropChange = median((fisherZDiffSim)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiffSim, na.rm = T),
            nTrue = sum(!is.na(significantSameDirection.r)), nValid = sum(!is.na(simData$significantSameDirection.r) & !is.na(simData$fis.r - simData$fis.o)))
# BF01 < 3 (exclude those with moderate or greater evidence for the null)
reductionSbf01LessThan3 <- simData %>% 
  filter(bf01Sim<3) %>% 
  summarise(meanPropChange = mean((fisherZDiffSim)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiffSim, na.rm = T), sdDiff = sd(fisherZDiffSim, na.rm = T),   
            medianPropChange = median((fisherZDiffSim)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiffSim, na.rm = T),
            nTrue = sum(!is.na(simStatisticallyEquiv.ro)), nValid = sum(!is.na(simData$bf01Sim) & !is.na(simData$fis.r - simData$fis.o)))
# bf10Sim > 3 (Only those with evidence for the alternative)
reductionSbf10SimMoreThan3 <- simData %>% 
  filter(bf10Sim>3) %>% 
  summarise(meanPropChange = mean((fisherZDiffSim)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) , 
            meanDiff = mean(fisherZDiffSim, na.rm = T), sdDiff = sd(fisherZDiffSim, na.rm = T),    
            medianPropChange = median((fisherZDiffSim)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiffSim, na.rm = T),
            nTrue = sum(!is.na(simStatisticallyEquiv.ro)), nValid = sum(!is.na(simData$bf10Sim) & !is.na(simData$fis.r - simData$fis.o)))
# bf0plusSim < 3 (exclude those with moderate evidence for the null, one sided alternative)
reductionSbf0plusSimLessThan3 <- simData %>% 
  filter(bf0plusSim<3) %>% 
  summarise(meanPropChange = mean((fisherZDiffSim)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiffSim, na.rm = T), sdDiff = sd(fisherZDiffSim, na.rm = T),   
            medianPropChange = median((fisherZDiffSim)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiffSim, na.rm = T),
            nTrue = sum(!is.na(simStatisticallyEquiv.ro)), nValid = sum(!is.na(simData$bf0plusSim)& !is.na(simData$fis.r - simData$fis.o)))
# bfplus0Sim (Exclude those without moderate evidence for the one sided alternative)
reductionSbfplus0SimMoreThan3 <- simData %>% 
  filter(bfplus0Sim>3) %>%
  summarise(meanPropChange = mean((fisherZDiffSim)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiffSim, na.rm = T), sdDiff = sd(fisherZDiffSim, na.rm = T),    
            medianPropChange = median((fisherZDiffSim)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiffSim, na.rm = T),
            nTrue = sum(!is.na(simStatisticallyEquiv.ro)), nValid = sum(!is.na(simData$bfplus0Sim) & !is.na(simData$fis.r - simData$fis.o)))
# bf0repSim < 3 (exclude those with moderate evidence for the null, replication alternative)
reductionSbf0repSimLessThan3 <- simData %>% 
  filter(bf0repSim<3) %>% 
  summarise(meanPropChange = mean((fisherZDiffSim)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiffSim, na.rm = T), sdDiff = sd(fisherZDiffSim, na.rm = T),  
            medianPropChange = median((fisherZDiffSim)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiffSim, na.rm = T),
            nTrue = sum(bf0repSim<3, na.rm = T), nValid = sum(!is.na(simData$bf0repSim) & !is.na(simData$fis.r - simData$fis.o)))
# bfrep0Sim (Exclude those without moderate evidence for the replication alternative)
reductionSbfrep0SimMoreThan3 <- simData %>% 
  filter(bfrep0Sim>3) %>%
  summarise(meanPropChange = mean((fisherZDiffSim)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiffSim, na.rm = T), sdDiff = sd(fisherZDiffSim, na.rm = T),  
            medianPropChange = median((fisherZDiffSim)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiffSim, na.rm = T),
            nTrue = sum(!is.na(bfrep0Sim>3)), nValid = sum(!is.na(simData$bfrep0Sim) & !is.na(simData$fis.r - simData$fis.o)))
# Significant TOST (Exclude those without moderate evidence for the one sided alternative)
reductionEquiv <- simData %>% 
  filter(!simStatisticallyEquiv.ro) %>%
  summarise(meanPropChange = mean((fisherZDiffSim)/fis.o, na.rm = TRUE), mean.r = mean(fis.r, na.rm = T) , mean.o = mean(fis.o, na.rm = T) ,
            meanDiff = mean(fisherZDiffSim, na.rm = T), sdDiff = sd(fisherZDiffSim, na.rm = T),   
            medianPropChange = median((fisherZDiffSim)/fis.o, na.rm = TRUE), median.r = median(fis.r, na.rm = T) , median.o = median(fis.o, na.rm = T) ,
            medianDiff = median(fisherZDiffSim, na.rm = T),
            nTrue = sum(!is.na(simStatisticallyEquiv.ro)), nValid = sum(!is.na(simData$simStatisticallyEquiv.ro)& !is.na(simData$fis.r - simData$fis.o)))


# Bringing all of the above together:
tableReductionsSim <- rbind(Overall = reductionsimData, StatisticalSignificance = reductionSignifcantR, Nonequivalence = reductionEquiv, bf0repSimBelow3 = reductionSbf0repSimLessThan3 , bfrep0SimAbove3 = reductionSbfrep0SimMoreThan3,  BF01Below3 = reductionSbf01LessThan3 , bf10SimAbove3 = reductionSbf10SimMoreThan3, BF0PBelow3 = reductionSbf0plusSimLessThan3, BFP0Above3 = reductionSbfplus0SimMoreThan3)

######################################
###### Multilevel meta-analysis ######
######################################
# Placeholder for error catching 
errorDF <- data.frame(k = NA, b = NA, ci.lb = NA, ci.ub = NA)

# Random effects model with random effects for authors nested within source
REMod <- tryCatch(rma.mv(yi = fisherZDiffSim, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o,  data = simData), error = function(e) errorDF)
# The first model but with only significant replications
REModOnlySigR <- tryCatch(rma.mv(yi = fisherZDiffSim, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = simData[simData$significantSameDirection.r==TRUE & !is.na(simData$significantSameDirection.r),]), error = function(e) errorDF)

# the first model removing all studies with BFs0Plus > 3 (i.e., moderate evidence for a null)
REModbf0plusSimGreaterThan3Excluded <- tryCatch(rma.mv(yi = fisherZDiffSim, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = simData[simData$bf0plusSim < 3 & !is.na(simData$bf0plusSim > 3),]), error = function(e) errorDF)

# the first model removing all studies with BFsPlus0 < 3 (i.e., those without evidence for the alternative)
REModbfplus0SimLessThan3Excluded <- tryCatch(rma.mv(yi = fisherZDiffSim, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = simData[simData$bfplus0Sim > 3 & !is.na(simData$bfplus0Sim > 3),]), error = function(e) errorDF)

# the first model removing all studies with BFs01 > 3 (i.e., moderate evidence for a null)
REModBF01GreaterThan3Excluded <- tryCatch(rma.mv(yi = fisherZDiffSim, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = simData[simData$bf01Sim < 3 & !is.na(simData$bf01Sim > 3),]), error = function(e) errorDF)

# the first model removing all studies with BFs10 < 3 (i.e., those without evidence for the alternative)
REModbf10SimLessThan3Excluded <- tryCatch(rma.mv(yi = fisherZDiffSim, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = simData[simData$bf10Sim > 3 & !is.na(simData$bf10Sim > 3),]), error = function(e) errorDF)

# the first model removing all studies with BFs0Rep < 3 (i.e., those without evidence for the for the null vs. origianl greater than 3)
REModbf0repSimGreaterThan3Excluded <-  tryCatch(rma.mv(yi = fisherZDiffSim, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = simData[simData$bf0repSim < 3 & !is.na(simData$bf0repSim < 3),]), error = function(e) errorDF)

# the first model removing all studies with BFs0Rep < 3 (i.e., those without evidence for the for the null vs. origianl greater than 3)
REModbfrep0SimLessThan3Excluded <- tryCatch(rma.mv(yi = fisherZDiffSim, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data =  simData[simData$bfrep0Sim > 3 & !is.na(simData$bfrep0Sim > 3),]), error = function(e) errorDF)

# The first model but with only non-equiv 
REModNonequiv <- tryCatch(rma.mv(yi = fisherZDiffSim, V = seDifference.ro^2, random =  ~ 1|source/authorsTitle.o, data = simData[simData$simStatisticallyEquiv.ro==FALSE & !is.na(simData$simStatisticallyEquiv.ro==FALSE),]), error = function(e) errorDF)

modResSim <- data.frame( modelN = REMod$k, modelEstimate = REMod$b, MLM95lb = REMod$ci.lb, MLM95ub = REMod$ci.ub, row.names = "Overall")
modResSim1 <- data.frame(modelN = REModOnlySigR$k, modelEstimate = REModOnlySigR$b, MLM95lb = REModOnlySigR$ci.lb, MLM95ub = REModOnlySigR$ci.ub, row.names = "StatisticalSignificance")
modResSim2 <- data.frame(modelN = REModbf0plusSimGreaterThan3Excluded$k, modelEstimate = REModbf0plusSimGreaterThan3Excluded$b, MLM95lb = REModbf0plusSimGreaterThan3Excluded$ci.lb, MLM95ub = REModbf0plusSimGreaterThan3Excluded$ci.ub, row.names = "BF0PBelow3")
modResSim3 <- data.frame(modelN = REModbfplus0SimLessThan3Excluded$k, modelEstimate = REModbfplus0SimLessThan3Excluded$b, MLM95lb = REModbfplus0SimLessThan3Excluded$ci.lb, MLM95ub = REModbfplus0SimLessThan3Excluded$ci.ub, row.names = "BFP0Above3")
modResSim4 <- data.frame(modelN = REModBF01GreaterThan3Excluded$k, modelEstimate = REModBF01GreaterThan3Excluded$b, MLM95lb = REModBF01GreaterThan3Excluded$ci.lb, MLM95ub = REModBF01GreaterThan3Excluded$ci.ub, row.names = "BF01Below3")
modResSim5 <- data.frame(modelN = REModbf10SimLessThan3Excluded$k, modelEstimate = REModbf10SimLessThan3Excluded$b, MLM95lb = REModbf10SimLessThan3Excluded$ci.lb, MLM95ub = REModbf10SimLessThan3Excluded$ci.ub, row.names = "bf10SimAbove3")
modResSim6 <- data.frame(modelN = REModbf0repSimGreaterThan3Excluded$k, modelEstimate = REModbf0repSimGreaterThan3Excluded$b, MLM95lb = REModbf0repSimGreaterThan3Excluded$ci.lb, MLM95ub = REModbf0repSimGreaterThan3Excluded$ci.ub, row.names = "bf0repSimBelow3")
modResSim7 <- data.frame(modelN = REModbfrep0SimLessThan3Excluded$k, modelEstimate = REModbfrep0SimLessThan3Excluded$b, MLM95lb = REModbfrep0SimLessThan3Excluded$ci.lb, MLM95ub = REModbfrep0SimLessThan3Excluded$ci.ub, row.names = "bfrep0SimAbove3")
modResSim8 <- data.frame(modelN = REModNonequiv$k, modelEstimate = REModNonequiv$b, MLM95lb = REModNonequiv$ci.lb, MLM95ub = REModNonequiv$ci.ub, row.names = "Nonequivalence")

modSumariesSim <- rbind(modResSim, modResSim1, modResSim2, modResSim3, modResSim4, modResSim5, modResSim6, modResSim7, modResSim8)

modSumariesSim$simNum <- i
modSumariesSim$propNull <- propNull
modSumariesSim$propAttenuation <- propAttenuation
modSumariesSim$trueMeanDifference <- trueMeanDiff
modSumariesSim$trueMeanDifferenceNo0s <- trueMeanDiffNo0s


if(i == 1 && add == FALSE) {
tableAllEstimatesSim <- merge.data.frame(tableReductionsSim, modSumariesSim, by = "row.names", sort = F)
} else{
simOutput <-  merge.data.frame(tableReductionsSim, modSumariesSim, by = "row.names", sort = F)
tableAllEstimatesSim <- data.table::rbindlist(list(tableAllEstimatesSim, simOutput), use.names = T, fill = T , idcol = FALSE)
print(paste(i, "out of", nSim))
}
}
write.csv(tableAllEstimatesSim, file = "Data/SimulationModelOutput.csv", row.names = FALSE)

tableAllEstimatesSim$squaredError <- (-tableAllEstimatesSim$propAttenuation-tableAllEstimatesSim$meanPropChange)^2
tableAllEstimatesSim$absoluteError <- abs(-tableAllEstimatesSim$propAttenuation-tableAllEstimatesSim$meanPropChange)

simulationSum <-as.tibble(tableAllEstimatesSim) %>%
group_by(propNull, propAttenuation, Row.names) %>%
  dplyr::summarise(meanMeanPropChange = mean(meanPropChange), sdMeanPropChange = sd(meanPropChange), meanChange = mean(meanDiff), 
                   meanModelEstimate = mean(modelEstimate), sdModelEstimate = sd(modelEstimate),  nSims = n(), MSE = mean(squaredError), 
                   RMSE = sqrt(mean(squaredError)),  MAE = mean(absoluteError), meanTrueDiff = mean(trueMeanDifference, na.rm = T),
                   meanTrueDiff = mean(trueMeanDifferenceNo0s, na.rm = T), errorMeanPropReduction = mean(-propAttenuation - meanPropChange))

vis <- simulationSum 

# Plotting distance between estimated mean change and true prop
ggplot(vis, aes(x = propAttenuation, y = propNull)) +
  geom_raster(aes(fill = errorMeanPropReduction), interpolate=F)  +
  scale_fill_gradient2(low="navy", mid="white", high="red", 
                       midpoint=0)+ facet_wrap(~ Row.names) + theme_bw()

# Plotting SD by method at different values 
ggplot(vis, aes(x = propAttenuation, y = propNull)) +
  geom_raster(aes(fill = sdMeanPropChange), interpolate=F)  +
  scale_fill_gradient2(low="navy", mid="white", high="red", 
                       midpoint=0)+ facet_wrap(~ Row.names) + theme_bw()

# Plotting MSE by method at different values 
ggplot(vis, aes(x = propAttenuation, y = propNull)) +
  geom_raster(aes(fill = MSE), interpolate=F)  +
  scale_fill_gradient2(low="navy", mid="white", high="red", 
                       midpoint=0)+ facet_wrap(~ Row.names) + theme_bw()

# Plotting MAE by method at different values 
ggplot(vis, aes(x = propAttenuation, y = propNull)) +
  geom_raster(aes(fill = MAE), interpolate=F)  +
  scale_fill_gradient2(low="navy", mid="white", high="red", 
                       midpoint=0)+ facet_wrap(~ Row.names) + theme_bw()

simulationSumByType <- as.tibble(tableAllEstimatesSim) %>%
  group_by(Row.names) %>%
  dplyr::summarise(meanMeanPropChange = mean(meanPropChange, na.rm = T), sdMeanPropChange = sd(meanPropChange, na.rm = T), 
                   meanModelEstimate = mean(modelEstimate, na.rm = T), sdModelEstimate = sd(modelEstimate, na.rm = T),  nSims = n(), MSE = mean(squaredError, na.rm = T), 
                   RMSE = sqrt(mean(squaredError, na.rm = T)),  MAE = mean(absoluteError, na.rm = T), AE_SD = sd(absoluteError, na.rm = T), Error_SD = sd(-propAttenuation - meanPropChange, na.rm = T))
View(simulationSumByType)



simulationSumByTypeLessThan75 <- as.tibble(tableAllEstimatesSim) %>%
  group_by(Row.names, below.8s = propNull < .8 & propAttenuation < .8) %>%
  dplyr::summarise(meanMeanPropChange = mean(meanPropChange, na.rm = T), sdMeanPropChange = sd(meanPropChange, na.rm = T), 
                   meanModelEstimate = mean(modelEstimate, na.rm = T), sdModelEstimate = sd(modelEstimate, na.rm = T),  nSims = n(), MSE = mean(squaredError, na.rm = T), 
                   RMSE = sqrt(mean(squaredError, na.rm = T)),  MAE = mean(absoluteError, na.rm = T), AE_SD = sd(absoluteError, na.rm = T), Error_SD = sd(-propAttenuation - meanPropChange, na.rm = T))
View(simulationSumByTypeLessThan75)


# Plotting mean distance between estimated mean change and true prop
ggplot(vis, aes(x = propAttenuation, y = propNull)) +
  geom_raster(aes(fill = meanModelEstimate/mean(allData$fis.o, na.rm =T)), interpolate=F)  +
  scale_fill_gradient2(low="navy", mid="white", high="red", 
                       midpoint=0)+ facet_wrap(~ Row.names) + theme_bw()


# Plotting mean distance between estimated mean change and true prop
ggplot(vis, aes(x = propAttenuation, y = propNull)) +
  geom_raster(aes(fill = meanMeanPropChange), interpolate=F)  +
  scale_fill_gradient2(low="navy", mid="white", high="red", 
                       midpoint=0, limits = c(-1,0))+ facet_wrap(~ Row.names) + theme_bw()

