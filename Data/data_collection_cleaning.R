library(tidyverse)
library(readxl)
library(MBESS)
library(compute.es)
library(readr)
library(metafor)
library(ggplot2)

#### First extracting data from the OSC's RPP  ####

########################################################################################
### Code from https://github.com/CenterForOpenScience/rpp/blob/master/masterscript.R####
########################################################################################

# source functions
# if(!require(httr)){install.packages('httr')}
# library(httr)
# info <- GET('https://osf.io/b2vn7/?action=download', write_disk('functions.r', overwrite = TRUE)) #downloads data file from the OSF
# source('functions.r')
if(!require(Hmisc)){install.packages('Hmisc')}
library(Hmisc)
if(!require(metafor)){install.packages('metafor')}
library(metafor)

# Read in Tilburg data for OSF projects 
# info <- GET('https://osf.io/fgjvw/?action=download', write_disk('rpp_data.csv', overwrite = TRUE)) #downloads data file from the OSF - has been done.
MASTER <- read.csv("Data/rpp_data.csv")[1:167, ]
colnames(MASTER)[1] <- "ID" # Change first column name to ID to be able to load .csv file

## Trimming master down
MASTER <- MASTER[(MASTER$Completion..R. == 1)&!is.na(MASTER$Completion..R.),]

### Transforming correlations to Fisher's z
RPP<-data.frame(ri.o  = MASTER$T_r..O.)

RPP$ri.o <- MASTER$T_r..O.
RPP$ri.r <- MASTER$T_r..R.
RPP$N.o <- MASTER$T_df2..O.+2
RPP$N.r <- MASTER$T_df2..R.+2


### Partial correlation, so degrees of freedom plus 2 in order to get N
RPP$N.o[MASTER$ID == 82] <- MASTER$T_df1..O.[MASTER$ID == 82]+2
RPP$N.r[MASTER$ID == 82] <- MASTER$T_df1..R.[MASTER$ID == 82]+2

### Correlation
RPP$N.o[MASTER$ID == 120] <- MASTER$T_N..O.[MASTER$ID == 120]
RPP$N.r[MASTER$ID == 120] <- MASTER$T_N..R.[MASTER$ID == 120]
RPP$N.o[MASTER$ID == 154] <- MASTER$T_N..O.[MASTER$ID == 154]
RPP$N.r[MASTER$ID == 154] <- MASTER$T_N..R.[MASTER$ID == 154]
RPP$N.o[MASTER$ID == 155] <- MASTER$T_N..O.[MASTER$ID == 155]
RPP$N.r[MASTER$ID == 155] <- MASTER$T_N..R.[MASTER$ID == 155]

### t
RPP$N.o[MASTER$ID == 121] <- MASTER$T_N..O.[MASTER$ID == 121]
RPP$N.r[MASTER$ID == 121] <- MASTER$T_N..R.[MASTER$ID == 121]

### Transform to Fisher's z
RPP$fis.o <- 0.5*log((1 + RPP$ri.o) / (1 - RPP$ri.o)) 
RPP$fis.r <- 0.5*log((1 + RPP$ri.r) / (1 - RPP$ri.r))

### Difference in Fisher's z scores
RPP$yi <- 1:length(RPP$fis.o)
for(i in 1:length(RPP$fis.o)) {
  
  if(is.na(RPP$fis.o[i]) == TRUE | is.na(RPP$fis.r[i]) == TRUE) { RPP$yi[i] <- NA }
  else if(RPP$fis.o[i] < 0 & RPP$fis.r[i] < 0) { RPP$yi[i] <- RPP$fis.o[i]*-1-RPP$fis.r[i]*-1 } 
  else if(RPP$fis.o[i] < 0 & RPP$fis.r[i] > 0) { RPP$yi[i] <- RPP$fis.o[i]*-1+RPP$fis.r[i] }
  else {  RPP$yi[i] <- RPP$fis.o[i]-RPP$fis.r[i] }
}

### Standard errors original and replication study
RPP$sei.o <- sqrt(1/(RPP$N.o-3))
RPP$sei.r <- sqrt(1/(RPP$N.r-3))

### p-values original and replication study
RPP$pval.o <- pnorm(RPP$fis.o, sd = RPP$sei.o, lower.tail = FALSE)
RPP$pval.r <- pnorm(RPP$fis.r, sd = RPP$sei.r, lower.tail = FALSE)

### Standard error of difference score
RPP$sei <- sqrt(1/(RPP$N.o-3) + 1/(RPP$N.r-3))

data <- data.frame(authorsTitle.o = paste0(MASTER$Authors..O., "-", MASTER$Study.Title..O.),
                   correlation.o = RPP$ri.o, # Pearson R
                   fis.o = RPP$fis.o, # Fishers z transform
                   seFish.o = RPP$sei.o, # Standard error Fisher trans correlation coef. 
                   n.o = RPP$N.o, # sample size 
                   pValFish.o = RPP$pval.o, # p value from fisher trans r value  
                   pVal.o = MASTER$Reported.P.value..O., # P value from same analysis as replicated study  
                   testStatistic.o  = MASTER$Test.statistic..O., # test stat from original study 
                   correlation.r = RPP$ri.r,
                   fis.r = RPP$fis.r, # Fishers z transform
                   n.r = RPP$N.r,
                   seFish.r = RPP$sei.r,
                   pValFish.r = RPP$pval.r,
                   pVal.r = MASTER$P.value..R., 
                   testStatistic.r = MASTER$Test.statistic..R.,
                   seDifference.ro = RPP$sei)

# Removing 3 missing articles (NS origs)
data <- data[!is.na(RPP$ri.o) & !is.na(RPP$ri.r),]

data$source <- "OSCRPP"

########## End RPP data recollection ########


# Removing everything apart from data from WS 
# rm(list = c("RPP","MASTER","info"))

##############################################
####### Importing dataset MANY LABS 1 ########
##############################################

ManyLabs1 <- read_excel("Data/ManyLabs1_Data.xlsx")
ManyLabs1ML_orig <- read_excel("Data/ManyLabs1_Original_ES95CI.xls")

ManyLabs1ML_orig<-ManyLabs1ML_orig[order(match(ManyLabs1ML_orig$Study.name, ManyLabs1$Effect)),]

ManyLabs1[c("95CIlb.O", "95CIub.O")]<- str_split(ManyLabs1$`95% CI Lower, Upper.o`, pattern = ", ", simplify = T)
ManyLabs1[c("95CIlb.r", "95CIub.r")]<- str_split(ManyLabs1$`99% CI Lower, Upper.r`, pattern = ", ", simplify = T)

# Estimating p values and SE from 95% CIs 


# Transforming into numerics
ManyLabs1[c("95CIlb.O", "95CIub.O","95CIlb.r", "95CIub.r", 'ES.o')] <- sapply(ManyLabs1[c("95CIlb.O", "95CIub.O","95CIlb.r", "95CIub.r", "ES.o")], as.numeric) 

# calculating sample size complete
ManyLabs1$n.o <- ManyLabs1ML_orig$N1 + ManyLabs1ML_orig$N2

# calculating r and fisher z from d 
es.o <- des(d = ManyLabs1$ES.o, n.1 = ManyLabs1$n.o/2, n.2 =  ManyLabs1$n.o/2)

# SE calculated following Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2011). Introduction to Meta-Analysis. West Sussex, United Kingdom: John Wiley & Sons.  equation 6.4
es.o$se.z <- ifelse(!is.na(es.o$r),  sqrt(1/(ManyLabs1$n.o -3)), NA) 

## Removing invalid SEs (studies which use CHI squre stats)
es.o$se.z[str_detect(string = ManyLabs1ML_orig$Result.used,  c("Chi"))] <- NA

# calculating r and fisher z from d 
es.r <- des(d = ManyLabs1$`Replication ES.r`, n.1 = ManyLabs1$N.r/2, n.2 =  ManyLabs1$N.r/2)

# SE calculated following Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2011). Introduction to Meta-Analysis. West Sussex, United Kingdom: John Wiley & Sons.  equation 6.4
es.r$se.z <- ifelse(!is.na(es.r$r),  sqrt(1/(ManyLabs1$N.r -3)), NA) 

## Removing invalid SEs (studies which use CHI square stats)
es.r$se.z[str_detect(string = ManyLabs1$`Key statistics.r`,  "X")] <- NA
es.r$pval.r[str_detect(string = ManyLabs1$`Key statistics.r`,  c("X"))] <- NA

# Extracting p values from original studies
ManyLabs1$pVal.o <- NA

# Extracting the p values from the reported test stat strings
ManyLabs1$pVal.o <- str_split(ManyLabs1$testStatistic.o, "p [=]|[<]", simplify = T)[,2] %>%
  str_extract("......") %>%
  str_remove_all(",|p| |") %>%
  as.numeric()

ManyLabs1$pVal.o[str_detect(ManyLabs1$testStatistic.o, "p < .01")] <- "< .01"
ManyLabs1$pVal.o[str_detect(ManyLabs1$testStatistic.o, "p < .05")] <- "< .05"

   
# Caclulating these for ease p values 
ManyLabs1$pVal.o[ManyLabs1$testStatistic.o == "r = .42"|ManyLabs1$testStatistic.o =="r = .42, n = 243"] <- 
  es.o$pval.r[ManyLabs1$testStatistic.o == "r = .42"|ManyLabs1$testStatistic.o =="r = .42, n = 243"] 
data.frame(ManyLabs1$pVal.o, ManyLabs1$testStatistic.o)


########## End Many labs 1 data recollection ########
data2 <- data.frame(authorsTitle.o = ManyLabs1ML_orig$Reference,
                    correlation.o = es.o$r,
                    cohenD.o = ManyLabs1$ES.o,
                    fis.o = es.o$fisher.z, # Fishers z transform
                    seFish.o = es.o$se.z,
                    pValFish.o = es.o$pval.r,
                    n.o = ManyLabs1$n.o,
                    seCohenD.o = NA, 
                    pVal.o = ManyLabs1$pVal.o,
                    testStatistic.o = ManyLabs1$testStatistic.o,
                    correlation.r = es.r$r,
                    cohenD.r = ManyLabs1$`Replication ES.r`,  
                    seCohenD.r = NA,
                    fis.r = es.r$fisher.z, # Fishers z transform
                    seFish.r = es.r$se.z,
                    n.r = ManyLabs1$N.r,
                    pValFish.r = es.r$pval.r,
                    pVal.r = ManyLabs1$p.r,
                    testStatistic.r = ManyLabs1$`Key statistics.r`,
                    seDifference.ro = NA)

data2$source <- "ManyLabs1"

### 
# Removing everything apart from data sets  
# rm(list = c("es","ManyLabs1","ManyLabs1ML_orig"))

##### Many labs 3 data recollection ####
ManyLabs3 <- read_csv("Data/ManyLabs3_Data_ManualAdditions.csv")

# converting effect sizes 
es.o <- des(d = ManyLabs3$ESOriginal, n.1 = ManyLabs3$n.o/2, n.2 = ManyLabs3$n.o/2, dig = 5)
es.o$seFish <-ifelse(!is.na(es.o$r),  sqrt(1/(ManyLabs3$N.r -3)), NA)

# Converting rep ESs 
es.r <- des(d = ManyLabs3$ReplicationES, n.1 = ManyLabs3$N.r/2,  n.2 = ManyLabs3$N.r/2, dig = 5)
# removing values based on eta-squared values
es.r$r[c(7,9)] <- NA 
es.r$seFish <- ifelse(!is.na(es.r$r),  sqrt(1/(ManyLabs3$n.o -3)), NA)

# removing Boroditsky, L. (2000). Metaphoric structuring: Understanding time through spatial metaphors. Cognition, 75(1), 1-28. who used a chi square test, making SEs for Fisher's z wrong
es.r$pval.r[str_detect(string = ManyLabs3$Effect,  c("MetaphoricRestructuring"))] <- NA
es.o$pval.r[str_detect(string = ManyLabs3$Effect,  c("MetaphoricRestructuring"))] <- NA

# Extracting test stats from rep study
testStats.r <- data.frame(str_split(ManyLabs3$KeyStatistics.r, " =", simplify = T)[,1],
           ManyLabs3$df.r,
           str_split(ManyLabs3$KeyStatistics.r, " =", simplify = T)[,2])
ManyLabs3$testStatistic.r <- str_glue("{testStats.r[,1]} ({testStats.r[,2]}) = {testStats.r[,3]}")



# Amalgomating
data3 <- data.frame(authorsTitle.o = ManyLabs3$originalEffects,
                    correlation.o = es.o$r, 
                    cohenD.o = ManyLabs3$ESOriginal, 
                    seCohenD.o =  NA,
                    fis.o = es.o$fisher.z, 
                    seFish.o = es.o$seFish,
                    n.o = ManyLabs3$n.o,
                    pVal.o = ManyLabs3$p.o,
                    testStatistic.o = ManyLabs3$testStatistic.o,
                    correlation.r = es.r$r,
                    cohenD.r = ManyLabs3$ReplicationES, 
                    seCohenD.r =  NA,
                    fis.r = es.r$fisher.z,
                    seFish.r = es.r$seFish,
                    n.r = ManyLabs3$N.r,
                    pVal.r =  es.r$pval.r,
                    seDifference.ro = NA,
                    testStatistic.r = ManyLabs3$testStatistic.r)

data3$source <- "ManyLabs3"


#### End of Many labs 3 recollection #### 
#### Social science experiments from nature and science Data recollection ####
# data from Camerer, C. F., Dreber, A., Holzmeister, F., Ho, T.-H., Huber, J., Johannesson, M., . . . Wu, H. (2018). Evaluating the replicability of social science experiments in Nature and Science between 2010 and 2015. Nature Human Behaviour, 2(9), 637-644. doi:10.1038/s41562-018-0399-z

natSci <- read_csv(file ="Data/socialScienceExperimentsInNatureAndScience.csv")
# View(natSci)

# calcualting z values and sampling variances 
r1Sums <- escalc(measure = "ZCOR", ri = natSci$r_rs1, ni = natSci$n_rs1)
r2Sums <- escalc(measure = "ZCOR", ri = natSci$r_rs2, ni = natSci$n_rs2)
rOSums <- escalc(measure = "ZCOR", ri = natSci$r_os, ni = natSci$n_os)

natSci$fis.o <- rOSums$yi
natSci$seFish.o <- sqrt(1/(natSci$n_os -3))



summary(r1Sums )

# intitialising data colls
ess <- cbind(yi = c(1,2), xi = c(1,2))
natSci <- natSci
natSci$nTotal.r <- ifelse(!is.na(natSci$n_rs2), natSci$n_rs1 + natSci$n_rs2, natSci$n_rs1)
natSci$fis.r <- NA
natSci$seFish.r <- NA
natSci$pVal.r <- NA

# Running a fixed effects meta-analysis on for each set of studies in the cases when a second study was run
for(i in 1:nrow(r1Sums)) {
  if(!is.na(r2Sums$yi[i])) {
    # Create dataframe for MA 
  ess <- rbind(r1Sums[i,], r2Sums[i,])
   modelout <- rma(yi = ess$yi, vi = ess$vi, method = "FE")
   natSci$fis.r[i] <- modelout$b
   natSci$seFish.r[i] <- modelout$se
   natSci$pVal.r[i]<- modelout$pval 
  } else {
    natSci$fis.r[i] <- r1Sums$yi[i]
    natSci$seFish.r[i] <- r1Sums$vi[i]
    natSci$pVal.r[i]<- natSci$p_rs1[i]
    }
}

# Recalculating corres from orig
natSci$correlation.r <- ztor(natSci$fis.r)

# Excluding SEs developed when original studs were performed using Chi Square tests 
  natSci$seFish.r[str_detect(natSci$type_os, "Chi")] <- NA
  natSci$seFish.o[str_detect(natSci$type_os, "Chi")] <- NA
# Checked that they were not used in the meta-analysis using "natSci$n_rs2[str_detect(natSci$type_os, "Chi")] " which gives 2 NAs, i.e., no second studies were performed for either 
  # cleaning up author names: 
  natSci$Authors <- str_remove_all(natSci$Authors, "\"|[0-9]?[0-9]\\s")
  
  data4 <- data.frame(authorsTitle.o = natSci$Authors,
                      correlation.o = natSci$r_os, 
                      cohenD.o = NA, 
                      fis.o = natSci$fis.o, 
                      seFish.o = natSci$seFish.o,
                      n.o = natSci$n_os,
                      seCohenD.o = NA, 
                      pVal.o = natSci$p_os,
                      testStatistic.o = str_c(natSci$type_os, "=", natSci$stat_os), 
                      correlation.r = natSci$correlation.r,
                      cohenD.r = NA, 
                      fis.r = natSci$fis.r,
                      seFish.r = natSci$seFish.r,
                      n.r = natSci$nTotal.r,
                      seCohenD.r = NA,
                      pVal.r = natSci$pVal.r,
                      seDifference.ro = NA)
  
  data4$source <- "NatSci"

### 
# Removing everything apart from data sets  
# rm(list = c("es.o", "es.r", "ManyLabs1","ManyLabs1ML_orig"))

# careful with coersion ~ especially of p values some of which are marked as <.001 for example


### ####
econ <- read_csv(file = "Data/evaluatingReplicabilityOfLaboratoryExperimentsInEconomics.csv")

es.o <- escalc(ri = econ$original_r, ni = econ$originalN, measure = "ZCOR")
es.r <- escalc(ri = econ$replication_r, ni = econ$replicationN, measure = "ZCOR")

### Standard errors original and replication study - None appear to use techniques for which SEs can be meaningfully extracted
es.o$seFish <- NA# sqrt(1/(econ$originalN-3))
es.r$seFish <- NA# sqrt(1/(econ$replicationN-3))

# View(econ)

# Amalgomating 
data5 <- data.frame(authorsTitle.o = econ$AuthorsJournalYear,
                    correlation.o = econ$original_r, 
                    cohenD.o = NA, 
                    seCohenD.o =  NA,
                    fis.o = es.o$yi, 
                    seFish.o = es.o$seFish,
                    n.o = econ$originalN,
                    pVal.o =  econ$originalPValue,
                    testStatistic.o = NA,
                    correlation.r = econ$replication_r,
                    cohenD.r = NA, 
                    seCohenD.r =  NA,
                    fis.r = es.r$yi,
                    seFish.r = es.r$seFish,
                    n.r = econ$replicationN,
                    pVal.r =  econ$replicationPValue,
                    seDifference.ro = NA)

data5$source <- "Econ"

##### Xphi data recollection #####
# Data from Cova, F., Strickland, B., Abatista, A., Allard, A., Andow, J., Attie, M., . . . Colombo, M. (2018). Estimating the reproducibility of experimental philosophy. Review of Philosophy and Psychology, 1-36. 
xPhi <- read_csv(file ="Data/XPhiReplicability_CompleteData.csv")
es.o <- escalc(ri = xPhi$OriginalRES, ni = xPhi$OriginalN_Effect, measure = "ZCOR")
es.r <- escalc(ri = xPhi$ReplicationRES, ni = xPhi$ReplicationN_Effect, measure = "ZCOR")

### Standard errors original and replication study - None appear to use techniques for which SEs can be meaningfully extracted
es.o$seFish <- sqrt(1/(xPhi$OriginalN_Effect-3))
es.r$seFish <- sqrt(1/(xPhi$ReplicationN_Effect-3))

# extracting original p values
resSplit <- str_split(xPhi$OriginalANALYSIS, "p", simplify = T)
resSplit[21,2] <- resSplit[21,3]
resSplit[36,2] <- .01
pvalues <- str_remove_all(resSplit[,2], " |=|\\)|,|B|1.62|,|Ex|orted") 
pvalues[25] <- NA

 # pnorm(q = 0.02384606, mean = 0, sd = 0.001482005, lower.tail = F)

xPhi$pVal.o <- pvalues
as.numeric(pvalues)

# extracting replication p values
resSplit.r <- str_split(xPhi$ReplicationANALYSIS, "p", simplify = T)
pvalues.r <- str_remove_all(resSplit.r[,2], " |=|, η2=0.007|, B = 0.42, Ex|,η20.007") 


xPhi$pVal.r <- pvalues.r

# Remove SEs from studies with Chi square stats or F degrees of freedom 1 > 1 
es.r$seFish[grepl(x =  xPhi$OriginalANALYSIS, pattern = "χ|X2|F\\(2,|F \\(2")] <- NA
es.o$seFish[grepl(x =  xPhi$OriginalANALYSIS, pattern = "χ|X2|F\\(2,|F \\(2")] <- NA

# Amalgomating 
data6 <- data.frame(authorsTitle.o = xPhi$PAPER_ID,
                    correlation.o = xPhi$OriginalRES, 
                    cohenD.o = NA, 
                    seCohenD.o =  NA,
                    fis.o = es.o$yi, 
                    seFish.o = es.o$seFish,
                    n.o = xPhi$OriginalN_Effect,
                    pVal.o =  xPhi$pVal.o,
                    testStatistic.o = xPhi$OriginalANALYSIS,
                    correlation.r = xPhi$ReplicationRES,
                    cohenD.r = NA, 
                    seCohenD.r =  NA,
                    fis.r = es.r$yi,
                    seFish.r = es.r$seFish,
                    n.r = xPhi$ReplicationN_Effect,
                    pVal.r =  xPhi$pVal.r,
                    seDifference.ro = NA)

data6$source <- "xPhi"

########## End xPhi data recollection #########

loopr <- read_excel('Data/EmbargoFolder/LOOPR_data.xlsx')

# Effect sizes from study
loopr$correlation.o <- as.numeric(loopr$OriginalEffect) 
loopr$correlation.r <- as.numeric(loopr$ReplicationEffect)
loopr$DisattenuatedCorrelation.r <- loopr$DisattenuatedReplicationEffect

# Making correlations positive in the first instance for easily comparibility with the other studies 
loopr$DisattenuatedCorrelationPositive.r <- ifelse(loopr$correlation.o<0, -loopr$DisattenuatedCorrelation.r, loopr$DisattenuatedCorrelation.r)
loopr$correlationPositive.o <- ifelse(loopr$correlation.o<0, -loopr$correlation.o, loopr$correlation.o)

# NOTE THIS USES THE DISATTENUATED CORRELATION for the replication - switched such that the original was always positive  
es.o <- escalc(ri = loopr$correlationPositive.o, ni = loopr$OriginalSampleSize, measure = "ZCOR")
es.r <- escalc(ri = loopr$DisattenuatedCorrelationPositive.r, ni = loopr$ReplicationSampleSize, measure = "ZCOR")

### Standard errors original and replication study - None appear to use techniques for which SEs can be meaningfully extracted
es.o$seFish <- sqrt(1/(loopr$OriginalSampleSize-3))
es.r$seFish <- sqrt(1/(loopr$ReplicationSampleSize-3))

# Removing those that were in fact beta coefficents from consideration 
loopr$DisattenuatedCorrelationPositive.r[loopr$ReplicationEffectType == "B"] <- loopr$correlationPositive.o[loopr$ReplicationEffectType == "B"] <- loopr$DisattenuatedCorrelation.r[loopr$ReplicationEffectType == "B"] <- loopr$correlation.o[loopr$ReplicationEffectType == "B"] <- loopr$correlation.r[(loopr$ReplicationEffectType == "B")] <- es.r[(loopr$ReplicationEffectType == "B"),] <- es.o[(loopr$ReplicationEffectType == "B"),] <- NA

# Amalgomating 
data7 <- data.frame(authorsTitle.o = paste0("Looper", as.numeric(gsub("([0-9]+).*$", "\\1", loopr$OutcomeNumber))),
                    correlation.o = loopr$correlationPositive.o, 
                    cohenD.o = NA, 
                    seCohenD.o =  NA,
                    fis.o = es.o$yi, 
                    seFish.o = es.o$seFish,
                    n.o = loopr$OriginalSampleSize,
                    pVal.o =  loopr$OriginalPValue,
                    testStatistic.o = NA,
                    correlation.r = loopr$DisattenuatedCorrelationPositive.r,
                    # correlationNotDisttenuated.r = loopr$DisattenuatedCorrelation.r, 
                    cohenD.r = NA, 
                    seCohenD.r =  NA,
                    fis.r = es.r$yi,
                    seFish.r = es.r$seFish,
                    n.r = loopr$ReplicationSampleSize,
                    pVal.r = loopr$ReplicationPValue,
                    seDifference.ro = NA,
                    testStatistic.r = NA)

data7$source <- "LOOPR"

# Bringing it all together
allData <- plyr::join_all(list(data, data2, data3, data4, data5, data6, data7), type = 'full')

# Replacing one value which seems to be have an effect size in it too
allData$pVal.r[ which(allData$pVal.r == "0.092,η20.007")] <- '0.092'

# Setting binary for significant / not replication
allData$significant.r  <- (as.numeric(allData$pVal.r) <.05 | is.na(as.numeric(allData$pVal.r)))

# Prepping for MLM - calculating SEs and variance
allData$fisherZDiff <- allData$fis.r - allData$fis.o
allData$seDifference.ro <- sqrt(1/(allData$n.o-3) + 1/(allData$n.r-3))

### Finding the 
### Finding minimium effect that would have been significant in the original study - 
minimumEffectDetectable <- qnorm(.05, mean = 0, sd = tmp$seFish.o, lower.tail = FALSE)
upper95.r <- tmp$fis.r + ( 1.96 * tmp$seFish.r )
lower95.r <- tmp$fis.r - ( 1.96 * tmp$seFish.r )

# metafor 
mean(qnorm(.05, mean = 0, sd = tmp$seFish.o, lower.tail = FALSE), na.rm =T )

