library(tidyverse)
library(readxl)
library(MBESS)
library(psych)
library(compute.es)
library(readr)
library(metafor)
library(ggplot2)
library(knitr)
# devtools::install_github('hollylkirk/ochRe')
library(ochRe)
# Setting up a dataframe of names of the source studies for easy plotting later

projectNames <- data_frame(c("OSC (2015) \n General Psychology",
                                                     "Klein et al. (2014)\n Many Labs 1", 
                                                     "Ebersole et al. (2016)\n Many Labs 3",
                                                     "Camerer, et al. (2018)\n Nature Science",
                                                     "Camerer et al. (2016)\n Economics",
                                                     "Cova, et al. (2018)\n Experimental Philosophy",
                                                     "Soto (2019)\n Personality Psychology",
                                                     "Klein et al. (2018)\n Many Labs 2"))


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

### ADDITION FELIX - Removing SEs for Fdf_1 > 2 and chi square
RPP$sei.o[str_detect(MASTER$Test.statistic..O., "F\\([2-9]|F\\(1[1-9]|X\\^2|b|z")] <- NA
RPP$sei.r[str_detect(MASTER$Test.statistic..R., "F\\([2-9]|F\\(1[1-9]|X\\^2|b|z")] <- NA

data <- data_frame(authorsTitle.o = paste0(MASTER$Authors..O., "-", MASTER$Study.Title..O.),
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
                   pVal.r = as.character(MASTER$P.value..R.), 
                   testStatistic.r = as.character(MASTER$Test.statistic..R.),
                   seDifference.ro = RPP$sei)

# Removing 3 articles with missing effect sizes, and two NS origs
data <- data[!is.na(RPP$ri.o) & !is.na(RPP$ri.r),]
data <- data[-which( data$authorsTitle.o=="PW Eastwick, EJ Finkel-Sex differences in mate preferences revisited: Do people know what they initially desire in a romantic partner?" | 
                       data$authorsTitle.o=='KA Ranganath, BA Nosek-Implicit attitude generalization occurs immediately; explicit attitude generalization takes time'),]

# pvalue for missing data replication study 
data$pVal.r[data$pVal.r=="2.2 x 10-16"] <- 2.2 * 10^16
data$pVal.r[data$pVal.r=="prep > .99"] <- "< .0005" # estimated using psych::p.rep(.0005)
data$pVal.r[data$pVal.r=="0"] <- "< .001"
data$pVal.r[data$pVal.r=="X"] <- data$pValFish.r[data$pVal.r=="X"] # this one is a correlation, this p value should be accurate

# Checking that all NAs are non-significant 
# data$pVal.r[is.na(as.numeric(data$pVal.r))]


data$source <- as.character(projectNames[1,1])
  #"OSC (2015)"
data$abrev <- "OSCRPP"
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

# Transforming into numerics
ManyLabs1[c("95CIlb.O", "95CIub.O","95CIlb.r", "95CIub.r", 'ES.o')] <- sapply(ManyLabs1[c("95CIlb.O", "95CIub.O","95CIlb.r", "95CIub.r", "ES.o")], as.numeric) 

# calculating sample size complete
ManyLabs1$n.o <- ManyLabs1ML_orig$N1 + ManyLabs1ML_orig$N2

# calculating r and fisher z from d (accounting for unequal sample sizes)
es.o.equal <- data.frame(des(d = ManyLabs1$ES.o, n.1 = ManyLabs1$n.o/2, n.2 =  ManyLabs1$n.o/2))
es.o.unequal <- data.frame(des(d = ManyLabs1$ES.o, n.1 = ManyLabs1ML_orig$N1, n.2 =  ManyLabs1ML_orig$N2))
es.o <- es.o.unequal 

# Bringing the above together
for(i in 1:nrow(es.o.equal)) { if(es.o$n.2[i] == 0 & !is.na(es.o$n.2[i]))  es.o[i,] <- es.o.equal[i,] }

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

data2$source <- as.character(projectNames[2,1])
data2$abrev <- "ML1"
  #"Klein et al. (2014), Many Labs 1"

# Checking that all NAs are non-significant 
# data2$pVal.r[is.na(as.numeric(as.character(data2$pVal.r)))]


### 
# Removing everything apart from data sets  
# rm(list = c("es","ManyLabs1","ManyLabs1ML_orig"))

##### Many labs 3 data recollection ####
ManyLabs3 <- read_csv("Data/ManyLabs3_Data_ManualAdditions.csv")

# converting effect sizes 
es.o <- des(d = ManyLabs3$ESOriginal, n.1 = ManyLabs3$n.o/2, n.2 = ManyLabs3$n.o/2, dig = 5)
es.o$seFish <-ifelse(!is.na(es.o$r),  sqrt(1/(ManyLabs3$n.o -3)), NA)

# Converting rep ESs 
es.r <- des(d = ManyLabs3$ReplicationES, n.1 = ManyLabs3$N.r/2,  n.2 = ManyLabs3$N.r/2, dig = 5)
# removing values based on eta-squared values
es.r$r[c(7,9)] <- NA 
es.r$seFish <- ifelse(!is.na(es.r$r),  sqrt(1/(ManyLabs3$N.r -3)), NA)

# removing Boroditsky, L. (2000). Metaphoric structuring: Understanding time through spatial metaphors. Cognition, 75(1), 1-28. who used a chi square test, making SEs for Fisher's z wrong
es.r$seFish[str_detect(string = ManyLabs3$Effect,  c("MetaphoricRestructuring"))] <- NA
es.o$seFish[str_detect(string = ManyLabs3$Effect,  c("MetaphoricRestructuring"))] <- NA

# Inputting an effect size from a t test for an original study 

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

data3$source <- as.character(projectNames[3,1])
  #"Ebersole et al. (2016), Many Labs 3"
data3$abrev <- "ML3"


# Checking that all NAs are non-significant 
# data3$pVal.r[is.na(as.numeric(as.character(data3$pVal.r)))]

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
  
data4$source <- as.character(projectNames[4,1])
  #"Camerer, et al. (2018), Nature Science"
data4$abrev <- "natSci"

  # Checking that all NAs are non-significant 
  # data4$pVal.r[is.na(as.numeric(as.character(data4$pVal.r)))]
  

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
                    testStatistic.o = econ$originalTest,
                    correlation.r = econ$replication_r,
                    cohenD.r = NA, 
                    seCohenD.r =  NA,
                    fis.r = es.r$yi,
                    seFish.r = es.r$seFish,
                    n.r = econ$replicationN,
                    pVal.r =  econ$replicationPValue,
                    seDifference.ro = NA)

data5$source <- as.character(projectNames[5,1])
  #  "Camerer et al. (2016), Economics"
data5$abrev <- "econ"

# Checking that all NAs are non-significant 
# data5$pVal.r[is.na(as.numeric(as.character(data5$pVal.r)))]


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


xPhi$pVal.r[str_detect(xPhi$pVal.r, "20.007")] <- 0.092


# EXCLUDING STUDIES WHICH WERE NOT SIGNIFICANT 
data.frame(xPhi$pVal.o,  as.numeric(xPhi$pVal.o) > .05 & !is.na(as.numeric(xPhi$pVal.o)))

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


data6$source <- as.character(projectNames[6,1])
  #"Cova, et al. (2018), Experimental Philosophy"
data6$abrev <- "xPhi"
# Checking that NA p values are significant 
data6$pVal.r[is.na(as.numeric(as.character(data6$pVal.r)))]

########## End xPhi data recollection #########
########## Begin LOOPR data collection ############
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

# Removing additional invalid SEs 

es.r$seFish[str_detect(loopr$OriginalAnalysis,"Structural equation model")] <- es.o$seFish[str_detect(loopr$OriginalAnalysis,"Structural equation model")] <- NA

# Amalgomating 
data7 <- data.frame(authorsTitle.o = loopr$OriginalStudyCitation,
                    correlation.o = loopr$correlationPositive.o, 
                    cohenD.o = NA, 
                    seCohenD.o =  NA,
                    fis.o = es.o$yi, 
                    seFish.o = es.o$seFish,
                    n.o = loopr$OriginalSampleSize,
                    pVal.o =  loopr$OriginalPValue,
                    testStatistic.o = loopr$OriginalAnalysis,
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

data7$source <- as.character(projectNames[7,1])
data7$abrev <- "loopr"
  #"Soto (2019), Personality Psychology"

# Checking that NA p values are significant 
# data7$pVal.r[is.na(as.numeric(as.character(data7$pVal.r)))]


####### end loopr data collection########
##### Begin ManyLabs  data collection #####

ml2 <- read_xlsx("Data/EmbargoFolder/ML2_Data_extracted_from_tables.xlsx")

## Removing six studies for which effect sizes cannot be derrived (i.e., 4 which are diffs b/w effects, only es for each reported, and 2 which were cohen's q's)
ml2 <- ml2[-which(is.na(as.numeric(ml2$Original_effect_size_d))),]

# estimating rs and fishers z values 
es.o <- des(as.numeric(ml2$Original_effect_size_d), n.1 = ml2$Original_study_sample_size/2, n.2 = ml2$Original_study_sample_size/2)
es.o$seFish <- sqrt(1/(ml2$Original_study_sample_size-3))

es.r <-  des(as.numeric(ml2$Replication_global_effect_size), n.1 = ml2$Replication_sample_size/2, n.2 = ml2$Replication_sample_size/2)
es.r$seFish  <- sqrt(1/(ml2$Replication_sample_size-3))

# Removing invalid SEs
es.r$seFish[str_detect(ml2$OriginalTest, 'Chi')] <- es.o$seFish[str_detect(ml2$OriginalTest, 'Chi')] <- NA

# Amalgomating
data8 <- data.frame(authorsTitle.o = ml2$ArticleRef,
                    correlation.o = es.o$r, 
                    cohenD.o = es.o$d, 
                    seCohenD.o =  NA,
                    fis.o = es.o$fisher.z, 
                    seFish.o = es.o$seFish,
                    n.o = ml2$Original_study_sample_size,
                    pVal.o =  NA,
                    testStatistic.o = ml2$OriginalTest,
                    correlation.r = es.r$r,
                    cohenD.r = es.r$d, 
                    seCohenD.r =  NA,
                    fis.r = es.r$fisher.z,
                    seFish.r = es.r$seFish,
                    n.r = es.r$N.total,
                    pVal.r = as.character(ml2$`replication_sample_size,_p_<_.0001`),
                    seDifference.ro = NA,
                    testStatistic.r = ml2$ReplicationTest)

data8$source <- as.character(projectNames[8,1])
  #"Klein et al. (2018), Many Labs 2"
data8$abrev <- "ML2"


# Checking that NA p values are significant 
# data8$pVal.r[is.na(as.numeric(as.character(data8$pVal.r)))]

# Bringing it all together
allData <- plyr::join_all(list(data, data2, data3, data4, data5, data6, data7, data8), type = 'full')

# Setting binary for significant / not replication
allData$significant.r  <- (as.numeric(allData$pVal.r) <.05 | is.na(as.numeric(allData$pVal.r)))

# Setting binary for significant and in the same direction / not replication
allData$significantSameDirection.r  <- (as.numeric(allData$pVal.r) <.05 | is.na(as.numeric(allData$pVal.r))) & (sign(allData$correlation.o) == sign(allData$correlation.r))

# Setting binary for significant and in the same direction / not replication
allData$significantSameDirection.r  <- (as.numeric(allData$pVal.r) <.05 | is.na(as.numeric(allData$pVal.r)))& sign(allData$correlation.o)==sign(allData$correlation.r)
allData$significantSameDirection.r[is.na(allData$significantSameDirection.r)] <- TRUE

# Prepping for MLM - calculating SEs and variance
allData$fisherZDiff <- allData$fis.r - allData$fis.o
allData$seDifference.ro <- sqrt(1/(allData$n.o-3) + 1/(allData$n.r-3))

# For plotting - correlation coefficent difference 
allData$correlationDifference.ro <- allData$correlation.r - allData$correlation.o

# percentage change from original to eventual effect size
allData$percentageChangeES.ro <-  (allData$fis.r - allData$fis.o)/allData$fis.o

projectNamesSingleLine <- data_frame(unique(allData$source), c("OSC (2015)",
                                                               "Klein et al. (2014), Many Labs 1", 
                                                               "Ebersole et al. (2016), Many Labs 3",
                                                               "Camerer, et al. (2018), Nature Science",
                                                               "Camerer et al. (2016), Economics",
                                                               "Cova, et al. (2018), Experimental Philosophy",
                                                               "Soto (2019), Personality Psychology",
                                                               "Klein et al. (2018), Many Labs 2"))
