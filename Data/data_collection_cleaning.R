library(tidyverse)
library(readxl)
library(MBESS)
library(compute.es)
library(readr)
library(metafor)

#source(file = "Data/Functions/any_to_any.R") # Not necessary, but possible later to aid translation b/w effect sizes with reasonable SEs being developed

#### First extracting data from the OSC's RPP  ####

########################################################################################
### Code from https://github.com/CenterForOpenScience/rpp/blob/master/masterscript.R####
########################################################################################

# source functions
if(!require(httr)){install.packages('httr')}
library(httr)
info <- GET('https://osf.io/b2vn7/?action=download', write_disk('functions.r', overwrite = TRUE)) #downloads data file from the OSF
source('functions.r')
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
RPP$N.o[MASTER$ID == 82] <- MASTER$T_df1..O.[82]+2
RPP$N.r[MASTER$ID == 82] <- MASTER$T_df1..R.[82]+2

### Correlation
RPP$N.o[MASTER$ID == 120] <- MASTER$T_N..O.[120]
RPP$N.r[MASTER$ID == 120] <- MASTER$T_N..R.[120]
RPP$N.o[MASTER$ID == 154] <- MASTER$T_N..O.[154]
RPP$N.r[MASTER$ID == 154] <- MASTER$T_N..R.[154]
RPP$N.o[MASTER$ID == 155] <- MASTER$T_N..O.[155]
RPP$N.r[MASTER$ID == 155] <- MASTER$T_N..R.[155]

### t
RPP$N.o[MASTER$ID == 121] <- MASTER$T_N..O.[121]
RPP$N.r[MASTER$ID == 121] <- MASTER$T_N..R.[121]

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
                   pVal.o = RPP$pval.o, # P value 
                   resultUsedInRep.o  = MASTER$Test.statistic..O.,
                   correlation.r = RPP$ri.o,
                   fis.r = RPP$fis.r, # Fishers z transform
                   n.r = RPP$N.r,
                   seFish.r = RPP$sei.r,
                   pVal.r = RPP$pval.r,
                   seDifference.ro = RPP$sei)

# Removing 3 missing articles (NS origs)
data <- data[!is.na(RPP$ri.o) & !is.na(RPP$ri.r),]

data$source <- "OSCRPP"

########## End RPP data recollection ########


# Removing everything apart from data from WS 
rm(list = c("RPP","MASTER","info"))

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

# calculating r and fisher z from d 
es.o <- des(d = ManyLabs1$ES.o, n.1 = ManyLabs1$n.o/2, n.2 =  ManyLabs1$n.o/2)

# SE calculated following Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2011). Introduction to Meta-Analysis. West Sussex, United Kingdom: John Wiley & Sons.  equation 6.4
es.o$se.z <- ifelse(!is.na(es.o$r),  sqrt(1/(ManyLabs1$n.o -3)), NA) 

## Removing invalid SEs (studies which use CHI squre stats)
es.o$se.z[str_detect(string =ManyLabs1ML_orig$Result.used,  c("Chi"))] <- NA

# calculating r and fisher z from d 
es.r <- des(d = ManyLabs1$`Replication ES.r`, n.1 = ManyLabs1$N.r/2, n.2 =  ManyLabs1$N.r/2)

# SE calculated following Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2011). Introduction to Meta-Analysis. West Sussex, United Kingdom: John Wiley & Sons.  equation 6.4
es.r$se.z <- ifelse(!is.na(es.r$r),  sqrt(1/(ManyLabs1$N.r -3)), NA) 

## Removing invalid SEs (studies which use CHI square stats)
es.r$se.z[str_detect(string = ManyLabs1$`Key statistics.r`,  "X")] <- NA
es.r$pval.r[str_detect(string = ManyLabs1$`Key statistics.r`,  c("X"))] <- NA

########## End Many labs 1 data recollection ########
data2 <- data.frame(authorsTitle.o = ManyLabs1ML_orig$Reference,
                    correlation.o = es.o$r,
                    cohenD.o = ManyLabs1$ES.o,
                    fis.o = es.o$fisher.z, # Fishers z transform
                    seFish.o = es.o$se.z,
                    n.o = ManyLabs1$n.o,
                    seCohenD.o = NA, 
                    pVal.o = es.o$pval.r,
                    resultUsedInRep.o = ManyLabs1ML_orig$Result.used,
                    correlation.r = es.r$r,
                    cohenD.r = ManyLabs1$`Replication ES.r`,  
                    seCohenD.r = NA,
                    fis.r = es.r$fisher.z, # Fishers z transform
                    seFish.r = es.r$se.z,
                    n.r = ManyLabs1$N.r,
                    pVal.r = es.r$pval.r,
                    seDifference.ro = NA)

data2$source <- "ManyLabs1"

### 
# Removing everything apart from data sets  
rm(list = c("es","ManyLabs1","ManyLabs1ML_orig"))

##### Many labs 3 data recollection ####
ManyLabs3 <- read_csv("Data/ManyLabs3_Data_ManualAdditions.csv")

# Removing partial eta^2eds  !!!
ManyLabs3 <-ManyLabs3[ManyLabs3$ESstat=="d",]

# converting effect sizes 
es.o <- des(d = ManyLabs3$ESOriginal, n.1 = ManyLabs3$n.o/2, n.2 = ManyLabs3$n.o/2, dig = 5)
es.o$seFish <-ifelse(!is.na(es.o$r),  sqrt(1/(ManyLabs3$N.r -3)), NA)

# Converting rep ESs 
es.r <- des(d = ManyLabs3$ReplicationES, n.1 = ManyLabs3$N.r/2,  n.2 = ManyLabs3$N.r/2, dig = 5)
es.r$seFish <- ifelse(!is.na(es.r$r),  sqrt(1/(ManyLabs3$n.o -3)), NA)

# removing Boroditsky, L. (2000). Metaphoric structuring: Understanding time through spatial metaphors. Cognition, 75(1), 1-28. who used a chi square test, making SEs for Fisher's z wrong
es.r$pval.r[str_detect(string = ManyLabs3$Effect,  c("MetaphoricRestructuring"))] <- NA
es.r$pval.r[str_detect(string = ManyLabs3$Effect,  c("MetaphoricRestructuring"))] <- NA
es.o$pval.r[str_detect(string = ManyLabs3$Effect,  c("MetaphoricRestructuring"))] <- NA
es.o$pval.r[str_detect(string = ManyLabs3$Effect,  c("MetaphoricRestructuring"))] <- NA

# Amalgomating 
data3 <- data.frame(authorsTitle.o = ManyLabs3$originalEffects,
                    correlation.o = es.o$r, 
                    cohenD.o = ManyLabs3$ESOriginal, 
                    seCohenD.o =  NA,
                    fis.o = es.o$fisher.z, 
                    seFish.o = es.o$seFish,
                    n.o = ManyLabs3$n.o,
                    pVal.o = es.o$pval.r,
                    resultUsedInRep.o = NA,
                    correlation.r = es.r$r,
                    cohenD.r = ManyLabs3$ReplicationES, 
                    seCohenD.r =  NA,
                    fis.r = es.r$fisher.z,
                    seFish.r = es.r$seFish,
                    n.r = ManyLabs3$N.r,
                    pVal.r =  es.r$pval.r,
                    seDifference.ro = NA)

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

natSciOutput$fis.o <- rOSums$yi
natSciOutput$seFish.o <- sqrt(rOSums$vi)

summary(r1Sums )

# intitialising data colls
ess <- cbind(yi = c(1,2), xi = c(1,2))
natSciOutput <- natSci
natSciOutput$nTotal.r <- ifelse(!is.na(natSci$n_rs2), natSci$n_rs1 + natSci$n_rs2, natSci$n_rs1)


# Running a fixed effects meta-analysis on for each set of studies in the cases when a second study was run
for(i in 1:nrow(r1Sums)) {
  if(!is.na(r2Sums$yi[i])) {
    # Create dataframe for MA 
  ess <- rbind(r1Sums[i,], r2Sums[i,])
   modelout <- rma(yi = ess$yi, vi = ess$vi, method = "FE")
   natSciOutput$fis.r[i] <- modelout$b
   natSciOutput$seFish.r[i] <- modelout$se
   natSciOutput$pVal.r[i]<- modelout$pval 
  } else {
    natSciOutput$fis.r[i] <- r1Sums$yi[i]
    natSciOutput$seFish.r[i] <- r1Sums$vi[i]
    natSciOutput$pVal.r[i]<- natSciOutput$p_rs1[i]
    }
}

# Recalculating corres from orig
natSciOutput$correlation.r <- ztor(natSciOutput$fis.r)

# Excluding SEs developed when original studs were performed using Chi Square tests 
  natSciOutput$seFish.r[str_detect(natSciOutput$type_os, "Chi")] <- NA
  natSciOutput$seFish.o[str_detect(natSciOutput$type_os, "Chi")] <- NA
# Checked that they were not used in the meta-analysis using "natSciOutput$n_rs2[str_detect(natSciOutput$type_os, "Chi")] " which gives 2 NAs, i.e., no second studies were performed for either 
  
  
  data4 <- data.frame(authorsTitle.o = NA,
                      correlation.o = natSciOutput$r_os, 
                      cohenD.o = NA, 
                      fis.o = natSciOutput$fis.o, 
                      seFish.o = natSciOutput$seFish.o,
                      n.o = natSciOutput$n_os,
                      seCohenD.o = NA, 
                      pVal.o = natSciOutput$p_os,
                      resultUsedInRep.o = str_c(natSciOutput$type_os, "=", natSciOutput$stat_os), 
                      correlation.r = natSciOutput$correlation.r,
                      cohenD.r = NA, 
                      fis.r = natSciOutput$fis.r,
                      seFish.r = natSciOutput$seFish.r,
                      n.r = natSciOutput$nTotal.r,
                      seCohenD.r = NA,
                      pVal.r = natSciOutput$pVal.r,
                      seDifference.ro = NA)

### 
# Removing everything apart from data sets  
# rm(list = c("es.o", "es.r", "ManyLabs1","ManyLabs1ML_orig"))

# careful with coersion ~ especially of p values some of which are marked as <.001 for example

# View(join_all(list(data, data2), type = 'full'))




  
  
  





### 



tmp <- join_all(list(data, data2, data3), type = 'full')

datax <- data.frame(authorsTitle.o = ,
                    correlation.o = , 
                    cohenD.o = , 
                    fis.o = , 
                    seFish.o = ,
                    n.o = ,
                    seCohenD.o = , 
                    pVal.o = ,
                    resultUsedInRep.o = ,
                    correlation.r = ,
                    cohenD.r = , 
                    fis.r = ,
                    seFish.r = ,
                    n.r = ,
                    seCohenD.r = ,
                    pVal.r = ,
                    seDifference.ro = )

# LATER - convert effect sizes and extract SEs - possibility of using 









## https://osf.io/z7aux/
## HAVE TO GO THROUGH AND REMOVE THOSE BASED ON Effect size statistics based on F(df1> 1, df2) and χ2(df) can be converted to correlations (see A3), but their standard errors cannot be computed. 






