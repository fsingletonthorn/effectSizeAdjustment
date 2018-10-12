library(tidyverse)
library(readxl)
library(MBESS)
library(compute.es)
library(readr)
#source(file = "Data/Functions/any_to_any.R") # Not necessary 

## First extracting data from the OSC's RPP 

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
info <- GET('https://osf.io/fgjvw/?action=download', write_disk('rpp_data.csv', overwrite = TRUE)) #downloads data file from the OSF
MASTER <- read.csv("rpp_data.csv")[1:167, ]
colnames(MASTER)[1] <- "ID" # Change first column name to ID to be able to load .csv file

## Trimming master down
MASTER <- MASTER[(MASTER$Completion..R. == 1)&!is.na(MASTER$Completion..R.),]

###########################################
### Transforming correlations to Fisher ###
###########################################
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

# Removing 3 missing articles 
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

# Calculating standard errors following: 
# http://handbook-5-1.cochrane.org/chapter_7/7_7_7_2_obtaining_standard_errors_from_confidence_intervals_and.htm
# SE = (upper limit – lower limit) / 3.92.

ManyLabs1[c("95CIlb.O", "95CIub.O")]<- str_split(ManyLabs1$`95% CI Lower, Upper.o`, pattern = ", ", simplify = T)
ManyLabs1[c("95CIlb.r", "95CIub.r")]<- str_split(ManyLabs1$`99% CI Lower, Upper.r`, pattern = ", ", simplify = T)
# Transforming into numerics
ManyLabs1[c("95CIlb.O", "95CIub.O","95CIlb.r", "95CIub.r", 'ES.o')] <- sapply(ManyLabs1[c("95CIlb.O", "95CIub.O","95CIlb.r", "95CIub.r", "ES.o")], as.numeric) 

# calculating SEs based on CIs
ManyLabs1['se.o'] <-  ((ManyLabs1["95CIub.O"] - ManyLabs1["95CIlb.O"]))/ 3.92
ManyLabs1['se.r'] <-  ((ManyLabs1["95CIub.r"] - ManyLabs1["95CIlb.r"]))/ 3.92

# calcualting p values following Altman, D. G., & Bland, J. M. (2011). How to obtain the P value from a confidence interval. BMJ, 343.  Retrieved from http://www.bmj.com/content/343/bmj.d2304.abstract
ManyLabs1['pVal.o'] <- exp(-0.717*(ManyLabs1$ES.o/ManyLabs1$se.o) - 0.416*(ManyLabs1$ES.o/ManyLabs1$se.o)^2)

# Extracting avaliable correlations - initializing 
ManyLabs1ML_orig$Correlation <- NA 
# All correlations reported in the paper were .42 (including 2 different studies)
ManyLabs1ML_orig$Correlation[grepl("Anchoring", ManyLabs1ML_orig$Study.name) | grepl("Relations between I and E math attitudes", ManyLabs1ML_orig$Study.name)] <- .42

es <- res(r =  ManyLabs1ML_orig$Correlation, n = ManyLabs1ML_orig$N1, dig = 5)

# SE calculated following Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2011). Introduction to Meta-Analysis. West Sussex, United Kingdom: John Wiley & Sons.  equation 6.4
es$se.z <- ifelse(!is.na(es$r),  sqrt(es$var.r), NA)

########## End Many labs 1 data recollection ########

data2 <- data.frame(authorsTitle.o = ManyLabs1ML_orig$Reference,
                   correlation.o = ManyLabs1ML_orig$Correlation, 
                   cohenD.o = ManyLabs1$ES.o, 
                   fis.o = es$fisher.z, # Fishers z transform
                   seFish.o = es$se.z,
                   n.o = ManyLabs1ML_orig$N1+ManyLabs1ML_orig$N2,
                   seCohenD.o = ManyLabs1$se.o, 
                   pVal.o = ManyLabs1['pVal.o'],
                   resultUsedInRep.o = ManyLabs1ML_orig$Result.used,
                   correlation.r = NA,
                   cohenD.r = ManyLabs1$`Replication ES.r`, 
                   fis.r = es$fisher.z, # Fishers z transform
                   seFish.r = NA,
                   n.r = ManyLabs1$N.r,
                   seCohenD.r = ManyLabs1$se.r,
                   pVal.r = ManyLabs1$p.r,
                   seDifference.ro = NA)

data2$source <- "ManyLabs1"

### 
# Removing everything apart from data sets  
rm(list = c("es","ManyLabs1","ManyLabs1ML_orig"))

##### Many labs 3 data recollection ####
ManyLabs3 <- read_csv("Data/ManyLabs3_Data_ManualAdditions.csv")

# Removing partial eta^2eds  
ManyLabs3 <-ManyLabs3[ManyLabs3$ESstat=="d",]
ManyLabs3$corr.o <- es 

# calcualting p values following Altman, D. G., & Bland, J. M. (2011). How to obtain the P value from a confidence interval. BMJ, 343.  Retrieved from http://www.bmj.com/content/343/bmj.d2304.abstract
ManyLabs3['pVal.o'] <- exp(-0.717*(ManyLabs3$ESOriginal/ManyLabs1$se.o) - 0.416*(ManyLabs1$ES.o/ManyLabs1$se.o)^2)

# converting effect sizes 
es.o <- des(d = ManyLabs3$ESOriginal, n.1 = ManyLabs3$n.o/2, n.2 = ManyLabs3$n.o/2, dig = 5)
es.o$seFish <-ifelse(!is.na(es$r),  sqrt(es$var.r), NA)

# Converting rep ESs 
es.r <- des(d = ManyLabs3$ReplicationES, n.1 = ManyLabs3$N.r/2,  n.2 = ManyLabs3$N.r/2, dig = 5)
es.r$seFish <- ifelse(!is.na(es$r),  sqrt(es$var.r), NA)

# calculating SEs based on CIs
ManyLabs3['se.o'] <-  ((ManyLabs3$ESOriginal95CIUpper -  ManyLabs3$ESOriginal95CILower))/ 3.92
ManyLabs3['se.r'] <-  ((ManyLabs3$Replication95CIUpper - ManyLabs3$Replication95CILower))/ 3.92

# calcualting p values following Altman, D. G., & Bland, J. M. (2011). How to obtain the P value from a confidence interval. BMJ, 343.  Retrieved from http://www.bmj.com/content/343/bmj.d2304.abstract
ManyLabs3['pVal.o'] <- exp(-0.717*(ManyLabs3$ESOriginal/ManyLabs3$se.o) - 0.416*(ManyLabs3$ESOriginal/ManyLabs3$se.o)^2)

data3 <- data.frame(authorsTitle.o = ManyLabs3$originalEffects,
                    correlation.o = es.o$r, 
                    cohenD.o = ManyLabs3$ESOriginal, 
                    fis.o = es.o$fisher.z, 
                    seFish.o = es.o$seFish,
                    n.o = ManyLabs3$n.o,
                    seCohenD.o = ,
                    pVal.o = ,
                    resultUsedInRep.o = ,
                    correlation.r = NA,
                    cohenD.r = ManyLabs3$MedianReplicationES, 
                    fis.r = ,
                    seFish.r = ,
                    n.r = ,
                    seCohenD.r = ,
                    pVal.r = ,
                    seDifference.ro = NA)













# careful with coercction ~ especially of p values some of which are marked as <.001 for example


# View(join_all(list(data, data2), type = 'full'))
### 

tmp <- join_all(list(data, data2), type = 'full')





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















