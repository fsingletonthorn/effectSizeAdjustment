library(tidyverse)
options(stringsAsFactors = F)

# Data and lables from table 2 of the supplmentary materials of Camerer, C. F., Dreber, A., Forsell, E., Ho, T.-H., Huber, J., Johannesson, M., . . . Wu, H. (2016). Evaluating replicability of laboratory experiments in economics. Science, 351(6280), 1433.  Retrieved from http://science.sciencemag.org/content/351/6280/1433.abstract


lables <- c("StudyRef", "originalPValue", "original_r", "originalN", "replicationPValue", "replication_r", "replicationN", "replicationOriginalRatio", "unstandardisedRatio")

table2 <- c("Abeler et al. (AER 2011) 33 0.046 0.18 120 0.160 0.08 318 No 0.43 (0.36)",
"Ambrus and Greiner (AER 2012) 34 0.057 0.31 117 0.012 0.23 357 Yes 0.74 (0.69)",
"Bartling et al. (AER 2012) 35 0.007 0.72 216 0.001 0.66 360 Yes 0.91 (1.21)",
"Charness and Dufwenberg (AER 2011) 36 0.010 0.38 162 0.003 0.36 264 Yes 0.95 (0.9)",
"Chen and Chen (AER 2011) 37 0.033 0.84 72 0.571 0.17 168 No 0.20 (0.22)",
"de Clippel et al. (AER 2014) 38 0.001 0.12 158 <0.001 0.27 156 Yes 2.26 (3.21)",
"Duffy and Puzzello (AER 2014) 39 0.010 0.76 54 0.674 0.12 96 No 0.15 (-0.19)",
"Dulleck et al. (AER 2011) 40 <0.001 0.72 168 0.001 0.73 128 Yes 1.01 (0.94)",
"Ericson and Fuster (QJE 2011) 41 0.030 0.21 112 0.055 0.12 262 No 0.58 (0.69)",
"Fehr et al. (AER 2013) 42 0.011 0.45 60 0.026 0.31 102 Yes 0.69 (0.84)",
"Friedman and Oprea (AER 2012) 43 <0.001 0.64 78 0.004 0.44 40 Yes 0.68 (0.68)",
"Fudenberg et al. (AER 2012) 44 0.001 0.30 124 <0.001 0.33 128 Yes 1.08 (0.96)",
"Huck et al. (AER 2011) 45 0.004 0.83 120 0.142 0.37 160 No 0.44 (0.43)",
"Ifcher and Zarghamee (AER 2011) 46 0.031 0.28 58 0.933 0.01 131 No 0.02 (-0.02)",
"Kessler and Roth (AER 2012) 47 <0.001 0.49 288 0.016 0.34 48 Yes 0.71 (0.62)",
"Kirchler et al (AER 2012) 48 0.016 0.66 120 0.010 0.53 220 Yes 0.80 (0.30)",
"Kogan et al. (AER 2011) 49 <0.001 0.32 126 0.001 0.30 90 Yes 0.94 (0.93)",
"Kuziemko et al. (QJE 2014) 50 0.070 0.28 42 0.154 0.12 144 No 0.42 (-0.39)")

# Removing whitesapce * 4 from the above
splitTable2 <- str_split(table2, pattern = "\\) ", n = 2, simplify = T)
splitTable2[,1] <- str_remove_all(splitTable2[,1], " ")
dat <- str_split(splitTable2[,2], " ", simplify = T)

data.frame(splitTable2[,1], dat)
