library(tidyverse)
colLabs <- c("Effect",
"ESstat",
"ESOriginal",
"ESOriginal95CILower",
"ESOriginal95CIUpper",
"MedianReplicationES",
"ReplicationES",
"Replication95CILower",
"Replication95CIUpper",
"ReplicationESMeta",
"ReplicationESCILowerMeta",
"ReplicationESCIUpperMeta",
"ProportionLessThan05",
"ProportionGreaterThan05",
"ProportionNS",
"KeyStatistics",
"df",
"N",
"p")
options(stringsAsFactors = F)

table3Data <- c("Stroop Task d 2.04 1.90, 2.18 0.89 0.91 0.84, 0.98 0.88 0.84, 0.92 0.00 1.00 0.00 t = 49.795 3336 3337 < .001",
"Metaphoric Restructuring d 0.63 0.07, 1.20 0.19 0.29 0.17, 0.42 0.26 0.15, 0.37 0.00 0.25 0.75 χ2 = 21.90 1 1335 < .001",
"Availability Heuristic d 0.82 0.47, 1.17 0.07 0.09 0.02, 0.16 0.09 0.02, 0.16 0.00 0.14 0.86 PSdep = .522 N/A 3088 0.015",
"Power and Perspective d 0.77 0.12, 1.41 0.03 0.03 -0.04, 0.11 0.03 -0.04, 0.10 0.00 0.05 0.95 t = .89 2967 2969 0.37",
"Weight Embodiment d 0.59 0.01, 1.16 0.05 0.03 -0.05, 0.11 0.03 -0.06, 0.11 0.00 0.00 1.00 t = .61 2283 2285 0.543",
"Warmth Perceptions d 0.86 0.40, 1.33 0.06 0.01 -0.08, 0.06 0.01 -0.06, 0.08 0.05 0.00 0.95 t = .22 3117 3119 0.827",
"Elaboration Likelihood ηp² 0.17 0.06, 0.29 0.000001 0.000001 0.000, 0.001 0.00005 0.000, 0.002 0.00 0.00 1.00 F = .129 1, 2361 2365 0.72",
"Self-Esteem and Subjective Distance ηp² n/a n/a 0.0001 0.0004 0.000, 0.005 0.001 0.000, 0.004 0.05 0.10 0.85 F = 1.98 1, 3131 3136 0.16",
"Credentials and Prejudice ηp² 0.04 0, 0.09 0.00 0.0003 0.000, 0.003 0.000 0.000, 0.000 0.00 0.10 0.90 F = .0004 1, 3130 3134 0.985")

table3Dat <- str_replace(table3Data, " ", "")
table3Dat[c(4, 8, 9)] <- str_replace(table3Dat[c(4, 8, 9)], " ", "")
table3Dat[c(8)] <- str_replace(table3Dat[c(8)], " ", "")
table3Dat <- str_replace(table3Dat, " = ", "")
table3Dat <- str_replace(table3Dat, "< ", "")
table3Dat <- str_replace(table3Dat, " 1, ", " df1=1,df2=")
table3Dat <- str_replace(table3Dat, ", ", " ")
table3Dat <- str_replace(table3Dat, "n/a", "NA")
table3Dat <- str_replace(table3Dat, "N/A", "NA")
table3Dat <- str_replace(table3Dat, "n/a", "NA NA")
table3Dat <- str_replace(table3Dat, ",", "")
table3 <- data.frame(str_split(table3Dat, " ", simplify = T))
names(table3) <- colLabs

write.csv(table2, file = "", na = "NA")
