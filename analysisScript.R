# source(file = 'Data/data_collection_cleaning.R')

# Calculating the number included in the meta-analysis (also the number included in the study)
nMeta <- sum((!is.na(allData$fisherZDiff)))

# Random effects mod with random effects for authors nested within source
REMod <- rma.mv(yi = fisherZDiff, V = data$seDifference.ro^2, random =  ~ authorsTitle.o|source,  data = allData)
summary(REMod)

REModFixSource <- rma.mv(yi = fisherZDiff, V = data$seDifference.ro^2, random =  ~ authorsTitle.o, mods = ~source,  data = allData)
REModOnlySigR <- rma.mv(yi = fisherZDiff, V = data$seDifference.ro^2, random =  ~ authorsTitle.o|source,  data = allData[allData$significant.r==TRUE,])


REAuthorMod <- rma(yi = fisherZDiff, V = seDifference.ro^2, random = ~ authorsTitle.o, mods = ~ source,  data = allData)

forest(REMod, xlim = c(-1, 1))

summary(REMod)
