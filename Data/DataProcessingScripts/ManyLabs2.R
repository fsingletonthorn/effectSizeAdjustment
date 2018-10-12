# many labs two data extraction

library( R.utils)
install.packages("Data\\DataProcessingScripts\\manylabRs_0.1.0.tar.gz", repos = NULL, type="source", verbose = T, clean = T)

sourceDirectory("Data/Functions")



ML2.key     <- get.GoogleSheet(data='ML2masteRkey')
source("Data/Functions/ML2_variable_functions.R") 
source("Data/Functions/manylabRs.R") 
source("Data/Functions/inIT.R") 
source("Data/Functions/fRedsRutils.R") 
source("Data/Functions/C-3PR_ASCII.R") 


load("~/PhD/Effect size adjustment testing paper/Data/EmbargoFolder/ML2_S1.rda")
load("~/PhD/Effect size adjustment testing paper/Data/EmbargoFolder/ML2_S2.rda")
View(get.analyses)
