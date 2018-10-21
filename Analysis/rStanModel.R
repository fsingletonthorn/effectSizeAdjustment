library(rstan)
library(brms)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

allData$study <- as.factor(allData$authorsTitle.o)

mod <- brm(fisherZDiff | se(seDifference.ro) ~ 1 + (1|source/study), data = allData[!is.na(allData$seDifference.ro) & !is.na(allData$fisherZDiff),], cores = 4, control =list(adapt_delta = .99, max_treedepth = 15) )


mod <- brm(fisherZDiff | se(seDifference.ro) ~ 1 + (1|study|source), data = allData[!is.na(allData$seDifference.ro) & !is.na(allData$fisherZDiff),], cores = 4, control =list(adapt_delta = .99, max_treedepth = 15) )

summary(mod)

pairs(mod)

fit <- stan_demo("eight_schools", refresh = 0)

vignette(package = "brms")



zinb <- read.csv("http://stats.idre.ucla.edu/stat/data/fish.csv")
zinb$camper <- factor(zinb$camper, labels = c("no", "yes"))
head(zinb)

fit_zinb1 <- brm(count ~ persons + child + camper, data = zinb,
                 family = zero_inflated_poisson("log"))

