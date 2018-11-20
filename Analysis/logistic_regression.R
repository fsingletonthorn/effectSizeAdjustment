# Logistic regression etc. 
# logisticoutcome <- glm(allData$significantSameDirection.r ~ allData$cleanedpVal.o + allData$n.o + allData$fis.o , family = binomial(link = "logit"))

# logisticoutcome <- brms::brm(as.numeric(significantSameDirection.r) ~ cleanedpVal.o +  n.o + fis.o, data = allData, family = bernoulli(link = "logit"), cores = 4)

# logisticoutcomesum <- summary(logisticoutcome)

# fitted.results <- predict(logisticoutcome, type = "response")

# logisticoutcome <- glm(allData$significantSameDirection.r ~ temp[["x.t"]] , family = binomial(link = "logit"))

# fitted.results <- predict(logisticoutcome, type='response')

# fitted.results <- ifelse(fitted.results > 0.5,1,0)

# accuracy <- mean(fitted.results == logisticoutcome$data[2])

### 
