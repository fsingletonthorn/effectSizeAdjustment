## Level (1-a) Wilson confidence interval for proportion
## WILSON, E. B. 1927. Probable inference, the law of succession,
## and statistical inference. Journal of the American Statistical
## Association 22: 209-212.
#adapted from https://www.r-bloggers.com/reference-chart-for-precision-of-wilson-binomial-proportion-confidence-interval/
# n is the absolute number of paricipants, p is the proportion

WilsonBinCI <-  function(n, p, a=0.05) {
  z <- qnorm(1-a/2,lower.tail=FALSE)
  l <- 1/(1+1/n*z^2)*(p + 1/2/n*z^2 + 
                        z*sqrt(1/n*p*(1-p) + 1/4/n^2*z^2))
  u <- 1/(1+1/n*z^2)*(p + 1/2/n*z^2 -
                        z*sqrt(1/n*p*(1-p) + 1/4/n^2*z^2))
list(lower=l, upper=u, width=u-l)
}

# this version requires you to provide a vector x, and specify which values indicate success
WilsonBinCIVec <-  function(x, markerForSuccess, a=0.05) {
  z <- qnorm(1-a/2,lower.tail=FALSE)
  n <- length(x)
  p <- mean(x == markerForSuccess)
  l <- 1/(1+1/n*z^2)*(p + 1/2/n*z^2 + 
                        z*sqrt(1/n*p*(1-p) + 1/4/n^2*z^2))
  u <- 1/(1+1/n*z^2)*(p + 1/2/n*z^2 -
                        z*sqrt(1/n*p*(1-p) + 1/4/n^2*z^2))
  list(lower=l, upper=u, width=u-l)
}


