#### Copied directly from http://www.blakemcshane.com/pcescode.R ####
## Seee McShane, B. B., & Böckenholt, U. (2016). Planning sample sizes when effect sizes are uncertain: The power-calibrated effect size approach. Psychological Methods, 21(1), 47-60. doi:10.1037/met0000036

########################
# Variable Definitions #
########################

# General:
#	alpha - The size of a one-sided test.
#			All formulae assume one-sided tests; for two-sided tests, set alpha = alpha / 2.
#	power - The desired level of power, 1-beta.

# Independent Means:
#	Delta - The difference between the means of the two conditions.
#	sigma - The standard deviation of the individual-level observations.
#	nu - The uncertainty in Delta expressed as a standard deviation.

# Dependent Means:
#	Delta - The difference between the means of the two conditions.
#	sigmaD - The standard deviation of the individual-level differences.
#	nu - The uncertainty in Delta expressed as a standard deviation.

# Independent Proportions:
#	p1 - The proportion in condition one.
#	p2 - The proportion in condition two.
#	Delta - The difference between the proportions of the two conditions, p2-p1.
#	nu - The uncertainty in Delta expressed as a standard deviation.

# Dependent Proportions:
#	p01 - The proportion of the total number of respondents switching from no to yes.
#	p10 - The proportion of the total number of respondents switching from yes to no.
#	p.change - The change proportion, p10/(p01+p10).
#	nu - The uncertainty in p.change expressed as a standard deviation.

# Correlation Coefficient:
#	r - The correlation coefficient.
#	nu - The uncertainty in the Fisher z-transformation of the correlation
#			 Z = arctanh(r) = (1/2)*log((1+r)/(1-r))
#		 which for correlations less than 0.5 in absolute value is essentially identical
#		 the uncertainty in the correlation.

###################################################
# Power-Calibrated Effect Size Formulae (Table 2) #
###################################################
pces.checks <- function(alpha, power, Delta, nu){
  if( alpha >= 0.50 ){ print("We require alpha < 0.50."); return(NA) }
  if( (1-power) >= 0.50 ){ print("We require power > 0.50."); return(NA) }
  if( !(alpha < (1-power) - 1E-08) ){ 
    print("We require alpha < (1-power) for one-sided tests.")
    return(NA)
  }
  if( nu >= abs(Delta)/abs(qnorm(1-power)) ){
    tmp <- "PCES does not exist as uncertainty is too large relative to the
    effect size and the desired level of power."
    print(gsub("[\t\n]", "", tmp))
    return(NA)
  }
  return(1)
}

pces.independent.means <- function(alpha, power, Delta, nu){
  if(Delta==0){ return(Inf) }
  if(Delta < 0){ Delta <- abs(Delta) }
  check <- pces.checks(alpha, power, Delta, nu)
  if(is.na(check)){ return(NA) }
  
  za <- qnorm(1-alpha)	# za = z_{1-\alpha} in the notation of the paper
  zb <- qnorm(1-power)	# zb = z_{\beta} in the notation of the paper	
  pces <- (za*Delta + zb*sqrt(Delta^2 + nu^2*(za^2-zb^2))) / (za+zb)
  pces
}

pces.dependent.means <- pces.independent.means

pces.independent.proportions <- pces.independent.means

pces.dependent.proportions <- function(alpha, power, p.change, nu){
  if(p.change < 0 | p.change > 1){ return(NA) }
  if(p.change==1/2){ return(Inf) }
  if(p.change < 1/2){ p.change <- 1 - p.change }
  
  Delta <- p.change - 0.5
  pces <- 1/2 + pces.independent.means(alpha, power, Delta, nu)
  pces
}

pces.correlation <- function(alpha, power, r, nu){
  if(r < -1 | r > 1){ return(NA) }
  if(r==0){ return(Inf) }
  if(r < 0){ r <- abs(r) }
  
  fisher.z <- atanh(r)
  pces <- pces.independent.means(alpha,power,fisher.z,nu)
  c("correlation"=tanh(pces), "fisher.z"=pces)
}

###########################################
# Textbook Sample Size Formulae (Table 1) #
###########################################
sampsize.independent.means <- function(alpha, power, Delta, sigma){
  if(Delta==0){ return(Inf) }
  if(Delta < 0){ Delta <- abs(Delta) }
  za <- qnorm(1-alpha)	# za = z_{1-\alpha} in the notation of the paper
  zb <- qnorm(1-power)	# zb = z_{\beta} in the notation of the paper
  
  n.per.cell <- 2 * sigma^2 * (za - zb)^2 / Delta^2
  ceiling(n.per.cell)
}

sampsize.dependent.means <- function(alpha, power, Delta, sigmaD){
  sampsize.independent.means(alpha, power, Delta, sigmaD/sqrt(2))
}

sampsize.independent.proportions <- function(alpha, power, p1, p2){
  Delta <- abs(p1 - p2)
  p.bar <- (p1 + p2) / 2
  sampsize.independent.means(alpha, power, Delta, sqrt(p.bar*(1-p.bar)))
}

sampsize.dependent.proportions <- function(alpha, power, p01, p10){
  if( p01 < 0 | p01 > 1 ){ return(NA) }
  if( p10 < 0 | p10 > 1 ){ return(NA) }
  if( p01 + p10 > 1 ){ return(NA) }
  za <- qnorm(1-alpha)	# za = z_{1-\alpha} in the notation of the paper
  zb <- qnorm(1-power)	# zb = z_{\beta} in the notation of the paper
  
  n.per.cell <- (za - zb)^2 * (p01 + p10) / (p01 - p10)^2
  ceiling(n.per.cell)
}

sampsize.correlation <- function(alpha, power, r){
  if(r < -1 | r>1){ return(NA) }
  if(r==0){ return(Inf) }
  if(r < 0){ r <- abs(r) }
  za <- qnorm(1-alpha)	# za = z_{1-\alpha} in the notation of the paper
  zb <- qnorm(1-power)	# zb = z_{\beta} in the notation of the paper
  
  fisher.z <- atanh(r)
  n.per.cell <- 3 + (za - zb)^2 / fisher.z^2
  ceiling(n.per.cell)
}

#########################################
# Power-Calibrated Sample Size Formulae #
#########################################
pces.sampsize.independent.means <- function(alpha, power, Delta, sigma, nu){
  pces <- pces.independent.means(alpha, power, Delta, nu)
  if(is.na(pces)){ return(NA) }
  sampsize.independent.means(alpha, power, pces, sigma)
}

pces.sampsize.dependent.means <- function(alpha, power, Delta, sigmaD, nu){
  pces <- pces.dependent.means(alpha, power, Delta, nu)
  if(is.na(pces)){ return(NA) }
  sampsize.dependent.means(alpha, power, pces, sigmaD)
}

pces.sampsize.independent.proportions <- function(alpha, power, p1, p2, nu){
  tmp <- c(p1,p2)
  p1 <- min(tmp)
  p2 <- max(tmp)
  Delta <- p2-p1
  pces <- pces.independent.proportions(alpha, power, Delta, nu)
  if(is.na(pces)){ return(NA) }
  tmp <- (Delta - pces)/2
  p1.tilde <- p1 + tmp
  p2.tilde <- p2 - tmp
  sampsize.independent.proportions(alpha, power, p1.tilde, p2.tilde)
}

pces.sampsize.dependent.proportions <- function(alpha, power, p01, p10, nu){
  tmp <- c(p01,p10)
  p01 <- min(tmp)
  p10 <- max(tmp)
  p.change <- p10/(p01+p10)	
  pces <- pces.dependent.proportions(alpha, power, p.change, nu)
  if(is.na(pces)){ return(NA) }
  p10.tilde <- (p01 + p10) * pces
  p01.tilde <- (p01 + p10) - p10.tilde
  sampsize.dependent.proportions(alpha, power, p01.tilde, p10.tilde)
}

pces.sampsize.correlation <- function(alpha, power, r, nu){
  pces <- as.numeric( pces.correlation(alpha, power, r, nu)["correlation"] )
  if(is.na(pces)){ return(NA) }
  sampsize.correlation(alpha, power, pces)
}

##############################
# Replicate website defaults #
##############################

# Note: These are also the examples in the Additional Examples subsection.

# General:
alpha <- 0.05
power <- 0.80
nu <- sqrt(0.01)

# Independent Means:
Delta <- 0.2
sigma <- sqrt(1)
pces.independent.means(alpha, power, Delta, nu)
sampsize.independent.means(alpha, power, Delta, sigma)
pces.sampsize.independent.means(alpha, power, Delta, sigma, nu)

# Dependent Means:
Delta <- 0.2
sigmaD <- sqrt(1)
pces.dependent.means(alpha, power, Delta, nu)
sampsize.dependent.means(alpha, power, Delta, sigmaD)
pces.sampsize.dependent.means(alpha, power, Delta, sigmaD, nu)

# Independent Proportions:
p1 <- 0.40
p2 <- 0.60
pces.independent.proportions(alpha, power, p2-p1, nu)
sampsize.independent.proportions(alpha, power, p1, p2)
pces.sampsize.independent.proportions(alpha, power, p1, p2, nu)

# Dependent Proportions:
p01 <- 0.10
p10 <- 0.20
pces.dependent.proportions(alpha, power, p10/(p01+p10), nu)
sampsize.dependent.proportions(alpha, power, p01, p10)
pces.sampsize.dependent.proportions(alpha, power, p01, p10, nu)

# Correlation Coefficient:
r <- 0.2
pces.correlation(alpha, power, r, nu)
sampsize.correlation(alpha, power, r)
pces.sampsize.correlation(alpha, power, r, nu)

#####################################
# Replicate choice overload example #
#####################################

# General:
library(metafor)
alpha <- 0.05
power <- 0.80

# Set data and conduct meta-analysis:
# Study 1	# Study 2								# Study 3
m1  <- c(8.09,		7.81, 									3.81)
sd1 <- c(1.05,		(7.81-7.40)/sqrt(4.18*(1/87+1/78)),		0.54)
n1  <- c(52,		78, 									32)	
m2  <- c(7.69,		7.40,									3.78)
sd2 <- c(0.82,		(7.81-7.40)/sqrt(4.18*(1/87+1/78)),		0.55)
n2  <- c(74,		87,										32)	
s.pool <- sqrt( ((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1+n2-2) )
ma <- rma(measure="SMD", m1i=m1, sd1i=sd1, n1i=n1, m2i=m2, sd2i=sd2, n2i=n2, method="FE")
ma

# Get effect sizes and sample sizes based on the first study:
Delta <- m1[1] - m2[1]
sigma <- s.pool[1]
nu <- sqrt(sigma^2/n1[1] + sigma^2/n2[1])

Delta														# Point estimate
Delta + qnorm(0.20)*nu										# Safeguard Power
pces.independent.means(alpha, power, Delta, nu)				# PCES

sampsize.independent.means(alpha, power, Delta, sigma)
sampsize.independent.means(alpha, power, Delta + qnorm(0.20)*nu, sigma)
sampsize.independent.means(alpha, power, pces.independent.means(alpha, power, Delta, nu), sigma)
pces.sampsize.independent.means(alpha, power, Delta, sigma, nu)

pces.independent.means(alpha, power, Delta, nu) / sigma		# PCES for G*Power

# Get effect sizes and sample sizes based on the meta-analysis:
# We set sigma to one because the meta-analysis is on the standardized mean difference scale.
Delta <- as.numeric(ma$b)
sigma <- 1
nu <- as.numeric(ma$se)

Delta														# Point estimate
Delta + qnorm(0.20)*nu										# Safeguard Power
pces.independent.means(alpha, power, Delta, nu)				# PCES

sampsize.independent.means(alpha, power, Delta, sigma)
sampsize.independent.means(alpha, power, Delta + qnorm(0.20)*nu, sigma)
sampsize.independent.means(alpha, power, pces.independent.means(alpha, power, Delta, nu), sigma)
pces.sampsize.independent.means(alpha, power, Delta, sigma, nu)