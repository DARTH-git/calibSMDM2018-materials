#########     Calibration of Sick-Sicker Markov model           #####################
#####################################################################################

#####################################################################################
# Please cite our publications when using this code
# - Alarid-Escudero F, Maclehose RF, Peralta Y, Kuntz KM, Enns EA. 
# Non-identifiability in model calibration and implications for medical decision making. 
# Med Decis Mak. 2018;Forthcoming. 
# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, Hunink MG. 
# An Overview of R in Health Decision Sciences. Med Decis Making. 2017; 37(3): 735-746. 
# - Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22. 

#####################################################################################

rm(list = ls())      # clear memory (removes all the variables from the workspace)

#### 01 Load packages ####
library(plotrix)
library(lhs)
library(scales)
library(psych)
library(MHadaptive)
library(matrixStats)
source("../Functions.R")

#### 02 Load targets ####
load("SickSickerTargets.RData")

#### 02.1 Plot targets ####
### Survival
plotrix::plotCI(x = 1:30, y = l.targets$Surv$value, 
                ui = l.targets$Surv$ub,
                li = l.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
### Prevalence
plotrix::plotCI(x = 1:30, y = l.targets$Prev$value, 
                ui = l.targets$Prev$ub,
                li = l.targets$Prev$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Proportion Sick and Sicker")
### Proportion Sick
plotrix::plotCI(x = c(10, 20, 30), y = l.targets$PropSick$value, 
                ui = l.targets$PropSick$ub,
                li = l.targets$PropSick$lb,
                xlab = "Time", ylab = "Proportion Sick")

#### 03 Function of Sick-Sicker Markov model ####
#### 03.1 Load Function ####
source("Markov_Sick-Sicker - Function.R")
#### 03.2 Test model ####
v.params0 <- c(p.S1S2 = 0.105, hr.S1 = 3, hr.S2 = 10)
markov_sick_sicker(v.params0) # It works!

#### 04 Weighted Sum of square errors function ####
w_sse <- function(params){
  if(is.null(dim(params))) { # If vector, change to matrix
    params <- t(params) 
  }
  wsse <- rep(0, nrow(params))
  for(j in 1:nrow(params)) { # j=1
    jj <- tryCatch( { 
      # Obtain model outputs from a set of parameters
      res_j <- markov_sick_sicker(v.params = params[j, ])
      # Survival => Normal likelihood
      wsse[j] <- sum(((l.targets$Surv$value - res_j$Surv)/l.targets$Surv$se)^2)
      # Disease prevalence => Normal likelihood
      wsse[j] <- wsse[j] + sum(((l.targets$Prev$value - res_j$Prev)/l.targets$Prev$se)^2)
      # Proportion Sick of total sick  => Normal likelihood
      wsse[j] <- wsse[j] + sum(((l.targets$PropSick$value - res_j$PropSick)/l.targets$PropSick$se)^2)
    }, error = function(e) NA) 
    if(is.na(jj)) { wsse[j] <- -Inf }
  } 
  return(wsse)
}
# Test it
w_sse(v.params0) # It works!

#### 05 Define Log-Likelihood function ####
l_likelihood <- function(params){
  if(is.null(dim(params))) { # If vector, change to matrix
    params <- t(params) 
  }
  llik <- rep(0, nrow(params))
  for(j in 1:nrow(params)) { # j=1
    jj <- tryCatch( { 
      # Obtain model outputs from a set of parameters
      res_j <- markov_sick_sicker(v.params = params[j, ])
      # Survival => Normal likelihood
      llik[j] <- sum(dnorm(x = l.targets$Surv$value, 
                           mean = res_j$Surv, 
                           sd = l.targets$Surv$se, 
                           log = T))
      # Disease prevalence => Normal likelihood
      llik[j] <- llik[j] + sum(dnorm(x = l.targets$Prev$value, 
                                     mean = res_j$Prev, 
                                     sd = l.targets$Prev$se, 
                                     log = T))
      # Proportion of sick in S1 state => Normal likelihood
      llik[j] <- llik[j] + sum(dnorm(x = l.targets$PropSick$value, 
                                     mean = res_j$PropSick, 
                                     sd = l.targets$PropSick$se, 
                                     log = T))
    }, error = function(e) NA) 
    if(is.na(jj)) { llik[j] <- -Inf }
  } 
  return(llik)
}
# Test it
l_likelihood(v.params0) # It works!

#### 06 Function to generate a Latin hypercube sample (LHS) ####
gen_lhs <- function(n.samp = 1000, seed = 071918){
  set.seed(seed)              # set a seed to be able to reproduce the same results
  ### Lower bounds
  lb <- c(p.S1S2 = 0.01, hr.S1 = 1.0, hr.S2 = 5)
  ### Upper bounds
  ub <- c(p.S1S2 = 0.50, hr.S1 = 4.5, hr.S2 = 15)
  
  ### Latin hypercube sampling
  lhs.samp <- randomLHS(n.samp, 3)
  
  df.lhs <- data.frame(
    p.S1S2 = qunif(lhs.samp[, 1],
                   min = lb["p.S1S2"],  # Parameters of desired distribution
                   max = ub["p.S1S2"]), # Parameters of desired distribution
    hr.S1 = qunif(lhs.samp[, 2],
                  min = lb["hr.S1"],  # Parameters of desired distribution
                  max = ub["hr.S1"]), # Parameters of desired distribution
    hr.S2 = qunif(lhs.samp[, 3],
                  min = lb["hr.S2"], # Parameters of desired distribution
                  max = ub["hr.S2"]) # Parameters of desired distribution
  )
  return(df.lhs)
}
# Test it
gen_lhs(n.samp = 10) # It works!
pairs.panels(gen_lhs(n.samp = 1000))

#### 07 Find MLE of LHS ###
# Number of samples to generate
n.samp <- 10000
# generate Latin hypercube sample
lhs.grid  <- gen_lhs(n.samp)
# loop over grid samples
v.llik <- l_likelihood(lhs.grid)
# Find MLE
v.mle.lhs <- lhs.grid[which.max(v.llik), ]
v.mle.lhs

#### 08 Find MLE via Directed-search ###
## Initial set of parameters to initialize Nelder-Mead
v.params.init <- c(p.S1S2 = 0.002, hr.S1 = 1, hr.S2 = 1)
## Run Directed Search
fit.nm.mle <- optim(v.params.init, l_likelihood, 
                     control = list(fnscale = -1, maxit = 1000),
                     lower = c(0.001, 0.1, 0.1),
                     upper = c(1, 10, 20),
                     method = "L-BFGS-B", hessian = T)
## Maximum Likelihood estimate
v.mle.nm <- fit.nm.mle$par
v.mle.nm
## Hessian matrix
m.hess <- fit.nm.mle$hessian
m.hess
## Covariance matrix
m.cov <- solve(-m.hess)
m.cov
## Correlation matrix
m.cor <- cov2cor(m.cov)
m.cor
## Standard error vectors
v.se <- sqrt(diag(m.cov))
v.se

#### 09 Plot model-predicted outputs at MLE ###
out.mle <- markov_sick_sicker(v.mle.lhs)
# real values
out.true <- markov_sick_sicker(v.params0)

par(mfrow = c(1, 3))
### Survival
plotrix::plotCI(x = 1:30, y = l.targets$Surv$value, 
                ui = l.targets$Surv$ub,
                li = l.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
lines(out.mle$Surv, col = "red")
points(out.mle$Surv, pch = 8, col = "red")
lines(out.true$Surv, col = "blue")

### Prevalence
plotrix::plotCI(x = 1:30, y = l.targets$Prev$value, 
                ui = l.targets$Prev$ub,
                li = l.targets$Prev$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Proportion Sick and Sicker")
lines(out.mle$Prev, col = "red")
points(out.mle$Prev, pch = 8, col = "red")
lines(out.true$Prev, col = "blue")

### Proportion Sick
plotrix::plotCI(x = c(10, 20, 30), y = l.targets$PropSick$value, 
                ui = l.targets$PropSick$ub,
                li = l.targets$PropSick$lb,
                ylim = c(0.24, 0.65),
                xlab = "Time", ylab = "Proportion Sick")
points(x = c(10, 20, 30), out.mle$PropSick, pch = 8, col = "red")
points(x = c(10, 20, 30), out.true$PropSick, pch = 8, col = "blue")
dev.off()

#### 10 Bayesian calibration ###
#### 10.1 Define Likelihood function ####
# Evaluates likelihood of a parameter set or sets
likelihood <- function(par_vector) { 
  exp(l_likelihood(par_vector)) 
}
# Test it
likelihood(v.params0) # It works!

#### 10.2 Define Log-Prior distribution ####
l_prior <- function(params){
  if(is.null(dim(params))) { # If vector, change to matrix
    params <- t(params) 
  }
  ## Number of parameters
  n.params <- ncol(params)
  ## Number of samples
  n.samp <- nrow(params)
  
  ### Create log-prior value
  lprior <- rep(0, n.samp)
  lprior <- lprior + dunif(params[, 1], 
                           min = 0.01, max = 0.5,
                           log = T) # p.S1S2
  lprior <- lprior + dunif(params[, 2], 
                           min = 1, max = 4.5,
                           log = T) # hr.S1
  lprior <- lprior + dunif(params[, 3], 
                           min = 5, max = 15,
                           log = T) # hr.S2
  return(lprior)
}
# Test it
l_prior(v.params0) # It works!

#### 10.3 Define Prior distribution ####
# Evaluates prior density of a parameter set or sets
prior <- function(params) { 
  exp(l_prior(params)) 
}
# Test it
prior(v.params0) # It works!

#### 10.4 Define Log-Posterior distribution ####
l_post <- function(params) { 
  log_post <- l_prior(params) + l_likelihood(params)
  return(log_post) 
}
# Test it
l_post(v.params0) # It works!

#### 10.5 MAP Estimate ####
## Run Nelder-Mead
fit.nm.map <- optim(v.params0, l_post, 
                     control = list(fnscale = -1, maxit = 1000),
                     hessian = T)
v.map <- fit.nm.map$par

#### 10.6 Run SIR ####
n.resamp <- 1000
### Compute likelihood
v.likelihood <- exp(v.llik)
### Estimate importance weights
weights <- v.likelihood/sum(v.likelihood)
### Sample from the posterior using importance weights
set.seed(31384)
post.index <- sample(x = n.samp, 
                     size = n.resamp, 
                     prob = weights,
                     replace = T)
post.sir <- lhs.grid[post.index, ]

### Compute posterior statistics
## Mean
mean.sir <- colMeans(post.sir)
mean.sir
## Median
median.sir <- colMedians(as.matrix(post.sir))
median.sir
## 95% Credible interval
quant.sir <- colQuantiles(as.matrix(post.sir), probs = c(0.025, 0.5, 0.975))
quant.sir

pairs.panels(post.sir)

#### 10.7 Run MCMC ####
fit.mcmc <- Metro_Hastings(li_func = l_post, 
                           pars = v.map,
                           par_names = c("p.S1S2", "hr.S1", "hr.S2"), 
                           iterations = 1e5L, # Default
                           burn_in = 1e4L)
## Chains
plot(fit.mcmc$trace[, 1], type = "l", ylim = range(fit.mcmc$trace))
lines(fit.mcmc$trace[, 2], col = "red")
lines(fit.mcmc$trace[, 3], col = "blue")

## Autocorrelation pre-thining
autocorr.plot(as.mcmc(fit.mcmc$trace))
ggs_autocorrelation(ggs(as.mcmc(fit.mcmc$trace)))

## Thining
fit.mcmc.thin <- mcmc_thin(fit.mcmc, thin = 20)
## Autocorrelation post-thining
autocorr.plot(as.mcmc(fit.mcmc.thin$trace))

## Caterpillar plots
ggs_caterpillar(ggs(as.mcmc(fit.mcmc.thin$trace)), horizontal = FALSE)

## Histogram
ggs_histogram(ggs(as.mcmc(fit.mcmc.thin$trace)), bins = 50)

ggs_pairs(ggs(as.mcmc(fit.mcmc.thin$trace)))

## Posterior estimates
post.mcmc.est <- BCI(fit.mcmc, interval = c(0.025, 0.5, 0.975))
## Tables of results
colnames(post.mcmc.est)[4] <- "Mean"
stargazer(post.mcmc.est, type = "text")

