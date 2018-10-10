#########     Calibration of three-state CRS Markov model       #####################
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
library(pbapply)
library(lhs)
library(scales)
library(psych)
library(MHadaptive)
library(matrixStats)
library(ggmcmc)
library(coda)
source("../Functions.R")

#### 02 Load targets ####
load("CRSTargets.RData")

#### 02.1 Plot targets ####
### Survival
plotrix::plotCI(x = 2:60, y = l.targets$Surv$value, 
                ui = l.targets$Surv$ub,
                li = l.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")

#### 03 Three-state CRS Markov model in a function ####
#### 03.1 Load Function ####
source("Markov_CRS - Function.R")
#### 03.2 Test model ####
v.params0 <- c(p.Mets = 0.10, p.DieMets = 0.05)
markov_crs(v.params0) # It works!

#### 04 Weighted Sum of square errors function ####
w_sse <- function(params){
  if(is.null(dim(params))) { # If vector, change to matrix
    params <- t(params) 
  }
  wsse <- rep(0, nrow(params))
  for(j in 1:nrow(params)) { # j=1
    jj <- tryCatch( { 
      # Obtain model outputs from a set of parameters
      res_j <- markov_crs(v.params = params[j, ])
      # Survival => Normal likelihood
      wsse[j] <- sum(((l.targets$Surv$value - res_j$Surv)/l.targets$Surv$se)^2)
    }, error = function(e) NA) 
    if(is.na(jj)) { wsse[j] <- Inf }
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
      res_j <- markov_crs(v.params = params[j, ])
      # Survival => Normal likelihood
      llik[j] <- sum(dnorm(x = l.targets$Surv$value, 
                           mean = res_j$Surv, 
                           sd = l.targets$Surv$se, 
                           log = T))
    }, error = function(e) NA) 
    if(is.na(jj)) { llik[j] <- -Inf }
  } 
  return(llik)
}
# Test it
l_likelihood(v.params0) # It works!

#### 05.1 Draw Log-Likelihood contour (possible because of ONLY two parameters!) ####
## Ceate grid for contour
xmin <- ymin <- 0.04
xmax <- ymax <- 0.16
lengthx <- 100
xdata <- seq(xmin,xmax, length.out = lengthx)
ydata <- seq(ymin,ymax, length.out = lengthx)
## Create combination of both parameters´ values (similar to TWSA!)
llk.cont <- expand.grid(p.Mets = xdata,
                        p.DieMets = ydata)
llk.cont$loglik <- pbapply(llk.cont, 1, l_likelihood)
llk.cont <- subset(llk.cont, is.finite(loglik))
summary(llk.cont$loglik)
zdata <- matrix(llk.cont$loglik, ncol=length(ydata))
## Levels for contour
levels.vec <- c(-1000, -400, -50, 66, 120, 140, 150, 155, 156)
contour(xdata,ydata,zdata,
        levels=levels.vec,
        # nlevels = 80,
        xlab='p.Mets',ylab='p.DieMets')

#### 06 Function to generate a Latin hypercube sample (LHS) ####
gen_lhs <- function(n.samp = 1000, seed = 072218){
  set.seed(seed)              # set a seed to be able to reproduce the same results
  ### Lower bounds
  lb <- c(p.Mets = 0.04, p.DieMets = 0.04)
  ### Upper bounds
  ub <- c(p.Mets = 0.16, p.DieMets = 0.16)
  
  ### Latin hypercube sampling
  lhs.samp <- randomLHS(n = n.samp, k = 2) # n: number of samples, k: number of parameters
  
  df.lhs <- data.frame(
    p.Mets = qunif(lhs.samp[, 1],
                   min = lb["p.Mets"],  # Parameters of desired distribution
                   max = ub["p.Mets"]), # Parameters of desired distribution
    p.DieMets = qunif(lhs.samp[, 2],
                      min = lb["p.DieMets"], # Parameters of desired distribution
                      max = ub["p.DieMets"]) # Parameters of desired distribution
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
#### 07.1 Show MLE of LHS in Likelihood Contour ###
contour(xdata,ydata,zdata,
        levels=levels.vec,
        # nlevels = 80,
        xlab='p.Mets',ylab='p.DieMets')
points(v.mle.lhs, pch = 2, col = "red", cex = 2)
legend("topright", c("LHS MLE"), col = c("red"), pch = c(2))

#### 08 Find MLE via Directed-search ###
## Initial sets of parameters to initialize Nelder-Mead
v.params.init  <- c(p.Mets = 0.07, p.DieMets = 0.07)
v.params.init2 <- lhs.grid[1, ] # 0.1470942 0.08993956
## Run Directed Search from user-defined initial parameters
fit.nm.mle <- optim(v.params.init, l_likelihood, 
                     control = list(fnscale = -1, maxit = 1000),
                     hessian = T)
## Run Directed Search from first set of LHS
fit.nm.mle2 <- optim(v.params.init2, l_likelihood, 
                     control = list(fnscale = -1, maxit = 1000),
                     hessian = T)

### Maximum Likelihood estimates
v.mle.nm <- fit.nm.mle$par
v.mle.nm 
v.mle.nm2 <- fit.nm.mle2$par
v.mle.nm2 
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
#### 08.3 Show MLE of NM in Likelihood Contour ###
contour(xdata,ydata,zdata,
        levels=levels.vec,
        # nlevels = 80,
        xlab='p.Mets',ylab='p.DieMets')
points(v.mle.lhs, pch = 2, col = "red", cex = 2)
points(v.mle.nm[1], v.mle.nm[2], pch = 8, col = "blue", cex = 2)
points(v.mle.nm2[1], v.mle.nm2[2], pch = 8, col = "blue", cex = 2)
legend("topright", c("LHS MLE", "NM MLE"), col = c("red", "blue"), pch = c(2, 8))

#### 09 Plot model-predicted outputs at MLE ###
out.mle <- markov_crs(v.mle.nm)
### Survival
plotrix::plotCI(x = 2:60, y = l.targets$Surv$value, 
                ui = l.targets$Surv$ub,
                li = l.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive",
                main = "Relative Survival Target")
lines(x = 2:60, out.mle$Surv, col = "red")
points(x = 2:60, out.mle$Surv, pch = 8, col = "red")

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
                           min = 0.04, max = 0.16,
                           log = T) # p.Mets
  lprior <- lprior + dunif(params[, 2], 
                           min = 0.04, max = 0.16,
                           log = T) # p.DieMets
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
fit.nm.map
v.map <- fit.nm.map$par
#### 10.5.1 Show MAP in Likelihood Contour ###
contour(xdata, ydata, zdata,
        levels = levels.vec,
        xlab='p.Mets',ylab='p.DieMets')
points(v.mle.lhs, pch = 2, col = "red", cex = 2)
points(v.mle.nm[1], v.mle.nm[2], pch = 8, col = "blue", cex = 2)
points(v.mle.nm2[1], v.mle.nm2[2], pch = 8, col = "blue", cex = 2)
points(v.map[1], v.map[2], pch = 1, col = "black", cex = 2)
legend("topright", c("LHS MLE", "NM MLE", "MAP"), col = c("red", "blue", "black"), pch = c(2, 8, 1))
#### 10.5.2 Plot model-predicted outputs at MAP ###
out.map <- markov_crs(v.map)
### Survival
plotrix::plotCI(x = 2:60, y = l.targets$Surv$value, 
                ui = l.targets$Surv$ub,
                li = l.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive",
                main = "Relative Survival Target")
lines(x = 2:60, out.map$Surv, col = "red")
points(x = 2:60, out.map$Surv, pch = 8, col = "red")

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
### Display posterior histograms
pairs.panels(post.sir)

#### 10.6.1 Show Posterior in Likelihood Contour ###
contour(xdata, ydata, zdata,
        levels = levels.vec,
        xlab='p.Mets',ylab='p.DieMets')
points(post.sir, cex=0.5)

#### 10.6.2 Plot model-predicted outputs at MAP ###
post.crs.sir <- matrix(unlist(pbapply(post.sir[1:1000, ], 1, markov_crs)), 
                       ncol = 59, byrow = T)
post.crs.sir
## Summarize posterior output
post.crs.sir.summ <- data_summary(melt(post.crs.sir), 
                                   varname = "value", 
                                   groupnames = "Var2")
plotrix::plotCI(x = 2:60, y = l.targets$Surv$value, 
                ui = l.targets$Surv$ub,
                li = l.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive",
                main = "Relative Survival Target")
apply(post.crs.sir[1:100, ], 1, lines, x = 2:60, col = "gray")
lines(2:60, post.crs.sir.summ$value, lty = 1, lwd = 2, col = "blue")
lines(2:60, post.crs.sir.summ$lb, lty = 2, lwd = 2, col = "blue")
lines(2:60, post.crs.sir.summ$ub, lty = 2, lwd = 2, col = "blue")
legend("topright", legend = c("Relative Survival (n = 200 95% CI)",
                              "Model output for one posterior parameter set",
                              "Posterior mean", 
                              "Posterior 95% Cred Int"),
       lty = c(1, 2, 1, 2), col= c("black", "gray", "blue", "blue"))

#### 10.7 Run MCMC ####
library(coda)
library(ggmcmc)
library(stargazer)

fit.mcmc <- Metro_Hastings(li_func = l_post, 
                           pars = v.map,
                           par_names = c("p.Mets", "p.DieMets"), 
                           iterations = 1e5L, # Default
                           burn_in = 1e4L)
### Plot posterior
plot(fit.mcmc$trace)
## Chains
plot(fit.mcmc$trace[, 1], type = "l", ylim = range(fit.mcmc$trace))
lines(fit.mcmc$trace[, 2], col = "red")

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
