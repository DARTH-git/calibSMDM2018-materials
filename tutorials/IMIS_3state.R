
###################  Calibration Specifications  ###################

# Model: Cancer Relative Survival (CRS) Markov Model
# Inputs to be calibrated: p.Mets, p.DieMets
# Targets: Surv

# Calibration method: Incremental mixture importance sampling (IMIS)
# Goodness-of-fit measure: Sum of Log-Likelihood

####################################################################


####################################################################
######  Load packages and function files  ####
####################################################################
# calibration functionality
library(lhs)
library(IMIS)
library(matrixStats) # package used for sumamry statistics
library(modeest)     # package used to estimate mode

# visualization
library(plotrix)
library(psych)


####################################################################
######  Load target data  ######
####################################################################
load("CRSTargets.RData")

# Plot the targets

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = CRS.targets$Surv$Time, y = CRS.targets$Surv$value, 
                ui = CRS.targets$Surv$ub,
                li = CRS.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")

# TARGET 2: (if you had more...)
# plotrix::plotCI(x = CRS.targets$Target2$Time, y = CRS.targets$Target2$value, 
#                 ui = CRS.targets$Target2$ub,
#                 li = CRS.targets$Target2$lb,
#                 ylim = c(0, 1), 
#                 xlab = "Time", ylab = "Target 2")


####################################################################
######  Load model as a function  ######
####################################################################
# - inputs are parameters to be estimated through calibration
# - outputs correspond to the target data

source("Markov_CRS - Function.R") # creates the function markov_crs()

# Check that it works
v.params.test <- c(p.Mets = 0.10, p.DieMets = 0.05)
markov_crs(v.params.test) # It works!


####################################################################
######  Specify calibration parameters  ######
####################################################################
# Specify seed (for reproducible sequence of random numbers)
set.seed(072218)

# number of random samples
n.resamp <- 1000

# names and number of input parameters to be calibrated
param.names <- c("p.Mets","p.DieMets")
n.param <- length(param.names)

# range on input search space
lb <- c(p.Mets = 0.04, p.DieMets = 0.04) # lower bound
ub <- c(p.Mets = 0.16, p.DieMets = 0.16) # upper bound

# number of calibration targets
target.names <- c("Surv")
n.target     <- length(target.names)


####################################################################
######  Calibrate!  ######
####################################################################


###  Write function to sample from prior ###

sample.prior <- function(n.samp){
  m.lhs.unit   <- randomLHS(n = n.samp, k = n.param)
  m.param.samp <- matrix(nrow = n.samp, ncol = n.param)
  colnames(m.param.samp) <- param.names
  for (i in 1:n.param){
    m.param.samp[, i] <- qunif(m.lhs.unit[,i],
                               min = lb[i],
                               max = ub[i])
    # ALTERNATIVE prior using beta (or other) distributions
    # m.param.samp[, i] <- qbeta(m.lhs.unit[,i],
    #                            min = 1,
    #                            max = 1)
  }
  return(m.param.samp)
}

# view resulting parameter set samples
pairs.panels(sample.prior(1000))

###  Write functions to evaluate log-prior and prior ###

f_log_prior <- function(v.params){
  if(is.null(dim(v.params))) { # If vector, change to matrix
    v.params <- t(v.params) 
  }
  n.samp <- nrow(v.params)
  colnames(v.params) <- param.names
  lprior <- rep(0, n.samp)
  for (i in 1:n.param){
    lprior <- lprior + dunif(v.params[, i],
                             min = lb[i],
                             max = ub[i], 
                             log = T)
    # ALTERNATIVE prior using beta distributions
    # lprior <- lprior + dbeta(v.params[, i],
    #                          min = 1,
    #                          max = 1, 
    #                          log = T)
  }
  return(lprior)
}
f_log_prior(v.params = v.params.test)

prior <- function(v.params) { 
  exp(f_log_prior(v.params)) 
}

###  Write log-likelihood and likelihood functions to pass to IMIS algorithm  ###

f_llik <- function(v.params){
  # par_vector: a vector (or matrix) of model parameters 
  if(is.null(dim(v.params))) { # If vector, change to matrix
    v.params <- t(v.params) 
  }
  n.samp <- nrow(v.params)
  v.llik <- matrix(0, nrow = n.samp, ncol = n.target) 
  llik.overall <- numeric(n.samp)
  for(j in 1:n.samp) { # j=1
    jj <- tryCatch( { 
      ###   Run model for parametr set "v.params" ###
      model.res <- markov_crs(v.params[j, ])
      
      ###  Calculate log-likelihood of model outputs to targets  ###
      # TARGET 1: Survival ("Surv")
      # log likelihood  
      v.llik[j, 1] <- sum(dnorm(x = CRS.targets$Surv$value,
                                mean = model.res$Surv,
                                sd = CRS.targets$Surv$se,
                                log = T))
      
      # TARGET 2: (if you had more...)
      # log likelihood
      # v.llik[j, 2] <- sum(dnorm(x = CRS.targets$Target2$value,
      #                        mean = model.res$Target2,
      #                        sd = CRS.targets$Target2$se,
      #                        log = T))
      
      # OVERALL 
      llik.overall[j] <- sum(v.llik[j, ])
    }, error = function(e) NA) 
    if(is.na(jj)) { llik.overall <- -Inf }
  } # End loop over sampled parameter sets
  # return LLIK
  return(llik.overall)
}
f_llik(v.params = v.params.test)

likelihood <- function(v.params){ 
  exp(f_llik(v.params)) 
}
likelihood(v.params = v.params.test)

###  Write function to evaluate log-posterior ###
f_log_post <- function(v.params) { 
  lpost <- f_log_prior(v.params) + f_llik(v.params)
  return(lpost) 
}
f_log_post(v.params = v.params.test)

###  Bayesian calibration using IMIS  ###

# number of resamples
n.resamp <- 1000

# run IMIS
fit.imis <- IMIS(B = 1000, # the incremental sample size at each iteration of IMIS.
                 B.re = n.resamp, # the desired posterior sample size
                 number_k = 10, # the maximum number of iterations in IMIS.
                 D = 0)

# obtain posterior
m.calib.post <- fit.imis$resample

####################################################################
######  Exploring posterior distribution  ######
####################################################################

# Plot the 1000 draws from the posterior
plot(m.calib.post,
     xlim = c(lb[1], ub[1]), ylim = c(lb[2], ub[2]),
     xlab = param.names[1], ylab = param.names[2])

# Plot the 1000 draws from the posterior with marginal histograms
pairs.panels(m.calib.post)

# Compute posterior mean
v.calib.post.mean <- colMeans(m.calib.post)

# Compute posterior median and 95% credible interval
m.calib.post.95cr <- colQuantiles(m.calib.post, probs = c(0.025, 0.5, 0.975))

# Compute posterior mode
v.calib.post.mode <- apply(m.calib.post, 2, function(x) as.numeric(mlv(x)[1]))


### Plot model-predicted output at mode vs targets ###
v.out.post.mode <- markov_crs(v.calib.post.mode)

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = CRS.targets$Surv$Time, y = CRS.targets$Surv$value, 
                ui = CRS.targets$Surv$ub,
                li = CRS.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
points(x = CRS.targets$Surv$Time, 
       y = v.out.post.mode$Surv, 
       pch = 8, col = "red")
legend("topright", 
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))

# TARGET 2: (if you had more...)
# plotrix::plotCI(x = CRS.targets$Target2$Time, y = CRS.targets$Target2$value, 
#                 ui = CRS.targets$Target2$ub,
#                 li = CRS.targets$Target2$lb,
#                 ylim = c(0, 1), 
#                 xlab = "Time", ylab = "Target 2")
# points(x = CRS.targets$Target2$Time, 
#        y = v.out.post.mode$Target2, 
#        pch = 8, col = "red")
# legend("topright", 
#        legend = c("Target", "Model-predicted output"),
#        col = c("black", "red"), pch = c(1, 8))