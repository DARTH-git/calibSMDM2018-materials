
###################  Calibration Specifications  ###################

# Model: Sick-Sicker 4-state Markov Model
# Inputs to be calibrated: p.S1S2, hr.S1, hr.S2
# Targets: Surv, Prev, PropSick

# Search method: Random search using Latin-Hypercube Sampling
# Goodness-of-fit measure: Sum of Log-Likelihood

####################################################################


####################################################################
######  Load packages and function files  ####
####################################################################
# calibration functionality
library(lhs)
library(IMIS)
library(modeest) # package used to estimate mode

# visualization
library(plotrix)
library(psych)
library(scatterplot3d) # now that we have three inputs to estimate, we'll need higher dimension visualization


####################################################################
######  Load target data  ######
####################################################################
load("SickSickerTargets.RData")

# Plot the targets

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = SickSicker.targets$Surv$Time, y = SickSicker.targets$Surv$value, 
                ui = SickSicker.targets$Surv$ub,
                li = SickSicker.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr(Alive)")

# TARGET 2: Prevalence ("Prev")
plotrix::plotCI(x = SickSicker.targets$Prev$Time, y = SickSicker.targets$Prev$value, 
                ui = SickSicker.targets$Prev$ub,
                li = SickSicker.targets$Prev$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr(Sick+Sicker)")

# TARGET 3: Proportion who are Sick ("PropSick"), among all those afflicted (Sick+Sicker)
plotrix::plotCI(x = SickSicker.targets$PropSick$Time, y = SickSicker.targets$PropSick$value, 
                ui = SickSicker.targets$PropSick$ub,
                li = SickSicker.targets$PropSick$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr(Sick | Sick+Sicker)")


####################################################################
######  Load model as a function  ######
####################################################################
# - inputs are parameters to be estimated through calibration
# - outputs correspond to the target data

source("Markov_Sick-Sicker - Function.R") # creates the function markov_sick_sicker()

# Check that it works
v.params.test = c(p.S1S2 = 0.105, hr.S1 = 3, hr.S2 = 10)
markov_sick_sicker(v.params0) # It works!


####################################################################
######  Specify calibration parameters  ######
####################################################################
# Specify seed (for reproducible sequence of random numbers)
set.seed(072218)

# number of random samples
n.resamp <- 1000

# names and number of input parameters to be calibrated
param.names <- c("p.S1S2", "hr.S1", "hr.S2")
n.param <- length(param.names)

# range on input search space
lb <- c(p.S1S2 = 0.01, hr.S1 = 1.0, hr.S2 = 5) # lower bound
ub <- c(p.S1S2 = 0.50, hr.S1 = 4.5, hr.S2 = 15) # upper bound

# number of calibration targets
target.names <- c("Surv", "Prev", "PropSick")
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
prior(v.params = v.params.test)

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
      model.res <- markov_sick_sicker(v.params[j, ])
      
      ###  Calculate log-likelihood of model outputs to targets  ###
      # TARGET 1: Survival ("Surv")
      # log likelihood  
      v.llik[j, 1] <- sum(dnorm(x = SickSicker.targets$Surv$value,
                                mean = model.res$Surv,
                                sd = SickSicker.targets$Surv$se,
                                log = T))
      
      # TARGET 2: Prevalence ("Prev")
      # log likelihood
      v.llik[j,2] = sum(dnorm(x = SickSicker.targets$Prev$value,
                             mean = model.res$Prev,
                             sd = SickSicker.targets$Prev$se,
                             log = T))
      
      # TARGET 3: Proportion Sick+Sicker who are Sick
      # log likelihood
      v.llik[j,3] = sum(dnorm(x = SickSicker.targets$PropSick$value,
                             mean = model.res$PropSick,
                             sd = SickSicker.targets$PropSick$se,
                             log = T))
      
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
scatterplot3d(x = m.calib.post[, 1],
              y = m.calib.post[, 2],
              z = m.calib.post[, 3],
              xlim = c(lb[1],ub[1]), ylim = c(lb[2],ub[2]), zlim = c(lb[3],ub[3]),
              xlab = param.names[1], ylab = param.names[2], zlab = param.names[3])

# Plot the 1000 draws from the posterior with marginal histograms
pairs.panels(m.calib.post)

# Compute posterior mean
v.calib.post.mean <- colMeans(m.calib.post)

# Compute posterior median and 95% credible interval
m.calib.post.95cr <- colQuantiles(m.calib.post, probs = c(0.025, 0.5, 0.975))

# Compute posterior mode
v.calib.post.mode <- apply(m.calib.post, 2, function(x) as.numeric(mlv(x)[1]))

### Plot model-predicted output at mode vs targets ###
v.out.post.mode <- markov_sick_sicker(v.calib.post.mode)

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = SickSicker.targets$Surv$Time, y = SickSicker.targets$Surv$value, 
                ui = SickSicker.targets$Surv$ub,
                li = SickSicker.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr(Alive)")
points(x = SickSicker.targets$Surv$Time, 
       y = v.out.post.mode$Surv, 
       pch = 8, col = "red")
legend("bottomright", 
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))

# TARGET 2: Prevalence ("Prev")
plotrix::plotCI(x = SickSicker.targets$Prev$Time, y = SickSicker.targets$Prev$value, 
                ui = SickSicker.targets$Prev$ub,
                li = SickSicker.targets$Prev$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr(Sick+Sicker)")
points(x = SickSicker.targets$Prev$Time, 
       y = v.out.post.mode$Prev, 
       pch = 8, col = "red")
legend("bottomright", 
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))

# TARGET 3: Proportion who are Sick ("PropSick"), among all those afflicted (Sick+Sicker)
plotrix::plotCI(x = SickSicker.targets$PropSick$Time, y = SickSicker.targets$PropSick$value, 
                ui = SickSicker.targets$PropSick$ub,
                li = SickSicker.targets$PropSick$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr(Sick | Sick+Sicker)")
points(x = SickSicker.targets$PropSick$Time, 
       y = v.out.post.mode$PropSick, 
       pch = 8, col = "red")
legend("bottomleft", 
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))
