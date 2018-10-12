## ----message=FALSE-------------------------------------------------------
# calibration functionality
library(lhs)
library(IMIS)
library(matrixStats) # package used for sumamry statistics
library(modeest)     # package used to estimate mode

# visualization
library(plotrix)
library(psych)

## ------------------------------------------------------------------------
load("CRSTargets.RData")

# Plot the targets

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = CRS.targets$Surv$Time, y = CRS.targets$Surv$value, 
                ui = CRS.targets$Surv$ub,
                li = CRS.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")


## ------------------------------------------------------------------------
source("markov_crs.R") # creates the function markov_crs()

## ------------------------------------------------------------------------
# names and number of input parameters to be calibrated
param.names <- c("p.Mets","p.DieMets")
n.param <- length(param.names)

# range on input search space
lb <- c(p.Mets = 0.04, p.DieMets = 0.04) # lower bound
ub <- c(p.Mets = 0.16, p.DieMets = 0.16) # upper bound

sample.prior <- function(n.samp){
  m.lhs.unit   <- lhs::randomLHS(n = n.samp, k = n.param)
  m.param.samp <- matrix(nrow = n.samp, ncol = n.param)
  colnames(m.param.samp) <- param.names
  for (i in 1:n.param){
    m.param.samp[, i] <- qunif(m.lhs.unit[,i],
                               min = lb[i],
                               max = ub[i])
  }
  return(m.param.samp)
}

# view resulting parameter set samples
psych::pairs.panels(sample.prior(1000))

## ------------------------------------------------------------------------
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
  }
  return(lprior)
}

## ------------------------------------------------------------------------
v.params.test <- c("p.Mets" = 0.1, "p.DieMets" = 0.1)
f_log_prior(v.params = v.params.test)

## ------------------------------------------------------------------------
prior <- function(v.params) { 
  exp(f_log_prior(v.params)) 
}

## ------------------------------------------------------------------------

# number of calibration targets
target.names <- c("Surv")
n.target     <- length(target.names)

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
      ###   Run model for parameter set "v.params" ###
      model.res <- markov_crs(v.params[j, ])
      
      ###  Calculate log-likelihood of model outputs to targets  ###
      # Survival ("Surv")
      v.llik[j, 1] <- sum(dnorm(x = CRS.targets$Surv$value,
                                mean = model.res$Surv,
                                sd = CRS.targets$Surv$se,
                                log = T))
      
      # OVERALL 
      llik.overall[j] <- sum(v.llik[j, ])
    }, error = function(e) NA) 
    if (is.na(jj)) { llik.overall <- -Inf }
  } # End loop over sampled parameter sets
  # return LLIK
  return(llik.overall)
}

## ------------------------------------------------------------------------
f_llik(v.params = v.params.test)

## ------------------------------------------------------------------------
likelihood <- function(v.params){ 
  exp(f_llik(v.params)) 
}
likelihood(v.params = v.params.test)

## ------------------------------------------------------------------------
f_log_post <- function(v.params) { 
  lpost <- f_log_prior(v.params) + f_llik(v.params)
  return(lpost) 
}
f_log_post(v.params = v.params.test)

## ------------------------------------------------------------------------
# number of resamples
n.resamp <- 1000

# run IMIS
fit.imis <- IMIS(B = 1000, # the incremental sample size at each iteration of IMIS.
                 B.re = n.resamp, # the desired posterior sample size
                 number_k = 10, # the maximum number of iterations in IMIS.
                 D = 0)

## ------------------------------------------------------------------------
# obtain posterior
m.calib.post <- fit.imis$resample
head(m.calib.post)

## ------------------------------------------------------------------------
# Plot the 1000 draws from the posterior
plot(m.calib.post,
     xlim = c(lb[1], ub[1]), ylim = c(lb[2], ub[2]),
     xlab = param.names[1], ylab = param.names[2])

# Plot the 1000 draws from the posterior with marginal histograms
psych::pairs.panels(m.calib.post)

## ----warning=FALSE-------------------------------------------------------
# Compute posterior mean
v.calib.post.mean <- colMeans(m.calib.post)
v.calib.post.mean

# Compute posterior median and 95% credible interval
m.calib.post.95cr <- matrixStats::colQuantiles(m.calib.post, probs = c(0.025, 0.5, 0.975))
m.calib.post.95cr

# Compute posterior mode
v.calib.post.mode <- apply(m.calib.post, 2, function(x) as.numeric(modeest::mlv(x)[1]))
v.calib.post.mode

# compute maximum a posteriori
v.calib.like <- likelihood(m.calib.post)
v.calib.post.map <- m.calib.post[which.max(v.calib.like), ]

## ------------------------------------------------------------------------
v.out.post.mode <- markov_crs(v.calib.post.mode)

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = CRS.targets$Surv$Time, y = CRS.targets$Surv$value, 
                ui = CRS.targets$Surv$ub,
                li = CRS.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
grid()
for (i in 1:nrow(m.calib.post)){
  mod_output <- markov_crs(m.calib.post[i, ])
  lines(x = CRS.targets$Surv$Time, 
        y = mod_output$Surv,
        col = "darkorange",
        lwd = 0.1)
}
lines(x = CRS.targets$Surv$Time, 
       y = v.out.post.mode$Surv,
      col = "dodgerblue",
      lwd = 2)

## ------------------------------------------------------------------------
v.out.post.map <- markov_crs(v.calib.post.map)

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = CRS.targets$Surv$Time, y = CRS.targets$Surv$value, 
                ui = CRS.targets$Surv$ub,
                li = CRS.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
grid()
for (i in 1:nrow(m.calib.post)){
  mod_output <- markov_crs(m.calib.post[i, ])
  lines(x = CRS.targets$Surv$Time, 
        y = mod_output$Surv,
        col = "darkorange",
        lwd = 0.1)
}
lines(x = CRS.targets$Surv$Time, 
       y = v.out.post.map$Surv,
      col = "dodgerblue",
      lwd = 2)

