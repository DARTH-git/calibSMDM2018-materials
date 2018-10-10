
###################  Calibration Specifications  ###################

# Model: Cancer Relative Survival (CRS) Markov Model
# Inputs to be calibrated: p.Mets, p.DieMets
# Targets: Surv

# Search method: Nelder-Mead
# Goodness-of-fit measure: Sum of Log-Likelihood

####################################################################


####################################################################
######  Load packages and function files  ####
####################################################################
# calibration functionality
library(lhs)

# visualization
library(plotrix)
library(psych)

# calibration functions
source("Calibration_Functions.R")


####################################################################
######  Load target data  ######
####################################################################
load("CRSTargets.RData")

# Plot the targets
plotrix::plotCI(x = 2:60, y = CRS.targets$Surv$value, 
                ui = CRS.targets$Surv$ub,
                li = CRS.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")


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
seed = 072218

# number of initial starting points
n.init = 100

# names and number of input parameters to be calibrated
param.names = c("p.Mets","p.DieMets")
n.param=length(param.names)

# range on input search space
lb <- c(p.Mets = 0.04, p.DieMets = 0.04) # lower bound
ub <- c(p.Mets = 0.16, p.DieMets = 0.16) # upper bound

# calibration targets
v.target.mean = CRS.targets$Surv$value
v.target.se = CRS.targets$Surv$se

####################################################################
######  Calibrate!  ######
####################################################################

# set seed
set.seed(seed)

###  Sample multiple random starting values for Nelder-Mead  ###
v.params.init = matrix(nrow=n.init,ncol=n.param)
for (i in 1:n.param){
  v.params.init[,i] = runif(n.init,min=lb[i],max=ub[i])
}
colnames(v.params.init) = param.names


###  Write goodness-of-fit function to pass to Nelder-Mead algorithm  ###

f_gof = function(v.params){
  # run model for parametr set "v.params"
  model.res = markov_crs(v.params)
  # extract model output to compare to targets
  v.res = model.res$Surv
  # calculate goodness-of-fit of this model output to the target
  # log likelihood
  GOF = sum(dnorm(x = v.target.mean,
                  mean = v.res,
                  sd = v.target.se,
                  log = T))
  # return GOF
  return(GOF)
}


###  Run Nelder-Mead for each starting point  ###
m.calib.res = matrix(nrow=n.init,ncol=n.param+1)
for (j in 1:n.init){

fit.nm = optim(v.params.init[j,], f_gof, 
               control = list(fnscale = -1, # fnscale = -1 switches from minimization to maximization
                              maxit = 1000), 
               hessian = T)

m.calib.res[j,] = c(fit.nm$par,fit.nm$value)

}


####################################################################
######   Exploring best-fitting input sets  ######
####################################################################


# The Nelder-Mead run that resulted in the best-fitting (max log-likelihood) set is the best
v.calib.best = m.calib.res[which.max[m.calib.res],]

print("And the best-fitting set is...")
print(v.calib.best)


