
###################  Calibration Specifications  ###################

# Model: Cancer Relative Survival (CRS) Markov Model
# Inputs to be calibrated: p.Mets, p.DieMets
# Targets: Surv

# Search method: Directed search using Nelder-mead
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
v.params.test = c(p.Mets = 0.10, p.DieMets = 0.05)
markov_crs(v.params.test) # It works!


####################################################################
######  Specify calibration parameters  ######
####################################################################
# Specify seed (for reproducible sequence of random numbers)
set.seed(072218)

# number of initial starting points
n.init = 100

# names and number of input parameters to be calibrated
param.names = c("p.Mets","p.DieMets")
n.param = length(param.names)

# range on input search space
lb = c(p.Mets = 0.04, p.DieMets = 0.04) # lower bound
ub = c(p.Mets = 0.16, p.DieMets = 0.16) # upper bound

# number of calibration targets
target.names = c("Surv")
n.target = length(target.names)


####################################################################
######  Calibrate!  ######
####################################################################


###  Sample multiple random starting values for Nelder-Mead  ###
v.params.init = matrix(nrow=n.init,ncol=n.param)
for (i in 1:n.param){
  v.params.init[,i] = runif(n.init,min=lb[i],max=ub[i])
}
colnames(v.params.init) = param.names


###  Write goodness-of-fit function to pass to Nelder-Mead algorithm  ###

f_gof = function(v.params){
  ###   Run model for parametr set "v.params" ###
  model.res = markov_crs(v.params)
  
  ###  Calculate goodness-of-fit of model outputs to targets  ###
  v.GOF = numeric(n.target)
  # TARGET 1: Survival ("Surv")
  # log likelihood  
  v.GOF[1] = sum(dnorm(x = CRS.targets$Surv$value,
                         mean = model.res$Surv,
                         sd = CRS.targets$Surv$se,
                         log = T))
  
  # TARGET 2: (if you had more...)
  # log likelihood
  # v.GOF[2] = sum(dnorm(x = CRS.targets$Target2$value,
  #                        mean = model.res$Target2,
  #                        sd = CRS.targets$Target2$se,
  #                        log = T))
  
  # OVERALL
  # can give different targets different weights
  v.weights = rep(1,n.target)
  # weighted sum
  GOF.overall = sum(v.GOF[1:n.target] * v.weights)
  
  # return GOF
  return(GOF.overall)
}


###  Run Nelder-Mead for each starting point  ###
m.calib.res = matrix(nrow = n.init, ncol = n.param+1)
colnames(m.calib.res) = c(param.names, "Overall_fit")
for (j in 1:n.init){
  
  fit.nm = optim(v.params.init[j,], f_gof, 
                 control = list(fnscale = -1, # fnscale = -1 switches from minimization to maximization
                                maxit = 1000), 
                 hessian = T)
  
  m.calib.res[j,] = c(fit.nm$par,fit.nm$value)
  
}


####################################################################
######  Exploring best-fitting input sets  ######
####################################################################

# Arrange parameter sets in order of fit
m.calib.res = m.calib.res[order(-m.calib.res[,"Overall_fit"]),]

# Examine the top 10 best-fitting sets
m.calib.res[1:10,]
# Plot the top 10 (top 10%)
plot(m.calib.res[1:10,1],m.calib.res[1:10,2],
     xlim=c(lb[1],ub[1]),ylim=c(lb[2],ub[2]),
     xlab = colnames(m.calib.res)[1],ylab = colnames(m.calib.res)[2])

###### COMMENTS: ######
## 1. Add the model predicted output at best-fitting set or sets
