
###################  Calibration Specifications  ###################

# Model: Sick-Sicker 4-state Markov Model
# Inputs to be calibrated: p.S1S2, hr.S1, hr.S2
# Targets: Surv, Prev, PropSick

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
v.params0 = c(p.S1S2 = 0.105, hr.S1 = 3, hr.S2 = 10)
markov_sick_sicker(v.params0) # It works!


####################################################################
######  Specify calibration parameters  ######
####################################################################
# Specify seed (for reproducible sequence of random numbers)
set.seed(072218)

# number of initial starting points
n.init = 100

# names and number of input parameters to be calibrated
param.names = c("p.S1S2", "hr.S1", "hr.S2")
n.param=length(param.names)

# range on input search space
lb = c(p.S1S2 = 0.01, hr.S1 = 1.0, hr.S2 = 5) # lower bound
ub = c(p.S1S2 = 0.50, hr.S1 = 4.5, hr.S2 = 15) # upper bound

# number of calibration targets
target.names = c("Surv", "Prev", "PropSick")
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
  model.res = markov_sick_sicker(v.params)
  
  ###  Calculate goodness-of-fit of model outputs to targets  ###
  v.GOF = numeric(n.target)
  # TARGET 1: Survival ("Surv")
  # log likelihood  
  v.GOF[1] = sum(dnorm(x = SickSicker.targets$Surv$value,
                       mean = model.res$Surv,
                       sd = SickSicker.targets$Surv$se,
                       log = T))
  
  # TARGET 2: Prevalence ("Prev")
  # log likelihood
  v.GOF[2] = sum(dnorm(x = SickSicker.targets$Prev$value,
                       mean = model.res$Prev,
                       sd = SickSicker.targets$Prev$se,
                       log = T))
  
  # TARGET 3: Proprotion Sick+Sicker who are Sick
  # log likelihood
  v.GOF[3] = sum(dnorm(x = SickSicker.targets$PropSick$value,
                       mean = model.res$PropSick,
                       sd = SickSicker.targets$PropSick$se,
                       log = T))
  
  # OVERALL
  # can give different targets different weights
  v.weights = rep(1,n.target)
  # weighted sum
  GOF.overall = sum(v.GOF[1:n.target] * v.weights)
  
  # return GOF
  return(GOF.overall)
}


###  Run Nelder-Mead for each starting point  ###
m.calib.res = matrix(nrow=n.init,ncol=n.param+1)
colnames(m.calib.res) = c(param.names,"Overall_fit")
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
scatterplot3d(x = m.calib.res[1:10, 1],
              y = m.calib.res[1:10, 2],
              z = m.calib.res[1:10, 3],
              xlim = c(lb[1],ub[1]), ylim = c(lb[2],ub[2]), zlim = c(lb[3],ub[3]),
              xlab = param.names[1], ylab = param.names[2], zlab = param.names[3])

