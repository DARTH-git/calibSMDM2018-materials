
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

# number of random samples
n.samp = 1000

# names and number of input parameters to be calibrated
param.names = c("p.S1S2", "hr.S1", "hr.S2")
n.param = length(param.names)

# range on input search space
lb = c(p.S1S2 = 0.01, hr.S1 = 1.0, hr.S2 = 5) # lower bound
ub = c(p.S1S2 = 0.50, hr.S1 = 4.5, hr.S2 = 15) # upper bound

# number of calibration targets
target.names = c("Surv", "Prev", "PropSick")
n.target = length(target.names)


####################################################################
######  Calibrate!  ######
####################################################################


###  Generate a random sample of input values  ###

# Sample unit Latin Hypercube
m.lhs.unit = randomLHS(n.samp, n.param)

# Rescale to min/max of each parameter
m.param.samp = matrix(nrow=n.samp,ncol=n.param)
for (i in 1:n.param){
  m.param.samp[,i] = qunif(m.lhs.unit[,i],
                           min = lb[i],
                           max = ub[i])
}
colnames(m.param.samp) = param.names

# view resulting parameter set samples
pairs.panels(m.param.samp)


###  Run the model for each set of input values ###

# initialize goodness-of-fit vector
m.GOF = matrix(nrow=n.samp,ncol=n.target)
colnames(m.GOF) = paste0(target.names,"_fit")

# loop through sampled sets of input values
for (j in 1:n.samp){
  
  ###  Run model for a given parameter set  ###
  model.res = markov_sick_sicker(v.params = m.param.samp[j, ])
  
  
  ###  Calculate goodness-of-fit of model outputs to targets  ###

  # TARGET 1: Survival ("Surv")
  # log likelihood  
  m.GOF[j,1] = sum(dnorm(x = SickSicker.targets$Surv$value,
                       mean = model.res$Surv,
                       sd = SickSicker.targets$Surv$se,
                       log = T))
  
  # weighted sum of squared errors (alternative to log likelihood)
  # w = 1/(CRS.targets$Surv$se^2)
  # m.GOF[j,1] = -sum(w*(CRS.targets$Surv$value - v.res)^2)
  
  
  # TARGET 2: Prevalence ("Prev")
  # log likelihood
  m.GOF[j,2] = sum(dnorm(x = SickSicker.targets$Prev$value,
                         mean = model.res$Prev,
                         sd = SickSicker.targets$Prev$se,
                         log = T))
  
  # TARGET 3: Proprotion Sick+Sicker who are Sick
  # log likelihood
  m.GOF[j,3] = sum(dnorm(x = SickSicker.targets$PropSick$value,
                         mean = model.res$PropSick,
                         sd = SickSicker.targets$PropSick$se,
                         log = T))
  
  
} # End loop over sampled parameter sets


###  Combine fits to the different targets into single GOF  ###
# can give different targets different weights
v.weights = matrix(1,nrow=n.target,ncol=1)
# matrix multiplication to calculate weight sum of each GOF matrix row
v.GOF.overall = c(m.GOF%*%v.weights)
# Store in GOF matrix with column name "Overall"
m.GOF = cbind(m.GOF,Overall_fit=v.GOF.overall)


####################################################################
######  Exploring best-fitting input sets  ######
####################################################################

# Arrange parameter sets in order of fit
m.calib.res = cbind(m.param.samp,m.GOF)
m.calib.res = m.calib.res[order(-m.calib.res[,"Overall_fit"]),]

# Examine the top 10 best-fitting sets
m.calib.res[1:10,]
# Plot the top 100 (top 10%)
scatterplot3d(x = m.calib.res[1:100, 1],
              y = m.calib.res[1:100, 2],
              z = m.calib.res[1:100, 3],
              xlim = c(lb[1],ub[1]), ylim = c(lb[2],ub[2]), zlim = c(lb[3],ub[3]),
              xlab = param.names[1], ylab = param.names[2], zlab = param.names[3])
