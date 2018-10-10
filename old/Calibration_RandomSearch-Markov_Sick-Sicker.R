
###################  Calibration Specifications  ###################

# Model: Sick-Sicker Markov Model
# Inputs to be calibrated: p.Mets, p.DieMets
# Targets: Surv

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

# calibration functions
source("Calibration_Functions.R")


####################################################################
######  Load target data  ######
####################################################################
load("SickSickerTargets.RData")

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

# number of random samples
n.samp = 10000

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

###  Generate a random sample of input values  ###

# use "gen_samp()" function 
m.param.samp = gen_samp(n=n.samp,k=n.param,lb=lb,ub=ub,method="LHS")

# view resulting lhs sample
pairs.panels(m.param.samp)


###  Run the model for each set of input values and calculate fit to targets  ###

# initialize goodness-of-fit vector
v.GOF = numeric(length=n.samp)*NA 

# loop through sampled sets of input values
for (j in 1:n.samp){
  
  # run model for parametr set m.param.samp[j,]
  model.res = markov_crs(v.params = m.param.samp[j, ])
  
  # extract model output corresponding to target values
  v.res = model.res$Surv
  
  # calculate goodness-of-fit of this model output to the target
  # log likelihood
  v.GOF[j] = sum(dnorm(x = v.target.mean,
                       mean = v.res,
                       sd = v.target.se,
                       log = T))
  
  # weighted sum of squared errors (alternative to log likelihood)
  #w = 1/(v.target.se^2)
  #v.GOF[j] = -sum(w*(v.target.mean - v.res)^2)
}

####################################################################
######  Identify best-fitting input sets  ######
####################################################################

# Arrange parameter sets in order of fit
m.calib.res = cbind(m.param.samp,GOF=v.GOF)
m.calib.res = m.calib.res[order(-m.calib.res[,"GOF"]),]

print("And the best-fitting set is...")
print(m.calib.res[which.max(m.calib.res[,"GOF"]),])


