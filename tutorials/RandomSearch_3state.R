
###################  Calibration Specifications  ###################

# Model: Cancer Relative Survival (CRS) Markov Model
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
n.samp = 1000

# names and number of input parameters to be calibrated
param.names = c("p.Mets","p.DieMets")
n.param = length(param.names)

# range on input search space
lb <- c(p.Mets = 0.04, p.DieMets = 0.04) # lower bound
ub <- c(p.Mets = 0.16, p.DieMets = 0.16) # upper bound

# number of calibration targets
target.names = c("Surv")
n.target = length(target.names)


####################################################################
######  Calibrate!  ######
####################################################################


###  Generate a random sample of input values  ###

# Sample unit Latin Hypercube
m.lhs.unit <- randomLHS(n.samp, n.param)

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
m.GOF = matrix(nrow = n.samp, ncol = n.target)
colnames(m.GOF) = paste0(target.names, "_fit")

# loop through sampled sets of input values
for (j in 1:n.samp){
  
  ###  Run model for a given parameter set  ###
  model.res = markov_crs(v.params = m.param.samp[j, ])
  
  
  ###  Calculate goodness-of-fit of model outputs to targets  ###

  # TARGET 1: Survival ("Surv")
  # log likelihood  
  m.GOF[j,1] = sum(dnorm(x = CRS.targets$Surv$value,
                       mean = model.res$Surv,
                       sd = CRS.targets$Surv$se,
                       log = T))
  
  # weighted sum of squared errors (alternative to log likelihood)
  # w = 1/(CRS.targets$Surv$se^2)
  # m.GOF[j,1] = -sum(w*(CRS.targets$Surv$value - v.res)^2)
  
  
  # TARGET 2: (if you had more...)
  # log likelihood
  # m.GOF[j,2] = sum(dnorm(x = CRS.targets$Target2$value,
  #                        mean = model.res$Target2,
  #                        sd = CRS.targets$Target2$se,
  #                        log = T))
  
  
} # End loop over sampled parameter sets


###  Combine fits to the different targets into single GOF  ###
# can give different targets different weights
v.weights = matrix(1, nrow = n.target, ncol = 1)
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
plot(m.calib.res[1:100,1],m.calib.res[1:100,2],
     xlim=c(lb[1],ub[1]),ylim=c(lb[2],ub[2]),
     xlab = colnames(m.calib.res)[1],ylab = colnames(m.calib.res)[2])

