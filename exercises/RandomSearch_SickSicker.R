
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




####################################################################
######  Load model as a function  ######
####################################################################
# - inputs are parameters to be estimated through calibration
# - outputs correspond to the target data




####################################################################
######  Specify calibration parameters  ######
####################################################################





####################################################################
######  Calibrate!  ######
####################################################################





####################################################################
######  Exploring best-fitting input sets  ######
####################################################################



