## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE)

## ------------------------------------------------------------------------
# calibration functionality
library(lhs)

# visualization
library(plotrix)
library(psych)

## ------------------------------------------------------------------------
load("data/CRSTargets.RData")

## ------------------------------------------------------------------------
plotrix::plotCI(x = CRS.targets$Surv$Time, y = CRS.targets$Surv$value, 
                ui = CRS.targets$Surv$ub,
                li = CRS.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")

## ------------------------------------------------------------------------
source("scripts/markov_crs.R") 

## ------------------------------------------------------------------------
v.params.test <- c(p.Mets = 0.1, p.DieMets = 0.2)
test_results <- markov_crs(v.params.test)

## ------------------------------------------------------------------------
str(test_results)
head(test_results$Surv)
plot(test_results$Surv)

## ------------------------------------------------------------------------
plotrix::plotCI(x = CRS.targets$Surv$Time, y = CRS.targets$Surv$value, 
                ui = CRS.targets$Surv$ub,
                li = CRS.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
points(test_results$Surv,
       col = "green", pch = 20)
legend("topright", legend = c("Targets", "Outputs"),
       col = c("black", "green"),
       pch = c(1, 20))

## ------------------------------------------------------------------------
set.seed(1010)

## ------------------------------------------------------------------------
param.names <- c("p.Mets","p.DieMets")
n.param <- length(param.names)
n.samp <- 1000

## ------------------------------------------------------------------------
# lower bound
lb <- c(p.Mets = 0.04, p.DieMets = 0.04) 

# upper bound
ub <- c(p.Mets = 0.16, p.DieMets = 0.16)

## ------------------------------------------------------------------------
# Sample unit Latin Hypercube
m.lhs.unit <- lhs::randomLHS(n.samp, n.param)
colnames(m.lhs.unit) <- param.names
head(m.lhs.unit)

## ----fig.width = 6-------------------------------------------------------
plot(x = m.lhs.unit[, param.names[1]],
     y = m.lhs.unit[, param.names[2]],
     xlab = param.names[1],
     ylab = param.names[2])

## ------------------------------------------------------------------------
# Rescale to min/max of each parameter
m.param.samp <- matrix(nrow=n.samp,ncol=n.param)
colnames(m.param.samp) <- param.names
for (i in 1:n.param){
  m.param.samp[,i] <- qunif(m.lhs.unit[,i],
                           min = lb[i],
                           max = ub[i])
}

## ------------------------------------------------------------------------
# view resulting parameter set samples
psych::pairs.panels(m.param.samp)

## ------------------------------------------------------------------------
# initialize goodness-of-fit vector
m.GOF <- rep(0, n.samp)

## ------------------------------------------------------------------------
gof_norm_loglike <- function(target_mean, target_sd, model_output){
  sum(dnorm(x = target_mean,
            mean = model_output,
            sd = target_sd,
            log = TRUE))
}

## ------------------------------------------------------------------------
for (j in 1:n.samp){
  
  ###  Run model for a given parameter set  ###
  model.res = markov_crs(v.params = m.param.samp[j, ])
  
  ###  Calculate goodness-of-fit of model outputs to targets  ###
  # log likelihood of the model output given the targets
  m.GOF[j] = gof_norm_loglike(model_output = model.res$Surv,
                              target_mean = CRS.targets$Surv$value,
                              target_sd = CRS.targets$Surv$se)
}

## ------------------------------------------------------------------------
# Arrange parameter sets in order of fit
m.calib.res <- cbind(m.param.samp,m.GOF)
m.calib.res <- m.calib.res[order(-m.calib.res[,"m.GOF"]),]

# Examine the top 10 best-fitting sets
m.calib.res[1:10,]

# Plot the top 100 (top 10%)
plot(x = m.calib.res[1:100,1], y = m.calib.res[1:100,2],
     xlim=c(lb[1],ub[1]),ylim=c(lb[2],ub[2]),
     xlab = colnames(m.calib.res)[1],ylab = colnames(m.calib.res)[2])

## ------------------------------------------------------------------------

best_fit_params <- c(m.calib.res[1, c("p.Mets",  "p.DieMets")])

best_fit_model <- markov_crs(best_fit_params)

plotrix::plotCI(x = CRS.targets$Surv$Time, y = CRS.targets$Surv$value, 
                ui = CRS.targets$Surv$ub,
                li = CRS.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")

points(x = names(best_fit_model$Surv), y = best_fit_model$Surv,
       col = "green",
       pch = 20)
legend("topright", legend = c("Targets", "Outputs"),
       col = c("black", "green"),
       pch = c(1, 20))

## ------------------------------------------------------------------------
gof_wsse <- function(model_output, target_mean, target_se){
    w = 1/(target_se^2)
    gof <-  (-1) * sum(w*(model_output - target_mean)^2) 
    return(gof)
}

