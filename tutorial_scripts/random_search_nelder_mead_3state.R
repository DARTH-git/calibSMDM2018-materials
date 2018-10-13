## ----message=FALSE-------------------------------------------------------
# calibration functionality
library(lhs)

# visualization
library(plotrix)
library(psych)

## ------------------------------------------------------------------------
load("CRSTargets.RData")

## ------------------------------------------------------------------------
plotrix::plotCI(x = CRS.targets$Surv$Time, y = CRS.targets$Surv$value, 
                ui = CRS.targets$Surv$ub,
                li = CRS.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")

## ------------------------------------------------------------------------
source("markov_crs.R")

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
rs.n.samp <- 1000

## ------------------------------------------------------------------------
# lower bound
lb <- c(p.Mets = 0.04, p.DieMets = 0.04) 

# upper bound
ub <- c(p.Mets = 0.16, p.DieMets = 0.16)

## ------------------------------------------------------------------------
# Sample unit Latin Hypercube
m.lhs.unit <- lhs::randomLHS(rs.n.samp, n.param)
colnames(m.lhs.unit) <- param.names
head(m.lhs.unit)

## ----fig.width = 6-------------------------------------------------------
plot(x = m.lhs.unit[, param.names[1]],
     y = m.lhs.unit[, param.names[2]],
     xlab = param.names[1],
     ylab = param.names[2])

## ------------------------------------------------------------------------
# Rescale to min/max of each parameter
rs.param.samp <- matrix(nrow=rs.n.samp,ncol=n.param)
colnames(rs.param.samp) <- param.names
for (i in 1:n.param){
  rs.param.samp[,i] <- qunif(m.lhs.unit[,i],
                           min = lb[i],
                           max = ub[i])
}

## ------------------------------------------------------------------------
# view resulting parameter set samples
psych::pairs.panels(rs.param.samp)

## ------------------------------------------------------------------------
# initialize goodness-of-fit vector
rs.GOF <- rep(0, rs.n.samp)

## ------------------------------------------------------------------------
gof_norm_loglike <- function(target_mean, target_sd, model_output){
  sum(dnorm(x = target_mean,
            mean = model_output,
            sd = target_sd,
            log = TRUE))
}

## ------------------------------------------------------------------------
for (j in 1:rs.n.samp){
  
  ###  Run model for a given parameter set  ###
  model.res = markov_crs(v.params = rs.param.samp[j, ])
  
  ###  Calculate goodness-of-fit of model outputs to targets  ###
  # log likelihood of the model output given the targets
  rs.GOF[j] = gof_norm_loglike(model_output = model.res$Surv,
                               target_mean = CRS.targets$Surv$value,
                               target_sd = CRS.targets$Surv$se)
}

## ------------------------------------------------------------------------
# Arrange parameter sets in order of fit
rs.calib.res <- cbind(rs.param.samp, rs.GOF)
rs.calib.res <- rs.calib.res[order(-rs.calib.res[,"rs.GOF"]),]

# Examine the top 10 best-fitting sets
rs.calib.res[1:10,]

# Plot the top 100 (top 10%)
plot(x = rs.calib.res[1:100,1], y = rs.calib.res[1:100,2],
     xlim=c(lb[1],ub[1]),ylim=c(lb[2],ub[2]),
     xlab = colnames(rs.calib.res)[1],ylab = colnames(rs.calib.res)[2])

## ------------------------------------------------------------------------
rs_best_fit_params <- c(rs.calib.res[1, c("p.Mets",  "p.DieMets")])
rs_best_fit_model <- markov_crs(rs_best_fit_params)

plotrix::plotCI(x = CRS.targets$Surv$Time, y = CRS.targets$Surv$value, 
                ui = CRS.targets$Surv$ub,
                li = CRS.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
points(x = names(rs_best_fit_model$Surv), y = rs_best_fit_model$Surv,
       col = "green",
       pch = 20)
legend("topright", legend = c("Targets", "Random Search"),
       col = c("black", "green"),
       pch = c(1, 20))

## ------------------------------------------------------------------------
n.init <- 100
nm.params.init <- matrix(nrow=n.init,ncol=n.param)
set.seed(101)
for (i in 1:n.param){
  nm.params.init[,i] <- runif(n.init,min=lb[i],max=ub[i])
}
colnames(nm.params.init) <- param.names

head(nm.params.init)

## ------------------------------------------------------------------------
nm_objective = function(v.params){
  ###   Run model for parametr set "v.params" ###
  model.res = markov_crs(v.params)

  # log likelihood  
  v.GOF = gof_norm_loglike(target_mean = CRS.targets$Surv$value,
                              model_output = model.res$Surv,
                              target_sd = CRS.targets$Surv$se)
  return(v.GOF)
}

## ------------------------------------------------------------------------
nm.calib.res <- matrix(nrow = n.init, ncol = n.param+1)
colnames(nm.calib.res) <- c(param.names, "Overall_fit")
for (j in 1:n.init){
  
  fit.nm <- optim(nm.params.init[j,], nm_objective, 
                 control = list(fnscale = -1, # fnscale = -1 switches from minimization to maximization
                                maxit = 1000), 
                 hessian = T)
  
  nm.calib.res[j,] <- c(fit.nm$par,fit.nm$value)
  
}

## ------------------------------------------------------------------------
# Arrange parameter sets in order of fit
nm.calib.res <- nm.calib.res[order(-nm.calib.res[,"Overall_fit"]),]

# Examine the top 10 best-fitting sets
nm.calib.res[1:10,]

# Plot the top 10 (top 10%)
plot(nm.calib.res[1:10,1],nm.calib.res[1:10,2],
     xlim=c(lb[1],ub[1]),ylim=c(lb[2],ub[2]),
     xlab = colnames(nm.calib.res)[1],ylab = colnames(nm.calib.res)[2])

## ------------------------------------------------------------------------
nm_best_fit_params <- c(nm.calib.res[1, c("p.Mets",  "p.DieMets")])
nm_best_fit_model <- markov_crs(nm_best_fit_params)

plotrix::plotCI(x = CRS.targets$Surv$Time, y = CRS.targets$Surv$value, 
                ui = CRS.targets$Surv$ub,
                li = CRS.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
points(x = names(nm_best_fit_model$Surv), y = nm_best_fit_model$Surv,
       col = "red",
       pch = 4,
       cex = 1.2)
legend("topright", legend = c("Targets", "Nelder-Mead"),
       col = c("black", "red"),
       pch = c(1, 2, 4))

## ------------------------------------------------------------------------
plotrix::plotCI(x = CRS.targets$Surv$Time, y = CRS.targets$Surv$value, 
                ui = CRS.targets$Surv$ub,
                li = CRS.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
points(x = names(rs_best_fit_model$Surv), y = rs_best_fit_model$Surv,
       col = "green",
       pch = 2,
       cex = 1.2)
points(x = names(nm_best_fit_model$Surv), y = nm_best_fit_model$Surv,
       col = "red",
       pch = 4,
       cex = 1.2)

legend("topright", legend = c("Targets", "Random Search", "Nelder-Mead"),
       col = c("black", "green", "red"),
       pch = c(1, 2, 4))

## ----code=readLines('markov_crs.R')--------------------------------------

