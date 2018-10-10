#### Sick-Sicker Markov model in a function ####
markov_sick_sicker <- function(v.params) {
  with(as.list(v.params), {
    ## Markov model parameters
    age     <- 25                   # age at baseline
    max.age <- 55                   # maximum age of follow up
    n.t  <- max.age - age           # time horizon, number of cycles
    v.n  <- c("H", "S1", "S2", "D") # the 4 states of the model: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
    n.s <- length(v.n)              # number of health states 
    
    # Transition probabilities and hazard ratios
    p.HD    = 0.005 # probability to die when healthy
    p.HS1   = 0.15  # probability to become sick when healthy
    p.S1H   = 0.5   # probability to become healthy when sick
    # p.S1S2  = 0.105         	# probability to become sicker when sick
    # hr.S1   = 3           	  # hazard ratio of death in sick vs healthy
    # hr.S2   = 10          	  # hazard ratio of death in sicker vs healthy
    # compute internal paramters as a function of external parameters
    r.HD    = - log(1 - p.HD) # rate of death in healthy
    r.S1D   = hr.S1 * r.HD 	  # rate of death in sick
    r.S2D   = hr.S2 * r.HD  	# rate of death in sicker
    p.S1D   = 1 - exp(-r.S1D) # probability to die in sick
    p.S2D   = 1 - exp(-r.S2D) # probability to die in sicker
    
    ####### INITIALIZATION ##########################################
    # create the cohort trace
    m.M <- matrix(NA, nrow = n.t + 1 , 
                  ncol = n.s,
                  dimnames = list(0:n.t, v.n))     # create Markov trace (n.t + 1 because R doesn't understand  Cycle 0)
    
    m.M[1, ] <- c(1, 0, 0, 0)                      # initialize Markov trace
    
    # create transition probability matrix for NO treatment
    m.P <- matrix(0,
                  nrow = n.s, 
                  ncol = n.s,
                  dimnames = list(v.n, v.n))
    # fill in the transition probability array
    ### From Healthy
    m.P["H", "H"]  <- 1 - (p.HS1 + p.HD)
    m.P["H", "S1"] <- p.HS1
    m.P["H", "D"]  <- p.HD
    ### From Sick
    m.P["S1", "H"]  <- p.S1H
    m.P["S1", "S1"] <- 1 - (p.S1H + p.S1S2 + p.S1D)
    m.P["S1", "S2"] <- p.S1S2
    m.P["S1", "D"]  <- p.S1D
    ### From Sicker
    m.P["S2", "S2"] <- 1 - p.S2D
    m.P["S2", "D"]  <- p.S2D
    ### From Dead
    m.P["D", "D"] <- 1
    
    # check rows add up to 1
    if (!isTRUE(all.equal(as.numeric(rowSums(m.P)), as.numeric(rep(1, n.s))))) {
      stop("This is not a valid transition Matrix")
    }
    
    ############# PROCESS ###########################################
    
    for (t in 1:n.t){                              # throughout the number of cycles
      m.M[t + 1, ] <- m.M[t, ] %*% m.P           # estimate the Markov trace for cycle the next cycle (t + 1)
    }
    
    ####### EPIDEMIOLOGICAL OUTPUT  ###########################################
    #### Overall Survival (OS) ####
    v.os <- 1 - m.M[, "D"]                # calculate the overall survival (OS) probability for no treatment
    
    #### Disease prevalence #####
    v.prev <- rowSums(m.M[, c("S1", "S2")])/v.os
    
    #### Proportion of sick in S1 state #####
    v.prop.S1 <- m.M[, "S1"] / v.prev
    
    ####### RETURN OUTPUT  ###########################################
    out <- list(Surv = v.os[-1],
                Prev = v.prev[-1],
                PropSick = v.prop.S1[c(11, 21, 31)])
    
    return(out)
  }
  )
}