# Functions used in calibration

#######################################################################
##  gen_samp:  Generates a random sample of input parameter values   ##
#######################################################################

gen_samp = function(n,k,lb,ub,method="LHS"){
  # n: number of samples
  # k: number of parameters to be sampled
  # ub: upper bound for each parameter (vector)
  # lb: lower bound for each parameter (vector)
  # method: Random ("rand") or Latin-Hypercube ("LHS")
  
  m.param.samp = matrix(nrow=n,ncol=k)
  
  # check sampling method has been correctly specified
  if (!(method %in% c("LHS","rand"))){
    stop("Error: unrecognized sampling method. Set method equal to \"LHS\" or \"rand\".")
  }
  
  # generate unit (0-1) sample using specified method
  # Latin Hypercube
  if (method=="LHS"){
    m.lhs.unit <- randomLHS(n, k)
  }
  
  # random
  if (method=="rand"){ 
    m.lhs.unit = matrix(runif(n*k),nrow=n,ncol=k)
  }
    
  # rescale sample based on lower and upper bounds of each parameter
  for (i in 1:n.param){
    m.param.samp[,i] = qunif(m.lhs.unit[,i],
                             min = lb[i],
                             max = ub[i])
  }
  colnames(m.param.samp) = names(lb)
  
  # return random sample of parameter values
  return(m.param.samp)
  
}



